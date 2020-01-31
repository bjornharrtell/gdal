/******************************************************************************
 *
 * Project:  CPL - Common Portability Library
 * Purpose:  Implement VSI large file api for lz4 files (.lz4).
 * Author:   Björn Harrtell <bjorn at wololo dot org>
 * Author:   Even Rouault, even.rouault at spatialys.com
 *
 ******************************************************************************
 * Copyright (c) 2020, Björn Harrtell <bjorn at wololo dot org>
 * Copyright (c) 2008-2014, Even Rouault <even dot rouault at spatialys.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/


#include "cpl_port.h"
#include "cpl_conv.h"
#include "cpl_vsi.h"

#include <cerrno>
#include <climits>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#if HAVE_FCNTL_H
#  include <fcntl.h>
#endif
#if HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif
#include <zlib.h>
#include <lz4.h>

#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "cpl_error.h"
#include "cpl_minizip_ioapi.h"
#include "cpl_minizip_unzip.h"
#include "cpl_multiproc.h"
#include "cpl_string.h"
#include "cpl_time.h"
#include "cpl_vsi_virtual.h"
#include "cpl_worker_thread_pool.h"

CPL_CVSID("$Id$")

constexpr int Z_BUFSIZE = 65536;  // Original size is 16384
constexpr int gz_magic[2] = {0x1f, 0x8b};  // gzip magic header

// gzip flag byte.
#define ASCII_FLAG   0x01  // bit 0 set: file probably ascii text
#define HEAD_CRC     0x02  // bit 1 set: header CRC present
#define EXTRA_FIELD  0x04  // bit 2 set: extra field present
#define ORIG_NAME    0x08  // bit 3 set: original file name present
#define COMMENT      0x10  // bit 4 set: file comment present
#define RESERVED     0xE0  // bits 5..7: reserved

#define ALLOC(size) malloc(size)
#define TRYFREE(p) {if (p) free(p);}

#define CPL_VSIL_LZ4_RETURN(ret)   \
        CPLError(CE_Failure, CPLE_AppDefined, \
                 "In file %s, at line %d, return %d", __FILE__, __LINE__, ret)

// #define ENABLE_DEBUG 1

/************************************************************************/
/* ==================================================================== */
/*                       VSILz4Handle                                  */
/* ==================================================================== */
/************************************************************************/

typedef struct
{
    vsi_l_offset  posInBaseHandle;
    z_stream      stream;
    uLong         crc;
    int           transparent;
    vsi_l_offset  in;
    vsi_l_offset  out;
} Lz4Snapshot;

class VSILz4Handle final : public VSIVirtualHandle
{
    VSIVirtualHandle* m_poBaseHandle = nullptr;
#ifdef DEBUG
    vsi_l_offset      m_offset = 0;
#endif
    vsi_l_offset      m_compressed_size = 0;
    vsi_l_offset      m_uncompressed_size = 0;
    vsi_l_offset      offsetEndCompressedData = 0;
    uLong             m_expected_crc = 0;
    char             *m_pszBaseFileName = nullptr; /* optional */
    bool              m_bWriteProperties = false;
    bool              m_bCanSaveInfo = false;

    /* Fields from gz_stream structure */
    z_stream stream;
    int      z_err = Z_OK;   /* error code for last stream operation */
    int      z_eof = 0;   /* set if end of input file (but not necessarily of the uncompressed stream ! "in" must be null too ) */
    Byte     *inbuf = nullptr;  /* input buffer */
    Byte     *outbuf = nullptr; /* output buffer */
    uLong    crc = 0;     /* crc32 of uncompressed data */
    int      m_transparent = 0; /* 1 if input file is not a .gz file */
    vsi_l_offset  startOff = 0;   /* startOff of compressed data in file (header skipped) */
    vsi_l_offset  in = 0;      /* bytes into deflate or inflate */
    vsi_l_offset  out = 0;     /* bytes out of deflate or inflate */
    vsi_l_offset  m_nLastReadOffset = 0;

    Lz4Snapshot* snapshots = nullptr;
    vsi_l_offset snapshot_byte_interval = 0; /* number of compressed bytes at which we create a "snapshot" */

    void check_header();
    int get_byte();
    int gzseek( vsi_l_offset nOffset, int nWhence );
    int gzrewind ();
    uLong getLong ();

    CPL_DISALLOW_COPY_ASSIGN(VSILz4Handle)

  public:

    VSILz4Handle( VSIVirtualHandle* poBaseHandle,
                   const char* pszBaseFileName,
                   vsi_l_offset offset = 0,
                   vsi_l_offset compressed_size = 0,
                   vsi_l_offset uncompressed_size = 0,
                   uLong expected_crc = 0,
                   int transparent = 0 );
    ~VSILz4Handle() override;

    bool              IsInitOK() const { return inbuf != nullptr; }

    int Seek( vsi_l_offset nOffset, int nWhence ) override;
    vsi_l_offset Tell() override;
    size_t Read( void *pBuffer, size_t nSize, size_t nMemb ) override;
    size_t Write( const void *pBuffer, size_t nSize, size_t nMemb ) override;
    int Eof() override;
    int Flush() override;
    int Close() override;

    VSILz4Handle*    Duplicate();
    bool              CloseBaseHandle();

    vsi_l_offset      GetLastReadOffset() { return m_nLastReadOffset; }
    const char*       GetBaseFileName() { return m_pszBaseFileName; }

    void              SetUncompressedSize( vsi_l_offset nUncompressedSize )
        { m_uncompressed_size = nUncompressedSize; }
    vsi_l_offset      GetUncompressedSize() { return m_uncompressed_size; }

    void              SaveInfo_unlocked();
    void              UnsetCanSaveInfo() { m_bCanSaveInfo = false; }
};

class VSILz4FilesystemHandler final : public VSIFilesystemHandler
{
    CPL_DISALLOW_COPY_ASSIGN(VSILz4FilesystemHandler)

    CPLMutex* hMutex = nullptr;
    VSILz4Handle* poHandleLastLz4File = nullptr;
    bool           m_bInSaveInfo = false;

public:
    VSILz4FilesystemHandler() = default;
    ~VSILz4FilesystemHandler() override;

    VSIVirtualHandle *Open( const char *pszFilename,
                            const char *pszAccess,
                            bool bSetError ) override;
    VSILz4Handle *OpenLz4ReadOnly( const char *pszFilename,
                                     const char *pszAccess );
    int Stat( const char *pszFilename, VSIStatBufL *pStatBuf,
              int nFlags ) override;
    int Unlink( const char *pszFilename ) override;
    int Rename( const char *oldpath, const char *newpath ) override;
    int Mkdir( const char *pszDirname, long nMode ) override;
    int Rmdir( const char *pszDirname ) override;
    char **ReadDirEx( const char *pszDirname, int nMaxFiles ) override;

    const char* GetOptions() override;

    void SaveInfo( VSILz4Handle* poHandle );
    void SaveInfo_unlocked( VSILz4Handle* poHandle );
};

/************************************************************************/
/*                            Duplicate()                               */
/************************************************************************/

VSILz4Handle* VSILz4Handle::Duplicate()
{
    CPLAssert (m_offset == 0);
    CPLAssert (m_compressed_size != 0);
    CPLAssert (m_pszBaseFileName != nullptr);

    VSIFilesystemHandler *poFSHandler =
        VSIFileManager::GetHandler( m_pszBaseFileName );

    VSIVirtualHandle* poNewBaseHandle =
        poFSHandler->Open( m_pszBaseFileName, "rb" );

    if( poNewBaseHandle == nullptr )
        return nullptr;

    VSILz4Handle* poHandle = new VSILz4Handle(poNewBaseHandle,
                                                m_pszBaseFileName,
                                                0,
                                                m_compressed_size,
                                                m_uncompressed_size);
    if( !(poHandle->IsInitOK()) )
    {
        delete poHandle;
        return nullptr;
    }

    poHandle->m_nLastReadOffset = m_nLastReadOffset;

    // Most important: duplicate the snapshots!

    for( unsigned int i=0;
         i < m_compressed_size / snapshot_byte_interval + 1;
         i++ )
    {
        if( snapshots[i].posInBaseHandle == 0 )
            break;

        poHandle->snapshots[i].posInBaseHandle = snapshots[i].posInBaseHandle;
        inflateCopy( &poHandle->snapshots[i].stream, &snapshots[i].stream);
        poHandle->snapshots[i].crc = snapshots[i].crc;
        poHandle->snapshots[i].transparent = snapshots[i].transparent;
        poHandle->snapshots[i].in = snapshots[i].in;
        poHandle->snapshots[i].out = snapshots[i].out;
    }

    return poHandle;
}

/************************************************************************/
/*                     CloseBaseHandle()                                */
/************************************************************************/

bool VSILz4Handle::CloseBaseHandle()
{
    bool bRet = true;
    if( m_poBaseHandle )
        bRet = VSIFCloseL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) == 0;
    m_poBaseHandle = nullptr;
    return bRet;
}

/************************************************************************/
/*                       VSILz4Handle()                                */
/************************************************************************/

VSILz4Handle::VSILz4Handle( VSIVirtualHandle* poBaseHandle,
                              const char* pszBaseFileName,
                              vsi_l_offset offset,
                              vsi_l_offset compressed_size,
                              vsi_l_offset uncompressed_size,
                              uLong expected_crc,
                              int transparent ) :
    m_poBaseHandle(poBaseHandle),
#ifdef DEBUG
    m_offset(offset),
#endif
    m_uncompressed_size(uncompressed_size),
    m_expected_crc(expected_crc),
    m_pszBaseFileName(pszBaseFileName ? CPLStrdup(pszBaseFileName) : nullptr),
    m_bWriteProperties(CPLTestBool(
        CPLGetConfigOption("CPL_VSIL_LZ4_WRITE_PROPERTIES", "YES"))),
    m_bCanSaveInfo(CPLTestBool(
        CPLGetConfigOption("CPL_VSIL_LZ4_SAVE_INFO", "YES"))),
    stream(),
    crc(0),
    m_transparent(transparent)
{
    if( compressed_size || transparent )
    {
        m_compressed_size = compressed_size;
    }
    else
    {
        if( VSIFSeekL(reinterpret_cast<VSILFILE*>(poBaseHandle), 0, SEEK_END) != 0 )
            CPLError(CE_Failure, CPLE_FileIO, "Seek() failed");
        m_compressed_size = VSIFTellL(reinterpret_cast<VSILFILE*>(poBaseHandle)) - offset;
        compressed_size = m_compressed_size;
    }
    offsetEndCompressedData = offset + compressed_size;

    if( VSIFSeekL(reinterpret_cast<VSILFILE*>(poBaseHandle), offset, SEEK_SET) != 0 )
        CPLError(CE_Failure, CPLE_FileIO, "Seek() failed");

    stream.zalloc = nullptr;
    stream.zfree = nullptr;
    stream.opaque = nullptr;
    stream.next_in = inbuf = nullptr;
    stream.next_out = outbuf = nullptr;
    stream.avail_in = stream.avail_out = 0;

    inbuf = static_cast<Byte *>(ALLOC(Z_BUFSIZE));
    stream.next_in = inbuf;

    int err = inflateInit2(&(stream), -MAX_WBITS);
    // windowBits is passed < 0 to tell that there is no zlib header.
    // Note that in this case inflate *requires* an extra "dummy" byte
    // after the compressed stream in order to complete decompression and
    // return Z_STREAM_END. Here the gzip CRC32 ensures that 4 bytes are
    // present after the compressed stream.
    if( err != Z_OK || inbuf == nullptr )
    {
        CPLError(CE_Failure, CPLE_NotSupported, "inflateInit2 init failed");
        TRYFREE(inbuf);
        inbuf = nullptr;
        return;
    }
    stream.avail_out = static_cast<uInt>(Z_BUFSIZE);

    if( offset == 0 ) check_header();  // Skip the .gz header.
    startOff = VSIFTellL(reinterpret_cast<VSILFILE*>(poBaseHandle)) - stream.avail_in;

    if( transparent == 0 )
    {
        snapshot_byte_interval = std::max(
            static_cast<vsi_l_offset>(Z_BUFSIZE), compressed_size / 100);
        snapshots = static_cast<Lz4Snapshot *>(
            CPLCalloc(sizeof(Lz4Snapshot),
                      static_cast<size_t>(
                          compressed_size / snapshot_byte_interval + 1)));
    }
}

/************************************************************************/
/*                      SaveInfo_unlocked()                             */
/************************************************************************/

void VSILz4Handle::SaveInfo_unlocked()
{
    if( m_pszBaseFileName && m_bCanSaveInfo )
    {
        VSIFilesystemHandler *poFSHandler =
            VSIFileManager::GetHandler( "/vsilz4/" );
        reinterpret_cast<VSILz4FilesystemHandler*>(poFSHandler)->
                                                    SaveInfo_unlocked(this);
        m_bCanSaveInfo = false;
    }
}

/************************************************************************/
/*                      ~VSILz4Handle()                                */
/************************************************************************/

VSILz4Handle::~VSILz4Handle()
{
    if( m_pszBaseFileName && m_bCanSaveInfo )
    {
        VSIFilesystemHandler *poFSHandler =
            VSIFileManager::GetHandler( "/vsilz4/" );
        reinterpret_cast<VSILz4FilesystemHandler*>(poFSHandler)->
            SaveInfo(this);
    }

    if( stream.state != nullptr )
    {
        inflateEnd(&(stream));
    }

    TRYFREE(inbuf);
    TRYFREE(outbuf);

    if( snapshots != nullptr )
    {
        for( size_t i=0;
             i < m_compressed_size / snapshot_byte_interval + 1;
             i++ )
        {
            if( snapshots[i].posInBaseHandle )
            {
                inflateEnd(&(snapshots[i].stream));
            }
        }
        CPLFree(snapshots);
    }
    CPLFree(m_pszBaseFileName);

    if( m_poBaseHandle )
        CPL_IGNORE_RET_VAL(VSIFCloseL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
}

/************************************************************************/
/*                      check_header()                                  */
/************************************************************************/

void VSILz4Handle::check_header()
{
    // Assure two bytes in the buffer so we can peek ahead -- handle case
    // where first byte of header is at the end of the buffer after the last
    // gzip segment.
    uInt len = stream.avail_in;
    if( len < 2 )
    {
        if( len ) inbuf[0] = stream.next_in[0];
        errno = 0;
        len = static_cast<uInt>(
            VSIFReadL(inbuf + len, 1, static_cast<size_t>(Z_BUFSIZE) >> len,
                      reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
#ifdef ENABLE_DEBUG
        CPLDebug("LZ4", CPL_FRMT_GUIB " " CPL_FRMT_GUIB,
                 VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)),
                 offsetEndCompressedData);
#endif
        if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) > offsetEndCompressedData )
        {
            len = len + static_cast<uInt>(
                offsetEndCompressedData - VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
            if( VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle),
                          offsetEndCompressedData, SEEK_SET) != 0 )
                z_err = Z_DATA_ERROR;
        }
        if( len == 0 )  // && ferror(file)
        {
            if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) !=
                offsetEndCompressedData )
                z_err = Z_ERRNO;
        }
        stream.avail_in += len;
        stream.next_in = inbuf;
        if( stream.avail_in < 2 )
        {
            m_transparent = stream.avail_in;
            return;
        }
    }

    // Peek ahead to check the gzip magic header.
    if( stream.next_in[0] != gz_magic[0] ||
        stream.next_in[1] != gz_magic[1]) {
        m_transparent = 1;
        return;
    }
    stream.avail_in -= 2;
    stream.next_in += 2;

    // Check the rest of the gzip header.
    const int method = get_byte();
    const int flags = get_byte();
    if( method != Z_DEFLATED || (flags & RESERVED) != 0 )
    {
        z_err = Z_DATA_ERROR;
        return;
    }

    // Discard time, xflags and OS code:
    for( len = 0; len < 6; len++ )
        CPL_IGNORE_RET_VAL(get_byte());

    if( (flags & EXTRA_FIELD) != 0 )
    {
        // Skip the extra field.
        len = static_cast<uInt>(get_byte()) & 0xFF;
        len += (static_cast<uInt>(get_byte()) & 0xFF) << 8;
        // len is garbage if EOF but the loop below will quit anyway.
        while( len != 0 && get_byte() != EOF )
        {
            --len;
        }
    }

    int c = 0;
    if( (flags & ORIG_NAME) != 0 )
    {
        // Skip the original file name.
        while( (c = get_byte()) != 0 && c != EOF ) {}
    }
    if( (flags & COMMENT) != 0 )
    {
        // skip the .gz file comment.
        while ((c = get_byte()) != 0 && c != EOF) {}
    }
    if( (flags & HEAD_CRC) != 0 )
    {
        // Skip the header crc.
        for( len = 0; len < 2; len++ )
            CPL_IGNORE_RET_VAL(get_byte());
    }
    z_err = z_eof ? Z_DATA_ERROR : Z_OK;
}

/************************************************************************/
/*                            get_byte()                                */
/************************************************************************/

int VSILz4Handle::get_byte()
{
    if( z_eof ) return EOF;
    if( stream.avail_in == 0 )
    {
        errno = 0;
        stream.avail_in = static_cast<uInt>(
            VSIFReadL(inbuf, 1, static_cast<size_t>(Z_BUFSIZE),
                      reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
#ifdef ENABLE_DEBUG
        CPLDebug("LZ4", CPL_FRMT_GUIB " " CPL_FRMT_GUIB,
                 VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)),
                 offsetEndCompressedData);
#endif
        if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) > offsetEndCompressedData )
        {
            stream.avail_in =
                stream.avail_in +
                static_cast<uInt>(
                    offsetEndCompressedData -
                    VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
            if( VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle),
                          offsetEndCompressedData, SEEK_SET) != 0 )
                return EOF;
        }
        if( stream.avail_in == 0 ) {
            z_eof = 1;
            if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) !=
                offsetEndCompressedData )
                z_err = Z_ERRNO;
            // if( ferror(file) ) z_err = Z_ERRNO;
            return EOF;
        }
        stream.next_in = inbuf;
    }
    stream.avail_in--;
    return *(stream.next_in)++;
}

/************************************************************************/
/*                            gzrewind()                                */
/************************************************************************/

int VSILz4Handle::gzrewind ()
{
    z_err = Z_OK;
    z_eof = 0;
    stream.avail_in = 0;
    stream.next_in = inbuf;
    crc = 0;
    if( !m_transparent )
        CPL_IGNORE_RET_VAL(inflateReset(&stream));
    in = 0;
    out = 0;
    return VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle), startOff, SEEK_SET);
}

/************************************************************************/
/*                              Seek()                                  */
/************************************************************************/

int VSILz4Handle::Seek( vsi_l_offset nOffset, int nWhence )
{
    /* The semantics of gzseek are different from ::Seek */
    /* It returns the current offset, where as ::Seek should return 0 */
    /* if successful */
    int ret = gzseek(nOffset, nWhence);
    return (ret >= 0) ? 0 : ret;
}

/************************************************************************/
/*                            gzseek()                                  */
/************************************************************************/

int VSILz4Handle::gzseek( vsi_l_offset offset, int whence )
{
    const vsi_l_offset original_offset = offset;
    const int original_nWhence = whence;

    z_eof = 0;
#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "Seek(" CPL_FRMT_GUIB ",%d)", offset, whence);
#endif

    if( m_transparent )
    {
        stream.avail_in = 0;
        stream.next_in = inbuf;
        if( whence == SEEK_CUR )
        {
            if( out + offset > m_compressed_size )
            {
                CPL_VSIL_LZ4_RETURN(-1);
                return -1L;
            }

            offset = startOff + out + offset;
        }
        else if( whence == SEEK_SET )
        {
            if( offset > m_compressed_size )
            {
                CPL_VSIL_LZ4_RETURN(-1);
                return -1L;
            }

            offset = startOff + offset;
        }
        else if( whence == SEEK_END )
        {
            // Commented test: because vsi_l_offset is unsigned (for the moment)
            // so no way to seek backward. See #1590 */
            if( offset > 0 ) // || -offset > compressed_size
            {
                CPL_VSIL_LZ4_RETURN(-1);
                return -1L;
            }

            offset = startOff + m_compressed_size - offset;
        }
        else
        {
            CPL_VSIL_LZ4_RETURN(-1);
            return -1L;
        }

        if( VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle), offset, SEEK_SET) < 0 )
        {
            CPL_VSIL_LZ4_RETURN(-1);
            return -1L;
        }

        out = offset - startOff;
        in = out;
#ifdef ENABLE_DEBUG
        CPLDebug("LZ4", "return " CPL_FRMT_GUIB, in);
#endif
        return in > INT_MAX ? INT_MAX : static_cast<int>(in);
    }

    // whence == SEEK_END is unsuppored in original gzseek.
    if( whence == SEEK_END )
    {
        // If we known the uncompressed size, we can fake a jump to
        // the end of the stream.
        if( offset == 0 && m_uncompressed_size != 0 )
        {
            out = m_uncompressed_size;
            return 1;
        }

        // We don't know the uncompressed size. This is unfortunate.
        // Do the slow version.
        static int firstWarning = 1;
        if( m_compressed_size > 10 * 1024 * 1024 && firstWarning )
        {
            CPLError(CE_Warning, CPLE_AppDefined,
                     "VSIFSeekL(xxx, SEEK_END) may be really slow "
                     "on Lz4 streams.");
            firstWarning = 0;
        }

        whence = SEEK_CUR;
        offset = 1024 * 1024 * 1024;
        offset *= 1024 * 1024;
    }

    // Rest of function is for reading only.

    // Compute absolute position.
    if( whence == SEEK_CUR )
    {
        offset += out;
    }

    // For a negative seek, rewind and use positive seek.
    if( offset >= out )
    {
        offset -= out;
    }
    else if( gzrewind() < 0 )
    {
            CPL_VSIL_LZ4_RETURN(-1);
            return -1L;
    }

    if( z_err != Z_OK && z_err != Z_STREAM_END )
    {
        CPL_VSIL_LZ4_RETURN(-1);
        return -1L;
    }

    for( unsigned int i = 0;
         i < m_compressed_size / snapshot_byte_interval + 1;
         i++ )
    {
        if( snapshots[i].posInBaseHandle == 0 )
            break;
        if( snapshots[i].out <= out + offset &&
            (i == m_compressed_size / snapshot_byte_interval ||
             snapshots[i+1].out == 0 || snapshots[i+1].out > out+offset) )
        {
            if( out >= snapshots[i].out )
                break;

#ifdef ENABLE_DEBUG
            CPLDebug(
                "SNAPSHOT", "using snapshot %d : "
                "posInBaseHandle(snapshot)=" CPL_FRMT_GUIB
                " in(snapshot)=" CPL_FRMT_GUIB
                " out(snapshot)=" CPL_FRMT_GUIB
                " out=" CPL_FRMT_GUIB
                " offset=" CPL_FRMT_GUIB,
                i, snapshots[i].posInBaseHandle, snapshots[i].in,
                snapshots[i].out, out, offset);
#endif
            offset = out + offset - snapshots[i].out;
            if( VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle),
                          snapshots[i].posInBaseHandle, SEEK_SET) != 0 )
                CPLError(CE_Failure, CPLE_FileIO, "Seek() failed");

            inflateEnd(&stream);
            inflateCopy(&stream, &snapshots[i].stream);
            crc = snapshots[i].crc;
            m_transparent = snapshots[i].transparent;
            in = snapshots[i].in;
            out = snapshots[i].out;
            break;
        }
    }

    // Offset is now the number of bytes to skip.

    if( offset != 0 && outbuf == nullptr )
    {
        outbuf = static_cast<Byte*>(ALLOC(Z_BUFSIZE));
        if( outbuf == nullptr )
        {
            CPL_VSIL_LZ4_RETURN(-1);
            return -1L;
        }
    }

    if( original_nWhence == SEEK_END && z_err == Z_STREAM_END )
    {
#ifdef ENABLE_DEBUG
        CPLDebug("LZ4", "gzseek return " CPL_FRMT_GUIB, out);
#endif
        return static_cast<int>(out);
    }

    while( offset > 0 )
    {
        int size = Z_BUFSIZE;
        if( offset < static_cast<vsi_l_offset>(Z_BUFSIZE) )
            size = static_cast<int>(offset);

        int read_size =
            static_cast<int>(Read(outbuf, 1, static_cast<uInt>(size)));
        if( read_size == 0 )
        {
            // CPL_VSIL_LZ4_RETURN(-1);
            return -1L;
        }
        if( original_nWhence == SEEK_END )
        {
            if( size != read_size )
            {
                z_err = Z_STREAM_END;
                break;
            }
        }
        offset -= read_size;
    }
#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "gzseek return " CPL_FRMT_GUIB, out);
#endif

    if( original_offset == 0 && original_nWhence == SEEK_END )
    {
        m_uncompressed_size = out;

        if( m_pszBaseFileName &&
            !STARTS_WITH_CI(m_pszBaseFileName, "/vsicurl/") &&
            m_bWriteProperties )
        {
            CPLString osCacheFilename (m_pszBaseFileName);
            osCacheFilename += ".properties";

            // Write a .properties file to avoid seeking next time.
            VSILFILE* fpCacheLength = VSIFOpenL(osCacheFilename.c_str(), "wb");
            if( fpCacheLength )
            {
                char szBuffer[32] = {};

                CPLPrintUIntBig(szBuffer, m_compressed_size, 31);
                char* pszFirstNonSpace = szBuffer;
                while( *pszFirstNonSpace == ' ' ) pszFirstNonSpace++;
                CPL_IGNORE_RET_VAL(
                    VSIFPrintfL(fpCacheLength,
                                "compressed_size=%s\n", pszFirstNonSpace));

                CPLPrintUIntBig(szBuffer, m_uncompressed_size, 31);
                pszFirstNonSpace = szBuffer;
                while( *pszFirstNonSpace == ' ' ) pszFirstNonSpace++;
                CPL_IGNORE_RET_VAL(
                    VSIFPrintfL(fpCacheLength,
                                "uncompressed_size=%s\n", pszFirstNonSpace));

                CPL_IGNORE_RET_VAL(VSIFCloseL(fpCacheLength));
            }
        }
    }

    return static_cast<int>(out);
}

/************************************************************************/
/*                              Tell()                                  */
/************************************************************************/

vsi_l_offset VSILz4Handle::Tell()
{
#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "Tell() = " CPL_FRMT_GUIB, out);
#endif
    return out;
}

/************************************************************************/
/*                              Read()                                  */
/************************************************************************/

size_t VSILz4Handle::Read( void * const buf, size_t const nSize,
                            size_t const nMemb )
{
#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "Read(%p, %d, %d)", buf,
             static_cast<int>(nSize),
             static_cast<int>(nMemb));
#endif

    if( (z_eof && in == 0) || z_err == Z_STREAM_END )
    {
        z_eof = 1;
        in = 0;
#ifdef ENABLE_DEBUG
        CPLDebug("LZ4", "Read: Eof");
#endif
        return 0;  /* EOF */
    }

    const unsigned len =
        static_cast<unsigned int>(nSize) * static_cast<unsigned int>(nMemb);
    Bytef *pStart = static_cast<Bytef*>(buf);  // Start off point for crc computation.
    // == stream.next_out but not forced far (for MSDOS).
    Byte *next_out = static_cast<Byte *>(buf);
    stream.next_out = static_cast<Bytef *>(buf);
    stream.avail_out = len;

    while( stream.avail_out != 0 )
    {
        if( m_transparent )
        {
            // Copy first the lookahead bytes:
            uInt nRead = 0;
            uInt n = stream.avail_in;
            if( n > stream.avail_out )
                n = stream.avail_out;
            if( n > 0 )
            {
                memcpy (stream.next_out, stream.next_in, n);
                next_out += n;
                stream.next_out = next_out;
                stream.next_in += n;
                stream.avail_out -= n;
                stream.avail_in -= n;
                nRead += n;
            }
            if( stream.avail_out > 0 )
            {
                const uInt nToRead = static_cast<uInt>(
                    std::min(m_compressed_size - (in + nRead),
                             static_cast<vsi_l_offset>(stream.avail_out)));
                uInt nReadFromFile = static_cast<uInt>(
                    VSIFReadL(next_out, 1, nToRead, reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
                stream.avail_out -= nReadFromFile;
                nRead += nReadFromFile;
            }
            in += nRead;
            out += nRead;
            if( nRead < len )
                z_eof = 1;
#ifdef ENABLE_DEBUG
            CPLDebug("LZ4", "Read return %d", static_cast<int>(nRead / nSize));
#endif
            return static_cast<int>(nRead) / nSize;
        }
        if( stream.avail_in == 0 && !z_eof )
        {
            vsi_l_offset posInBaseHandle =
                VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle));
            if( posInBaseHandle - startOff > m_compressed_size )
            {
                // If we reach here, file size has changed (because at
                // construction time startOff + m_compressed_size marked the
                // end of file).
                // We should probably have a better fix than that, by detecting
                // at open time that the saved snapshot is not valid and
                // discarding it.
                CPLError(CE_Failure, CPLE_AppDefined,
                         "File size of underlying /vsilz4/ file has changed");
                z_eof = 1;
                in = 0;
                CPL_VSIL_LZ4_RETURN(0);
                return 0;
            }
            Lz4Snapshot* snapshot =
                &snapshots[(posInBaseHandle - startOff) /
                           snapshot_byte_interval];
            if( snapshot->posInBaseHandle == 0 )
            {
                snapshot->crc =
                    crc32(crc, pStart,
                          static_cast<uInt>(stream.next_out - pStart));
#ifdef ENABLE_DEBUG
                CPLDebug("SNAPSHOT",
                         "creating snapshot %d : "
                         "posInBaseHandle=" CPL_FRMT_GUIB
                         " in=" CPL_FRMT_GUIB
                         " out=" CPL_FRMT_GUIB
                         " crc=%X",
                         static_cast<int>((posInBaseHandle - startOff) /
                                          snapshot_byte_interval),
                         posInBaseHandle, in, out,
                         static_cast<unsigned int>(snapshot->crc));
#endif
                snapshot->posInBaseHandle = posInBaseHandle;
                inflateCopy(&snapshot->stream, &stream);
                snapshot->transparent = m_transparent;
                snapshot->in = in;
                snapshot->out = out;

                if( out > m_nLastReadOffset )
                    m_nLastReadOffset = out;
            }

            errno = 0;
            stream.avail_in = static_cast<uInt>(
                VSIFReadL(inbuf, 1, Z_BUFSIZE, reinterpret_cast<VSILFILE*>(m_poBaseHandle)));
#ifdef ENABLE_DEBUG
            CPLDebug("LZ4", CPL_FRMT_GUIB " " CPL_FRMT_GUIB,
                     VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)),
                     offsetEndCompressedData);
#endif
            if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) > offsetEndCompressedData )
            {
#ifdef ENABLE_DEBUG
                CPLDebug("LZ4", "avail_in before = %d", stream.avail_in);
#endif
                stream.avail_in =
                    stream.avail_in -
                    static_cast<uInt>(VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) -
                                      offsetEndCompressedData);
                if( VSIFSeekL(reinterpret_cast<VSILFILE*>(m_poBaseHandle),
                              offsetEndCompressedData, SEEK_SET) != 0 )
                    CPLError(CE_Failure, CPLE_FileIO, "Seek() failed");
#ifdef ENABLE_DEBUG
                CPLDebug("LZ4", "avail_in after = %d", stream.avail_in);
#endif
            }
            if( stream.avail_in == 0 )
            {
                z_eof = 1;
                if( VSIFTellL(reinterpret_cast<VSILFILE*>(m_poBaseHandle)) !=
                    offsetEndCompressedData )
                {
                    z_err = Z_ERRNO;
                    break;
                }
            }
            stream.next_in = inbuf;
        }
        in += stream.avail_in;
        out += stream.avail_out;
        z_err = inflate(& (stream), Z_NO_FLUSH);
        in -= stream.avail_in;
        out -= stream.avail_out;

        if( z_err == Z_STREAM_END && m_compressed_size != 2 )
        {
            // Check CRC and original size.
            crc = crc32(crc, pStart,
                        static_cast<uInt>(stream.next_out - pStart));
            pStart = stream.next_out;
            if( m_expected_crc )
            {
#ifdef ENABLE_DEBUG
                CPLDebug(
                    "LZ4",
                    "Computed CRC = %X. Expected CRC = %X",
                    static_cast<unsigned int>(crc),
                    static_cast<unsigned int>(m_expected_crc));
#endif
            }
            if( m_expected_crc != 0 && m_expected_crc != crc )
            {
                CPLError(CE_Failure, CPLE_FileIO,
                         "CRC error. Got %X instead of %X",
                         static_cast<unsigned int>(crc),
                         static_cast<unsigned int>(m_expected_crc));
                z_err = Z_DATA_ERROR;
            }
            else if( m_expected_crc == 0 )
            {
                const uLong read_crc =
                    static_cast<unsigned long>(getLong());
                if( read_crc != crc )
                {
                    CPLError(CE_Failure, CPLE_FileIO,
                             "CRC error. Got %X instead of %X",
                             static_cast<unsigned int>(crc),
                             static_cast<unsigned int>(read_crc));
                    z_err = Z_DATA_ERROR;
                }
                else
                {
                    CPL_IGNORE_RET_VAL(getLong());
                    // The uncompressed length returned by above getlong() may
                    // be different from out in case of concatenated .gz files.
                    // Check for such files:
                    check_header();
                    if( z_err == Z_OK )
                    {
                        inflateReset(& (stream));
                        crc = 0;
                    }
                }
            }
        }
        if( z_err != Z_OK || z_eof )
            break;
    }
    crc = crc32(crc, pStart, static_cast<uInt>(stream.next_out - pStart));

    size_t ret = (len - stream.avail_out) / nSize;
    if( z_err != Z_OK && z_err != Z_STREAM_END )
    {
        z_eof = 1;
        in = 0;
        CPLError(CE_Failure, CPLE_AppDefined,
                 "In file %s, at line %d, decompression failed with "
                 "z_err = %d, return = %d",
                 __FILE__, __LINE__, z_err, static_cast<int>(ret));
    }

#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "Read return %d (z_err=%d, z_eof=%d)",
             static_cast<int>(ret), z_err, z_eof);
#endif
    return ret;
}

/************************************************************************/
/*                              getLong()                               */
/************************************************************************/

uLong VSILz4Handle::getLong ()
{
    uLong x = static_cast<uLong>(get_byte()) & 0xFF;

    x += (static_cast<uLong>(get_byte()) & 0xFF) << 8;
    x += (static_cast<uLong>(get_byte()) & 0xFF) << 16;
    const int c = get_byte();
    if( c == EOF )
    {
        z_err = Z_DATA_ERROR;
        return 0;
    }
    x += static_cast<uLong>(c) << 24;
    // coverity[overflow_sink]
    return x;
}

/************************************************************************/
/*                              Write()                                 */
/************************************************************************/

size_t VSILz4Handle::Write( const void * /* pBuffer */,
                             size_t /* nSize */,
                             size_t /* nMemb */ )
{
    CPLError(CE_Failure, CPLE_NotSupported,
             "VSIFWriteL is not supported on Lz4 streams");
    return 0;
}

/************************************************************************/
/*                               Eof()                                  */
/************************************************************************/

int VSILz4Handle::Eof()
{
#ifdef ENABLE_DEBUG
    CPLDebug("LZ4", "Eof()");
#endif
    return z_eof && in == 0;
}

/************************************************************************/
/*                              Flush()                                 */
/************************************************************************/

int VSILz4Handle::Flush()
{
    return 0;
}

/************************************************************************/
/*                              Close()                                 */
/************************************************************************/

int VSILz4Handle::Close()
{
    return 0;
}

/************************************************************************/
/* ==================================================================== */
/*                       VSILz4WriteHandleMT                           */
/* ==================================================================== */
/************************************************************************/

class VSILz4WriteHandleMT final : public VSIVirtualHandle
{
    CPL_DISALLOW_COPY_ASSIGN(VSILz4WriteHandleMT)

    VSIVirtualHandle*  poBaseHandle_ = nullptr;
    vsi_l_offset       nCurOffset_ = 0;
    uLong              nCRC_ = 0;
    int                nDeflateType_ = CPL_DEFLATE_TYPE_ZLIB;
    bool               bAutoCloseBaseHandle_ = false;
    int                nThreads_ = 0;
    std::unique_ptr<CPLWorkerThreadPool> poPool_{};
    std::list<std::string*> aposBuffers_{};
    std::string*       pCurBuffer_ = nullptr;
    std::mutex         sMutex_{};
    int                nSeqNumberGenerated_ = 0;
    int                nSeqNumberExpected_ = 0;
    int                nSeqNumberExpectedCRC_ = 0;
    size_t             nChunkSize_ = 0;
    bool               bHasErrored_ = false;

    struct Job
    {
        VSILz4WriteHandleMT *pParent_ = nullptr;
        std::string*       pBuffer_ = nullptr;
        int                nSeqNumber_ = 0;
        bool               bFinish_ = false;
        bool               bInCRCComputation_ = false;

        std::string        sCompressedData_{};
        uLong              nCRC_ = 0;
    };
    std::list<Job*>  apoFinishedJobs_{};
    std::list<Job*>  apoCRCFinishedJobs_{};
    std::list<Job*>  apoFreeJobs_{};

    static void DeflateCompress(void* inData);
    static void CRCCompute(void* inData);
    bool ProcessCompletedJobs();
    Job* GetJobObject();
#ifdef DEBUG_VERBOSE
    void DumpState();
#endif

  public:
    VSILz4WriteHandleMT( VSIVirtualHandle* poBaseHandle,
                        int nThreads,
                        int nDeflateType,
                        bool bAutoCloseBaseHandleIn );

    ~VSILz4WriteHandleMT() override;

    int Seek( vsi_l_offset nOffset, int nWhence ) override;
    vsi_l_offset Tell() override;
    size_t Read( void *pBuffer, size_t nSize, size_t nMemb ) override;
    size_t Write( const void *pBuffer, size_t nSize, size_t nMemb ) override;
    int Eof() override;
    int Flush() override;
    int Close() override;
};

/************************************************************************/
/*                        VSILz4WriteHandleMT()                        */
/************************************************************************/

VSILz4WriteHandleMT::VSILz4WriteHandleMT(  VSIVirtualHandle* poBaseHandle,
                        int nThreads,
                        int nDeflateType,
                        bool bAutoCloseBaseHandleIn ):
    poBaseHandle_(poBaseHandle),
    nDeflateType_(nDeflateType),
    bAutoCloseBaseHandle_(bAutoCloseBaseHandleIn),
    nThreads_(nThreads)
{
    const char* pszChunkSize = CPLGetConfigOption
        ("CPL_VSIL_DEFLATE_CHUNK_SIZE", "1024K");
    nChunkSize_ = static_cast<size_t>(atoi(pszChunkSize));
    if( strchr(pszChunkSize, 'K') )
        nChunkSize_ *= 1024;
    else if( strchr(pszChunkSize, 'M') )
        nChunkSize_ *= 1024 * 1024;
    nChunkSize_ = std::max(static_cast<size_t>(32 * 1024),
                    std::min(static_cast<size_t>(UINT_MAX), nChunkSize_));

    for( int i = 0; i < 1 + nThreads_; i++ )
        aposBuffers_.emplace_back( new std::string() );

    if( nDeflateType == CPL_DEFLATE_TYPE_ZLIB )
    {
        char header[11] = {};

        // Write a very simple .gz header:
        snprintf( header, sizeof(header),
                    "%c%c%c%c%c%c%c%c%c%c", gz_magic[0], gz_magic[1],
                Z_DEFLATED, 0 /*flags*/, 0, 0, 0, 0 /*time*/, 0 /*xflags*/,
                0x03 );
        poBaseHandle_->Write( header, 1, 10 );
    }
}

/************************************************************************/
/*                       ~VSILz4WriteHandleMT()                        */
/************************************************************************/

VSILz4WriteHandleMT::~VSILz4WriteHandleMT()

{
    VSILz4WriteHandleMT::Close();
    for( auto& psJob: apoFinishedJobs_ )
    {
        delete psJob->pBuffer_;
        delete psJob;
    }
    for( auto& psJob: apoCRCFinishedJobs_ )
    {
        delete psJob->pBuffer_;
        delete psJob;
    }
    for( auto& psJob: apoFreeJobs_ )
    {
        delete psJob->pBuffer_;
        delete psJob;
    }
    for( auto& pstr: aposBuffers_ )
    {
        delete pstr;
    }
    delete pCurBuffer_;
}

/************************************************************************/
/*                               Close()                                */
/************************************************************************/

int VSILz4WriteHandleMT::Close()

{
    if( !poBaseHandle_ )
        return 0;

    int nRet = 0;

    if( !pCurBuffer_ )
        pCurBuffer_ = new std::string();

    {
        auto psJob = GetJobObject();
        psJob->bFinish_ = true;
        psJob->pParent_ = this;
        psJob->pBuffer_ = pCurBuffer_;
        pCurBuffer_ = nullptr;
        psJob->nSeqNumber_ = nSeqNumberGenerated_;
        VSILz4WriteHandleMT::DeflateCompress( psJob );
    }

    if( poPool_ )
    {
        poPool_->WaitCompletion(0);
    }
    if( !ProcessCompletedJobs() )
    {
        nRet = -1;
    }
    else
    {
        CPLAssert(apoFinishedJobs_.empty());
        if( nDeflateType_ == CPL_DEFLATE_TYPE_ZLIB )
        {
            if( poPool_ )
            {
                poPool_->WaitCompletion(0);
            }
            ProcessCompletedJobs();
        }
        CPLAssert(apoCRCFinishedJobs_.empty());
    }

    if( nDeflateType_ == CPL_DEFLATE_TYPE_ZLIB )
    {
        const GUInt32 anTrailer[2] = {
            CPL_LSBWORD32(static_cast<GUInt32>(nCRC_)),
            CPL_LSBWORD32(static_cast<GUInt32>(nCurOffset_))
        };

        if( poBaseHandle_->Write( anTrailer, 1, 8 ) < 8 )
        {
            nRet = -1;
        }
    }

    if( bAutoCloseBaseHandle_ )
    {
        int nRetClose = poBaseHandle_->Close();
        if( nRet == 0 )
            nRet = nRetClose;

        delete poBaseHandle_;
    }
    poBaseHandle_ = nullptr;

    return nRet;
}

/************************************************************************/
/*                                Read()                                */
/************************************************************************/

size_t VSILz4WriteHandleMT::Read( void * /* pBuffer */,
                                 size_t /* nSize */,
                                 size_t /* nMemb */ )
{
    CPLError(CE_Failure, CPLE_NotSupported,
             "VSIFReadL is not supported on Lz4 write streams");
    return 0;
}

/************************************************************************/
/*                        DeflateCompress()                             */
/************************************************************************/

void VSILz4WriteHandleMT::DeflateCompress(void* inData)
{
    Job* psJob = static_cast<Job*>(inData);

    CPLAssert( psJob->pBuffer_);

    z_stream           sStream;
    sStream.zalloc = nullptr;
    sStream.zfree = nullptr;
    sStream.opaque = nullptr;

    sStream.avail_in = static_cast<uInt>(psJob->pBuffer_->size());
    sStream.next_in = reinterpret_cast<Bytef*>(&(*psJob->pBuffer_)[0]);

    int ret = deflateInit2( &sStream, Z_DEFAULT_COMPRESSION,
        Z_DEFLATED,
        (psJob->pParent_->nDeflateType_ == CPL_DEFLATE_TYPE_ZLIB) ?
            MAX_WBITS : -MAX_WBITS, 8,
        Z_DEFAULT_STRATEGY );
    CPLAssertAlwaysEval( ret == Z_OK );

    size_t nRealSize = 0;

    while( sStream.avail_in > 0 )
    {
        psJob->sCompressedData_.resize(nRealSize + Z_BUFSIZE);
        sStream.avail_out = static_cast<uInt>(Z_BUFSIZE);
        sStream.next_out = reinterpret_cast<Bytef*>(
            &psJob->sCompressedData_[0]) + nRealSize;

        const int zlibRet = deflate( &sStream, Z_NO_FLUSH );
        CPLAssertAlwaysEval( zlibRet == Z_OK );

        nRealSize += static_cast<uInt>(Z_BUFSIZE) - sStream.avail_out;
    }

    psJob->sCompressedData_.resize(nRealSize + Z_BUFSIZE);
    sStream.avail_out = static_cast<uInt>(Z_BUFSIZE);
    sStream.next_out = reinterpret_cast<Bytef*>(
        &psJob->sCompressedData_[0]) + nRealSize;

    // Do a Z_SYNC_FLUSH and Z_FULL_FLUSH, so as to have two markers when
    // independent as pigz 2.3.4 or later. The following 9 byte sequence will be
    // found: 0x00 0x00 0xff 0xff 0x00 0x00 0x00 0xff 0xff
    // Z_FULL_FLUSH only is sufficient, but it is not obvious if a
    // 0x00 0x00 0xff 0xff marker in the codestream is just a SYNC_FLUSH (
    // without dictionary reset) or a FULL_FLUSH (with dictionary reset)
    {
        const int zlibRet = deflate( &sStream, Z_SYNC_FLUSH );
        CPLAssertAlwaysEval( zlibRet == Z_OK );
    }

    {
        const int zlibRet = deflate( &sStream, Z_FULL_FLUSH );
        CPLAssertAlwaysEval( zlibRet == Z_OK );
    }

    if( psJob->bFinish_ )
    {
        const int zlibRet = deflate( &sStream, Z_FINISH );
        CPLAssertAlwaysEval( zlibRet == Z_STREAM_END );
    }

    nRealSize += static_cast<uInt>(Z_BUFSIZE) - sStream.avail_out;
    psJob->sCompressedData_.resize(nRealSize);

    deflateEnd( &sStream );

    {
        std::lock_guard<std::mutex> oLock(psJob->pParent_->sMutex_);
        psJob->pParent_->apoFinishedJobs_.push_back(psJob);
    }
}

/************************************************************************/
/*                          CRCCompute()                                */
/************************************************************************/

void VSILz4WriteHandleMT::CRCCompute(void* inData)
{
    Job* psJob = static_cast<Job*>(inData);
    psJob->bInCRCComputation_ = true;
    psJob->nCRC_ = crc32(0U,
        reinterpret_cast<const Bytef*>(psJob->pBuffer_->data()),
        static_cast<uInt>(psJob->pBuffer_->size()));

    {
        std::lock_guard<std::mutex> oLock(psJob->pParent_->sMutex_);
        psJob->pParent_->apoCRCFinishedJobs_.push_back(psJob);
    }
}

/************************************************************************/
/*                                DumpState()                           */
/************************************************************************/

#ifdef DEBUG_VERBOSE
void VSILz4WriteHandleMT::DumpState()
{
    fprintf(stderr, "Finished jobs (expected = %d):\n", nSeqNumberExpected_); // ok
    for(const auto* psJob: apoFinishedJobs_ )
    {
        fprintf(stderr,  "seq number=%d, bInCRCComputation = %d\n",  // ok
                psJob->nSeqNumber_, psJob->bInCRCComputation_ ? 1 : 0);
    }
    fprintf(stderr, "Finished CRC jobs (expected = %d):\n",  // ok
            nSeqNumberExpectedCRC_);
    for(const auto* psJob: apoFinishedJobs_ )
    {
        fprintf(stderr,  "seq number=%d\n",  // ok
                psJob->nSeqNumber_);
    }
    fprintf(stderr, "apoFreeJobs_.size() = %d\n",  // ok
            static_cast<int>(apoFreeJobs_.size()));
    fprintf(stderr, "aposBuffers_.size() = %d\n",  // ok
            static_cast<int>(aposBuffers_.size()));
}
#endif

/************************************************************************/
/*                         ProcessCompletedJobs()                       */
/************************************************************************/

bool VSILz4WriteHandleMT::ProcessCompletedJobs()
{
    std::lock_guard<std::mutex> oLock(sMutex_);
    bool do_it_again = true;
    while( do_it_again )
    {
        do_it_again = false;
        if( nDeflateType_ == CPL_DEFLATE_TYPE_ZLIB )
        {
            for( auto iter = apoFinishedJobs_.begin();
                    iter != apoFinishedJobs_.end(); ++iter )
            {
                auto psJob = *iter;

                if( !psJob->bInCRCComputation_ )
                {
                    psJob->bInCRCComputation_ = true;
                    sMutex_.unlock();
                    if( poPool_ )
                    {
                        poPool_->SubmitJob( VSILz4WriteHandleMT::CRCCompute,
                                            psJob );
                    }
                    else
                    {
                        CRCCompute(psJob);
                    }
                    sMutex_.lock();
                }
            }
        }

        for( auto iter = apoFinishedJobs_.begin();
                iter != apoFinishedJobs_.end(); ++iter )
        {
            auto psJob = *iter;
            if( psJob->nSeqNumber_ == nSeqNumberExpected_ )
            {
                apoFinishedJobs_.erase(iter);

                sMutex_.unlock();

                const size_t nToWrite = psJob->sCompressedData_.size();
                bool bError =
                    poBaseHandle_->Write( psJob->sCompressedData_.data(), 1,
                                          nToWrite) < nToWrite;
                sMutex_.lock();
                nSeqNumberExpected_ ++;

                if( nDeflateType_ != CPL_DEFLATE_TYPE_ZLIB )
                {
                    aposBuffers_.push_back(psJob->pBuffer_);
                    psJob->pBuffer_ = nullptr;

                    apoFreeJobs_.push_back(psJob);
                }

                if( bError )
                {
                    return false;
                }

                do_it_again = true;
                break;
            }
        }

        if( nDeflateType_ == CPL_DEFLATE_TYPE_ZLIB )
        {
            for( auto iter = apoCRCFinishedJobs_.begin();
                    iter != apoCRCFinishedJobs_.end(); ++iter )
            {
                auto psJob = *iter;
                if( psJob->nSeqNumber_ == nSeqNumberExpectedCRC_ )
                {
                    apoCRCFinishedJobs_.erase(iter);

                    nCRC_ = crc32_combine(nCRC_, psJob->nCRC_,
                                    static_cast<uLong>(psJob->pBuffer_->size()));

                    nSeqNumberExpectedCRC_ ++;

                    aposBuffers_.push_back(psJob->pBuffer_);
                    psJob->pBuffer_ = nullptr;

                    apoFreeJobs_.push_back(psJob);
                    do_it_again = true;
                    break;
                }
            }
        }
    }
    return true;
}

/************************************************************************/
/*                           GetJobObject()                             */
/************************************************************************/

VSILz4WriteHandleMT::Job* VSILz4WriteHandleMT::GetJobObject()
{
    {
        std::lock_guard<std::mutex> oLock(sMutex_);
        if( !apoFreeJobs_.empty() )
        {
            auto job = apoFreeJobs_.back();
            apoFreeJobs_.pop_back();
            job->sCompressedData_.clear();
            job->bInCRCComputation_ = false;
            return job;
        }
    }
    return new Job();
}

/************************************************************************/
/*                               Write()                                */
/************************************************************************/

size_t VSILz4WriteHandleMT::Write( const void * const pBuffer,
                                    size_t const nSize, size_t const nMemb )

{
    if( bHasErrored_ )
        return 0;

    const char* pszBuffer = static_cast<const char*>(pBuffer);
    size_t nBytesToWrite = nSize * nMemb;
    while( nBytesToWrite > 0 )
    {
        if( pCurBuffer_ == nullptr )
        {
            while(true)
            {
                {
                    std::lock_guard<std::mutex> oLock(sMutex_);
                    if( !aposBuffers_.empty() )
                    {
                        pCurBuffer_ = aposBuffers_.back();
                        aposBuffers_.pop_back();
                        break;
                    }
                }
                if( poPool_ )
                {
                    poPool_->WaitEvent();
                }
                if( !ProcessCompletedJobs() )
                {
                    bHasErrored_ = true;
                    return 0;
                }
            }
            pCurBuffer_->clear();
        }
        size_t nConsumed = std::min( nBytesToWrite,
                                     nChunkSize_ - pCurBuffer_->size() );
        pCurBuffer_->append(pszBuffer, nConsumed);
        nCurOffset_ += nConsumed;
        pszBuffer += nConsumed;
        nBytesToWrite -= nConsumed;
        if( pCurBuffer_->size() == nChunkSize_ )
        {
            if( poPool_ == nullptr )
            {
                poPool_.reset(new CPLWorkerThreadPool());
                if( !poPool_->Setup(nThreads_, nullptr, nullptr, false) )
                {
                    bHasErrored_ = true;
                    poPool_.reset();
                    return 0;
                }
            }

            auto psJob = GetJobObject();
            psJob->pParent_ = this;
            psJob->pBuffer_ = pCurBuffer_;
            psJob->nSeqNumber_ = nSeqNumberGenerated_;
            nSeqNumberGenerated_ ++;
            pCurBuffer_ = nullptr;
            poPool_->SubmitJob( VSILz4WriteHandleMT::DeflateCompress, psJob );
        }
    }

    return nMemb;
}

/************************************************************************/
/*                               Flush()                                */
/************************************************************************/

int VSILz4WriteHandleMT::Flush()

{
    // we *could* do something for this but for now we choose not to.

    return 0;
}

/************************************************************************/
/*                                Eof()                                 */
/************************************************************************/

int VSILz4WriteHandleMT::Eof()

{
    return 1;
}

/************************************************************************/
/*                                Seek()                                */
/************************************************************************/

int VSILz4WriteHandleMT::Seek( vsi_l_offset nOffset, int nWhence )

{
    if( nOffset == 0 && (nWhence == SEEK_END || nWhence == SEEK_CUR) )
        return 0;
    else if( nWhence == SEEK_SET && nOffset == nCurOffset_ )
        return 0;
    else
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Seeking on writable compressed data streams not supported.");

        return -1;
    }
}

/************************************************************************/
/*                                Tell()                                */
/************************************************************************/

vsi_l_offset VSILz4WriteHandleMT::Tell()

{
    return nCurOffset_;
}

/************************************************************************/
/* ==================================================================== */
/*                       VSILz4WriteHandle                             */
/* ==================================================================== */
/************************************************************************/

class VSILz4WriteHandle final : public VSIVirtualHandle
{
    CPL_DISALLOW_COPY_ASSIGN(VSILz4WriteHandle)

    VSIVirtualHandle*  m_poBaseHandle = nullptr;
    z_stream           sStream;
    Byte              *pabyInBuf = nullptr;
    Byte              *pabyOutBuf = nullptr;
    bool               bCompressActive = false;
    vsi_l_offset       nCurOffset = 0;
    uLong              nCRC = 0;
    int                nDeflateType = CPL_DEFLATE_TYPE_ZLIB;
    bool               bAutoCloseBaseHandle = false;

  public:
    VSILz4WriteHandle( VSIVirtualHandle* poBaseHandle, int nDeflateType,
                        bool bAutoCloseBaseHandleIn );

    ~VSILz4WriteHandle() override;

    int Seek( vsi_l_offset nOffset, int nWhence ) override;
    vsi_l_offset Tell() override;
    size_t Read( void *pBuffer, size_t nSize, size_t nMemb ) override;
    size_t Write( const void *pBuffer, size_t nSize, size_t nMemb ) override;
    int Eof() override;
    int Flush() override;
    int Close() override;
};

/************************************************************************/
/*                         VSILz4WriteHandle()                         */
/************************************************************************/

VSILz4WriteHandle::VSILz4WriteHandle( VSIVirtualHandle *poBaseHandle,
                                        int nDeflateTypeIn,
                                        bool bAutoCloseBaseHandleIn ) :
    m_poBaseHandle(poBaseHandle),
    sStream(),
    pabyInBuf(static_cast<Byte *>(CPLMalloc( Z_BUFSIZE ))),
    pabyOutBuf(static_cast<Byte *>(CPLMalloc( Z_BUFSIZE ))),
    nCRC(crc32(0L, nullptr, 0)),
    nDeflateType(nDeflateTypeIn),
    bAutoCloseBaseHandle(bAutoCloseBaseHandleIn)
{
    sStream.zalloc = nullptr;
    sStream.zfree = nullptr;
    sStream.opaque = nullptr;
    sStream.next_in = nullptr;
    sStream.next_out = nullptr;
    sStream.avail_in = sStream.avail_out = 0;

    sStream.next_in = pabyInBuf;

    if( deflateInit2( &sStream, Z_DEFAULT_COMPRESSION,
                      Z_DEFLATED,
                      ( nDeflateType == CPL_DEFLATE_TYPE_ZLIB) ?
                        MAX_WBITS : -MAX_WBITS, 8,
                      Z_DEFAULT_STRATEGY ) != Z_OK )
    {
        bCompressActive = false;
    }
    else
    {
        if( nDeflateType == CPL_DEFLATE_TYPE_ZLIB )
        {
            char header[11] = {};

            // Write a very simple .gz header:
            snprintf( header, sizeof(header),
                      "%c%c%c%c%c%c%c%c%c%c", gz_magic[0], gz_magic[1],
                    Z_DEFLATED, 0 /*flags*/, 0, 0, 0, 0 /*time*/, 0 /*xflags*/,
                    0x03 );
            m_poBaseHandle->Write( header, 1, 10 );
        }

        bCompressActive = true;
    }
}

/************************************************************************/
/*                       VSICreateLz4Writable()                        */
/************************************************************************/

VSIVirtualHandle* VSICreateLz4Writable( VSIVirtualHandle* poBaseHandle,
                                         int nDeflateTypeIn,
                                         int bAutoCloseBaseHandle )
{
    const char* pszThreads = CPLGetConfigOption("GDAL_NUM_THREADS", nullptr);
    if( pszThreads )
    {
        int nThreads = 0;
        if( EQUAL(pszThreads, "ALL_CPUS") )
            nThreads = CPLGetNumCPUs();
        else
            nThreads = atoi(pszThreads);
        nThreads = std::max(1, std::min(128, nThreads));
        if( nThreads > 1 )
        {
            // coverity[tainted_data]
            return new VSILz4WriteHandleMT( poBaseHandle,
                                                nThreads,
                                                nDeflateTypeIn,
                                                CPL_TO_BOOL(bAutoCloseBaseHandle) );
        }
    }
    return new VSILz4WriteHandle( poBaseHandle,
                                   nDeflateTypeIn,
                                   CPL_TO_BOOL(bAutoCloseBaseHandle) );
}

/************************************************************************/
/*                        ~VSILz4WriteHandle()                         */
/************************************************************************/

VSILz4WriteHandle::~VSILz4WriteHandle()

{
    if( bCompressActive )
        VSILz4WriteHandle::Close();

    CPLFree( pabyInBuf );
    CPLFree( pabyOutBuf );
}

/************************************************************************/
/*                               Close()                                */
/************************************************************************/

int VSILz4WriteHandle::Close()

{
    int nRet = 0;
    if( bCompressActive )
    {
        sStream.next_out = pabyOutBuf;
        sStream.avail_out = static_cast<uInt>(Z_BUFSIZE);

        const int zlibRet = deflate( &sStream, Z_FINISH );
        CPLAssertAlwaysEval( zlibRet == Z_STREAM_END );

        const size_t nOutBytes =
            static_cast<uInt>(Z_BUFSIZE) - sStream.avail_out;

        if( m_poBaseHandle->Write( pabyOutBuf, 1, nOutBytes ) < nOutBytes )
            return EOF;

        deflateEnd( &sStream );

        if( nDeflateType == CPL_DEFLATE_TYPE_ZLIB )
        {
            const GUInt32 anTrailer[2] = {
                CPL_LSBWORD32(static_cast<GUInt32>(nCRC)),
                CPL_LSBWORD32(static_cast<GUInt32>(nCurOffset))
            };

            m_poBaseHandle->Write( anTrailer, 1, 8 );
        }

        if( bAutoCloseBaseHandle )
        {
            nRet = m_poBaseHandle->Close();

            delete m_poBaseHandle;
        }

        bCompressActive = false;
    }

    return nRet;
}

/************************************************************************/
/*                                Read()                                */
/************************************************************************/

size_t VSILz4WriteHandle::Read( void * /* pBuffer */,
                                 size_t /* nSize */,
                                 size_t /* nMemb */ )
{
    CPLError(CE_Failure, CPLE_NotSupported,
             "VSIFReadL is not supported on Lz4 write streams");
    return 0;
}

/************************************************************************/
/*                               Write()                                */
/************************************************************************/

size_t VSILz4WriteHandle::Write( const void * const pBuffer,
                                  size_t const nSize, size_t const nMemb )

{
    size_t nBytesToWrite = nSize * nMemb;

    {
        size_t nOffset = 0;
        while( nOffset < nBytesToWrite )
        {
            uInt nChunk = static_cast<uInt>(
                std::min(static_cast<size_t>(UINT_MAX),
                         nBytesToWrite - nOffset));
            nCRC = crc32(nCRC,
                         reinterpret_cast<const Bytef *>(pBuffer) + nOffset,
                         nChunk);
            nOffset += nChunk;
        }
    }

    if( !bCompressActive )
        return 0;

    size_t nNextByte = 0;
    while( nNextByte < nBytesToWrite )
    {
        sStream.next_out = pabyOutBuf;
        sStream.avail_out = static_cast<uInt>(Z_BUFSIZE);

        if( sStream.avail_in > 0 )
            memmove( pabyInBuf, sStream.next_in, sStream.avail_in );

        const uInt nNewBytesToWrite = static_cast<uInt>(std::min(
            static_cast<size_t>(Z_BUFSIZE-sStream.avail_in),
            nBytesToWrite - nNextByte));
        memcpy( pabyInBuf + sStream.avail_in,
                reinterpret_cast<const Byte *>(pBuffer) + nNextByte,
                nNewBytesToWrite );

        sStream.next_in = pabyInBuf;
        sStream.avail_in += nNewBytesToWrite;

        const int zlibRet = deflate( &sStream, Z_NO_FLUSH );
        CPLAssertAlwaysEval( zlibRet == Z_OK );

        const size_t nOutBytes =
            static_cast<uInt>(Z_BUFSIZE) - sStream.avail_out;

        if( nOutBytes > 0 )
        {
            if( m_poBaseHandle->Write( pabyOutBuf, 1, nOutBytes ) < nOutBytes )
                return 0;
        }

        nNextByte += nNewBytesToWrite;
        nCurOffset += nNewBytesToWrite;
    }

    return nMemb;
}

/************************************************************************/
/*                               Flush()                                */
/************************************************************************/

int VSILz4WriteHandle::Flush()

{
    // we *could* do something for this but for now we choose not to.

    return 0;
}

/************************************************************************/
/*                                Eof()                                 */
/************************************************************************/

int VSILz4WriteHandle::Eof()

{
    return 1;
}

/************************************************************************/
/*                                Seek()                                */
/************************************************************************/

int VSILz4WriteHandle::Seek( vsi_l_offset nOffset, int nWhence )

{
    if( nOffset == 0 && (nWhence == SEEK_END || nWhence == SEEK_CUR) )
        return 0;
    else if( nWhence == SEEK_SET && nOffset == nCurOffset )
        return 0;
    else
    {
        CPLError(CE_Failure, CPLE_NotSupported,
                 "Seeking on writable compressed data streams not supported.");

        return -1;
    }
}

/************************************************************************/
/*                                Tell()                                */
/************************************************************************/

vsi_l_offset VSILz4WriteHandle::Tell()

{
    return nCurOffset;
}

/************************************************************************/
/* ==================================================================== */
/*                       VSILz4FilesystemHandler                       */
/* ==================================================================== */
/************************************************************************/

/************************************************************************/
/*                  ~VSILz4FilesystemHandler()                         */
/************************************************************************/

VSILz4FilesystemHandler::~VSILz4FilesystemHandler()
{
    if( poHandleLastLz4File )
    {
        poHandleLastLz4File->UnsetCanSaveInfo();
        delete poHandleLastLz4File;
    }

    if( hMutex != nullptr )
        CPLDestroyMutex( hMutex );
    hMutex = nullptr;
}

/************************************************************************/
/*                            SaveInfo()                                */
/************************************************************************/

void VSILz4FilesystemHandler::SaveInfo( VSILz4Handle* poHandle )
{
    CPLMutexHolder oHolder(&hMutex);
    SaveInfo_unlocked(poHandle);
}

void VSILz4FilesystemHandler::SaveInfo_unlocked( VSILz4Handle* poHandle )
{
    if( m_bInSaveInfo )
        return;
    m_bInSaveInfo = true;

    CPLAssert( poHandle != poHandleLastLz4File );
    CPLAssert(poHandle->GetBaseFileName() != nullptr);

    if( poHandleLastLz4File == nullptr ||
        strcmp(poHandleLastLz4File->GetBaseFileName(),
               poHandle->GetBaseFileName()) != 0 ||
        poHandle->GetLastReadOffset() >
            poHandleLastLz4File->GetLastReadOffset() )
    {
        VSILz4Handle* poTmp = poHandleLastLz4File;
        poHandleLastLz4File = nullptr;
        if( poTmp )
        {
            poTmp->UnsetCanSaveInfo();
            delete poTmp;
        }
        CPLAssert(poHandleLastLz4File == nullptr);
        poHandleLastLz4File = poHandle->Duplicate();
        if( poHandleLastLz4File )
            poHandleLastLz4File->CloseBaseHandle();
    }
    m_bInSaveInfo = false;
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

VSIVirtualHandle* VSILz4FilesystemHandler::Open( const char *pszFilename,
                                                  const char *pszAccess,
                                                  bool /* bSetError */ )
{
    if( !STARTS_WITH_CI(pszFilename, "/vsilz4/") )
        return nullptr;

    VSIFilesystemHandler *poFSHandler =
        VSIFileManager::GetHandler( pszFilename + strlen("/vsilz4/"));

/* -------------------------------------------------------------------- */
/*      Is this an attempt to write a new file without update (w+)      */
/*      access?  If so, create a writable handle for the underlying     */
/*      filename.                                                       */
/* -------------------------------------------------------------------- */
    if( strchr(pszAccess, 'w') != nullptr )
    {
        if( strchr(pszAccess, '+') != nullptr )
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Write+update (w+) not supported for /vsilz4, "
                     "only read-only or write-only.");
            return nullptr;
        }

        VSIVirtualHandle* poVirtualHandle =
            poFSHandler->Open( pszFilename + strlen("/vsilz4/"), "wb" );

        if( poVirtualHandle == nullptr )
            return nullptr;

        return VSICreateLz4Writable( poVirtualHandle,
                                       strchr(pszAccess, 'z') != nullptr,
                                       TRUE );
    }

/* -------------------------------------------------------------------- */
/*      Otherwise we are in the read access case.                       */
/* -------------------------------------------------------------------- */

    VSILz4Handle* poLZ4Handle = OpenLz4ReadOnly(pszFilename, pszAccess);
    if( poLZ4Handle )
        // Wrap the VSILz4Handle inside a buffered reader that will
        // improve dramatically performance when doing small backward
        // seeks.
        return VSICreateBufferedReaderHandle(poLZ4Handle);

    return nullptr;
}

/************************************************************************/
/*                          OpenLz4ReadOnly()                          */
/************************************************************************/

VSILz4Handle* VSILz4FilesystemHandler::OpenLz4ReadOnly(
    const char *pszFilename, const char *pszAccess)
{
    VSIFilesystemHandler *poFSHandler =
        VSIFileManager::GetHandler( pszFilename + strlen("/vsilz4/"));

    CPLMutexHolder oHolder(&hMutex);

#ifndef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    // Disable caching in fuzzing mode as the /vsilz4/ file is likely to
    // change very often
    // TODO: filename-based logic isn't enough. We should probably check
    // timestamp and/or file size.
    if( poHandleLastLz4File != nullptr &&
        strcmp(pszFilename + strlen("/vsilz4/"),
               poHandleLastLz4File->GetBaseFileName()) == 0 &&
        EQUAL(pszAccess, "rb") )
    {
        VSILz4Handle* poHandle = poHandleLastLz4File->Duplicate();
        if( poHandle )
            return poHandle;
    }
#else
    CPL_IGNORE_RET_VAL(pszAccess);
#endif

    VSIVirtualHandle* poVirtualHandle =
        poFSHandler->Open( pszFilename + strlen("/vsilz4/"), "rb" );

    if( poVirtualHandle == nullptr )
        return nullptr;

    unsigned char signature[2] = { '\0', '\0' };
    if( VSIFReadL(signature, 1, 2, reinterpret_cast<VSILFILE*>(poVirtualHandle)) != 2 ||
        signature[0] != gz_magic[0] || signature[1] != gz_magic[1] )
    {
        poVirtualHandle->Close();
        delete poVirtualHandle;
        return nullptr;
    }

    if( poHandleLastLz4File )
    {
        poHandleLastLz4File->UnsetCanSaveInfo();
        delete poHandleLastLz4File;
        poHandleLastLz4File = nullptr;
    }

    VSILz4Handle* poHandle =
        new VSILz4Handle(poVirtualHandle, pszFilename + strlen("/vsilz4/"));
    if( !(poHandle->IsInitOK()) )
    {
        delete poHandle;
        return nullptr;
    }
    return poHandle;
}

/************************************************************************/
/*                                Stat()                                */
/************************************************************************/

int VSILz4FilesystemHandler::Stat( const char *pszFilename,
                                    VSIStatBufL *pStatBuf,
                                    int nFlags )
{
    if( !STARTS_WITH_CI(pszFilename, "/vsilz4/") )
        return -1;

    CPLMutexHolder oHolder(&hMutex);

    memset(pStatBuf, 0, sizeof(VSIStatBufL));

    if( poHandleLastLz4File != nullptr &&
        strcmp(pszFilename+strlen("/vsilz4/"),
               poHandleLastLz4File->GetBaseFileName()) == 0 )
    {
        if( poHandleLastLz4File->GetUncompressedSize() != 0 )
        {
            pStatBuf->st_mode = S_IFREG;
            pStatBuf->st_size = poHandleLastLz4File->GetUncompressedSize();
            return 0;
        }
    }

    // Begin by doing a stat on the real file.
    int ret = VSIStatExL(pszFilename+strlen("/vsilz4/"), pStatBuf, nFlags);

    if( ret == 0 && (nFlags & VSI_STAT_SIZE_FLAG) )
    {
        CPLString osCacheFilename(pszFilename + strlen("/vsilz4/"));
        osCacheFilename += ".properties";

        // Can we save a bit of seeking by using a .properties file?
        VSILFILE* fpCacheLength = VSIFOpenL(osCacheFilename.c_str(), "rb");
        if( fpCacheLength )
        {
            const char* pszLine;
            GUIntBig nCompressedSize = 0;
            GUIntBig nUncompressedSize = 0;
            while( (pszLine = CPLReadLineL(fpCacheLength)) != nullptr )
            {
                if( STARTS_WITH_CI(pszLine, "compressed_size=") )
                {
                    const char* pszBuffer =
                        pszLine + strlen("compressed_size=");
                    nCompressedSize =
                        CPLScanUIntBig(pszBuffer,
                                       static_cast<int>(strlen(pszBuffer)));
                }
                else if( STARTS_WITH_CI(pszLine, "uncompressed_size=") )
                {
                    const char* pszBuffer =
                        pszLine + strlen("uncompressed_size=");
                    nUncompressedSize =
                        CPLScanUIntBig(pszBuffer,
                                       static_cast<int>(strlen(pszBuffer)));
                }
            }

            CPL_IGNORE_RET_VAL(VSIFCloseL(fpCacheLength));

            if( nCompressedSize == static_cast<GUIntBig>(pStatBuf->st_size) )
            {
                // Patch with the uncompressed size.
                pStatBuf->st_size = nUncompressedSize;

                VSILz4Handle* poHandle =
                    VSILz4FilesystemHandler::OpenLz4ReadOnly(pszFilename,
                                                               "rb");
                if( poHandle )
                {
                    poHandle->SetUncompressedSize(nUncompressedSize);
                    SaveInfo_unlocked(poHandle);
                    delete poHandle;
                }

                return ret;
            }
        }

        // No, then seek at the end of the data (slow).
        VSILz4Handle* poHandle =
                VSILz4FilesystemHandler::OpenLz4ReadOnly(pszFilename, "rb");
        if( poHandle )
        {
            poHandle->Seek(0, SEEK_END);
            const GUIntBig uncompressed_size =
                static_cast<GUIntBig>(poHandle->Tell());
            poHandle->Seek(0, SEEK_SET);

            // Patch with the uncompressed size.
            pStatBuf->st_size = uncompressed_size;

            delete poHandle;
        }
        else
        {
            ret = -1;
        }
    }

    return ret;
}

/************************************************************************/
/*                               Unlink()                               */
/************************************************************************/

int VSILz4FilesystemHandler::Unlink( const char * /* pszFilename */ )
{
    return -1;
}

/************************************************************************/
/*                               Rename()                               */
/************************************************************************/

int VSILz4FilesystemHandler::Rename( const char * /* oldpath */,
                                      const char * /* newpath */ )
{
    return -1;
}

/************************************************************************/
/*                               Mkdir()                                */
/************************************************************************/

int VSILz4FilesystemHandler::Mkdir( const char * /* pszDirname */,
                                     long /* nMode */ )
{
    return -1;
}
/************************************************************************/
/*                               Rmdir()                                */
/************************************************************************/

int VSILz4FilesystemHandler::Rmdir( const char * /* pszDirname */ )
{
    return -1;
}

/************************************************************************/
/*                             ReadDirEx()                                */
/************************************************************************/

char** VSILz4FilesystemHandler::ReadDirEx( const char * /*pszDirname*/,
                                            int /* nMaxFiles */ )
{
    return nullptr;
}

/************************************************************************/
/*                           GetOptions()                               */
/************************************************************************/

const char* VSILz4FilesystemHandler::GetOptions()
{
    return
    "<Options>"
    "  <Option name='GDAL_NUM_THREADS' type='string' "
        "description='Number of threads for compression. Either a integer or ALL_CPUS'/>"
    "  <Option name='CPL_VSIL_DEFLATE_CHUNK_SIZE' type='string' "
        "description='Chunk of uncompressed data for parallelization. "
        "Use K(ilobytes) or M(egabytes) suffix' default='1M'/>"
    "</Options>";
}


//! @endcond
/************************************************************************/
/*                   VSIInstallLz4FileHandler()                        */
/************************************************************************/

/**
 * \brief Install lz4 file system handler.
 *
 * A special file handler is installed that allows reading on-the-fly and
 * writing in lz4 (.lz4) files.
 *
 * All portions of the file system underneath the base
 * path "/vsilz4/" will be handled by this driver.
 *
 * Additional documentation is to be found at:
 * http://trac.osgeo.org/gdal/wiki/UserDocs/ReadInZip
 *
 * @since GDAL 3.1.0
 */

void VSIInstallLz4FileHandler()
{
    VSIFileManager::InstallHandler( "/vsilz4/", new VSILz4FilesystemHandler );
}
//! @cond Doxygen_Suppress
