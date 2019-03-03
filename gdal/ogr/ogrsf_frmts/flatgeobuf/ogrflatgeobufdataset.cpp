#include "ogr_flatgeobuf.h"

#include "flatgeobuf_generated.h"

static int OGRFlatGeobufDriverIdentify(GDALOpenInfo* poOpenInfo){
    if( STARTS_WITH_CI(poOpenInfo->pszFilename, "FGB:") )
        return TRUE;

    if( STARTS_WITH(poOpenInfo->pszFilename, "/vsicurl") )
    {
        if( CPLGetValueType(CPLGetFilename(poOpenInfo->pszFilename)) ==
                CPL_VALUE_INTEGER )
        {
            return TRUE;
        }
    }

    if( poOpenInfo->bIsDirectory )
    {
        if( CPLGetValueType(CPLGetFilename(poOpenInfo->pszFilename)) ==
                                                        CPL_VALUE_INTEGER )
        {
            // TODO: what is this?
            return FALSE;
        }
        return FALSE;
    }

    const auto nHeaderBytes = poOpenInfo->nHeaderBytes;
    const auto pabyHeader = poOpenInfo->pabyHeader;

    if(nHeaderBytes < 4)
        return FALSE;

    if( pabyHeader[0] == 0x66 &&
        pabyHeader[1] == 0x67 &&
        pabyHeader[2] == 0x62 ) {
        if (pabyHeader[3] == 0x00) {
            CPLDebug("FlatGeobuf", "Verified magicbytes");
            return TRUE;
        } else {
            CPLError(CE_Failure, CPLE_OpenFailed,
                "Unsupported FlatGeobuf version %d.\n",
                poOpenInfo->pabyHeader[3]);
        }

    }

    return FALSE;
}

void RegisterOGRFlatGeobuf()
{
    if( GDALGetDriverByName("FlatGeobuf") != NULL )
        return;

    GDALDriver *poDriver = new GDALDriver();
    poDriver->SetDescription("FlatGeobuf");
    poDriver->SetMetadataItem(GDAL_DCAP_VECTOR, "YES");
    poDriver->SetMetadataItem(GDAL_DMD_LONGNAME, "FlatGeobuf");
    poDriver->SetMetadataItem(GDAL_DMD_EXTENSION, "fgb");
    poDriver->SetMetadataItem(GDAL_DMD_HELPTOPIC, "drv_flatgeobuf.html");
    poDriver->SetMetadataItem(GDAL_DCAP_VIRTUALIO, "YES");
    poDriver->SetMetadataItem( GDAL_DMD_CREATIONFIELDDATATYPES,
                               "Integer Integer64 Real String IntegerList "
                               "Integer64List RealList StringList" );
    poDriver->SetMetadataItem( GDAL_DMD_CREATIONFIELDDATASUBTYPES, "Boolean" );

    poDriver->pfnOpen = OGRFlatGeobufDataset::Open;
    poDriver->pfnCreate = OGRFlatGeobufDataset::Create;
    poDriver->pfnIdentify = OGRFlatGeobufDriverIdentify;

    GetGDALDriverManager()->RegisterDriver(poDriver);
}


/************************************************************************/
/*                          OGRFlatGeobufDataset()                          */
/************************************************************************/

OGRFlatGeobufDataset::OGRFlatGeobufDataset()
{

}

OGRFlatGeobufDataset::OGRFlatGeobufDataset(const char *pszName)
{
    CPLDebug("FlatGeobuf", "Request to create dataset %s", pszName);
    m_create = true;
    m_pszName = pszName;
}

/************************************************************************/
/*                         ~OGRFlatGeobufDataset()                          */
/************************************************************************/

OGRFlatGeobufDataset::~OGRFlatGeobufDataset()
{
}

/************************************************************************/
/*                                Open()                                */
/************************************************************************/

GDALDataset *OGRFlatGeobufDataset::Open(GDALOpenInfo* poOpenInfo)
{
    if( !OGRFlatGeobufDriverIdentify(poOpenInfo) || poOpenInfo->eAccess == GA_Update )
        return nullptr;

    VSILFILE *fp = poOpenInfo->fpL;
    CPLString osFilename(poOpenInfo->pszFilename);

    uint64_t offset = 4;
    CPLDebug("FlatGeobuf", "Start at offset (%d)", 4);
    VSIFSeekL(fp, offset, SEEK_SET);
    uint32_t headerSize;
    VSIFReadL(&headerSize, 4, 1, fp);
    GByte* buf = static_cast<GByte*>(VSI_MALLOC_VERBOSE(headerSize));
    VSIFReadL(buf, 1, headerSize, fp);
    auto header = GetHeader(buf);
    auto featuresCount = header->features_count();
    auto treeSize = PackedRTree::size(featuresCount);
    offset += 4 + headerSize;
    CPLDebug("FlatGeobuf", "Add headerSize to offset (%d)", 4 + headerSize);
    offset += treeSize;
    CPLDebug("FlatGeobuf", "Add treeSize to offset (%zu)", treeSize);
    offset += featuresCount * 8;
    CPLDebug("FlatGeobuf", "Add featuresCount * 8 to offset (%zu)", featuresCount * 8);

    auto poDS = new OGRFlatGeobufDataset();
    poDS->SetDescription(osFilename);
    poDS->m_apoLayers.push_back(
        std::unique_ptr<OGRLayer>(new OGRFlatGeobufLayer(header, osFilename, offset))
    );

    return poDS;
}

GDALDataset *OGRFlatGeobufDataset::Create( const char *pszName,
                                        CPL_UNUSED int nBands,
                                        CPL_UNUSED int nXSize,
                                        CPL_UNUSED int nYSize,
                                        CPL_UNUSED GDALDataType eDT,
                                        char **papszOptions )
{
    // First, ensure there isn't any such file yet.
    VSIStatBufL sStatBuf;

    if (strcmp(pszName, "/dev/stdout") == 0)
        pszName = "/vsistdout/";

    if( VSIStatL(pszName, &sStatBuf) == 0 )
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "It seems a file system object called '%s' already exists.",
                 pszName);

        return nullptr;
    }

    // If the target is not a simple .fgb then create it as a directory.
    CPLString osDirName;

    if( EQUAL(CPLGetExtension(pszName), "fgb") )
    {
        osDirName = CPLGetPath(pszName);
        if( osDirName == "" )
            osDirName = ".";

        // HACK: CPLGetPath("/vsimem/foo.fgb") = "/vsimem", but this is not
        // recognized afterwards as a valid directory name.
        if( osDirName == "/vsimem" )
            osDirName = "/vsimem/";
    }
    else
    {
        if( STARTS_WITH(pszName, "/vsizip/"))
        {
            // Do nothing.
        }
        else if( !EQUAL(pszName, "/vsistdout/") &&
                 VSIMkdir(pszName, 0755) != 0 )
        {
            CPLError(CE_Failure, CPLE_AppDefined,
                     "Failed to create directory %s:\n%s",
                     pszName, VSIStrerror(errno));
            return nullptr;
        }
        osDirName = pszName;
    }

    if( EQUAL(CPLGetExtension(pszName), "fgb") )
    {
        return new OGRFlatGeobufDataset(pszName);
    }

    CPLError(CE_Failure, CPLE_AppDefined,
                     "Creating empty dataset not yet implemented");

    return nullptr;
}

OGRLayer* OGRFlatGeobufDataset::GetLayer( int iLayer ) {
    if( iLayer < 0 || iLayer >= GetLayerCount() )
        return nullptr;
    return m_apoLayers[iLayer].get();
}

int OGRFlatGeobufDataset::TestCapability( const char * pszCap )
{
    if (EQUAL(pszCap, ODrCCreateDataSource))
        return m_create;
    else if (EQUAL(pszCap, ODsCCreateLayer))
        return m_create;
    else if (EQUAL(pszCap, OLCSequentialWrite))
        return m_create;
    else if (EQUAL(pszCap, OLCCreateGeomField))
        return m_create;
    else if (EQUAL(pszCap, OLCIgnoreFields))
        return true;
    else if (EQUAL(pszCap, ODsCMeasuredGeometries))
        return true;
    else if (EQUAL(pszCap, OLCFastFeatureCount))
        return true;
    else if (EQUAL(pszCap, OLCFastGetExtent))
        return true;
    else if (EQUAL(pszCap, OLCFastSpatialFilter))
        return true;
    else
        return false;
}

GeometryType OGRFlatGeobufDataset::toGeometryType(OGRwkbGeometryType eGType)
{
    switch (eGType) {
        case OGRwkbGeometryType::wkbPoint: return GeometryType::Point;
        case OGRwkbGeometryType::wkbMultiPoint: return GeometryType::MultiPoint;
        case OGRwkbGeometryType::wkbLineString: return GeometryType::LineString;
        case OGRwkbGeometryType::wkbMultiLineString: return GeometryType::MultiLineString;
        case OGRwkbGeometryType::wkbPolygon: return GeometryType::Polygon;
        default:
            throw std::runtime_error("Unknown geometry type");
    }
}

OGRwkbGeometryType OGRFlatGeobufDataset::toOGRwkbGeometryType(GeometryType geometryType)
{
    switch (geometryType) {
        case GeometryType::Point: return OGRwkbGeometryType::wkbPoint;
        case GeometryType::MultiPoint: return OGRwkbGeometryType::wkbMultiPoint;
        case GeometryType::LineString: return OGRwkbGeometryType::wkbLineString;
        case GeometryType::MultiLineString: return OGRwkbGeometryType::wkbMultiLineString;
        case GeometryType::Polygon: return OGRwkbGeometryType::wkbPolygon;
        default:
            throw std::runtime_error("Unknown geometry type");
    }
}

OGRLayer* OGRFlatGeobufDataset::ICreateLayer( const char *pszLayerName,
                                OGRSpatialReference *poSpatialRef,
                                OGRwkbGeometryType eGType,
                                char **papszOptions )
{
    // Verify we are in update mode.
    if( !m_create )
    {
        CPLError(CE_Failure, CPLE_NoWriteAccess,
                 "Data source %s opened read-only.\n"
                 "New layer %s cannot be created.",
                 m_pszName, pszLayerName);

        return nullptr;
    }

    // Verify that the datasource is a directory.
    VSIStatBufL sStatBuf;

    /*if( STARTS_WITH(m_pszName, "/vsizip/"))
    {
        // Do nothing.
    }
    else if( !EQUAL(m_pszName, "/vsistdout/") &&
             (VSIStatL(m_pszName, &sStatBuf) != 0 ||
              !VSI_ISDIR(sStatBuf.st_mode)) )
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to create FlatGeobuf layer (file) against a "
                 "non-directory datasource.");
        return nullptr;
    }*/

    // What filename would we use?
    CPLString osFilename;

    CPLDebug("FlatGeobuf", "m_pszName: %s", m_pszName);
    CPLDebug("FlatGeobuf", "pszLayerName: %s", pszLayerName);
    osFilename = CPLFormFilename(m_pszName, pszLayerName, "fgb");
    CPLDebug("FlatGeobuf", "osFilename: %s", osFilename.c_str());

    // Does this directory/file already exist?
    if( VSIStatL(osFilename, &sStatBuf) == 0 )
    {
        CPLError(CE_Failure, CPLE_AppDefined,
                 "Attempt to create layer %s, but %s already exists.",
                 pszLayerName, osFilename.c_str());
        return nullptr;
    }

    // Create a layer.
    OGRFlatGeobufLayer *poLayer = new OGRFlatGeobufLayer(pszLayerName, m_pszName, poSpatialRef, eGType);

    m_apoLayers.push_back(
        std::unique_ptr<OGRLayer>(poLayer)
    );

    return poLayer;
}
