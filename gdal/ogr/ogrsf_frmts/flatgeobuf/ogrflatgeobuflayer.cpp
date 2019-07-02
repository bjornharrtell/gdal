#include "ogrsf_frmts.h"
#include "cpl_vsi_virtual.h"
#include "cpl_conv.h"
#include "cpl_json.h"
#include "cpl_http.h"
#include "ogr_p.h"

#include "ogr_flatgeobuf.h"

using namespace flatbuffers;
using namespace FlatGeobuf;

OGRFlatGeobufLayer::OGRFlatGeobufLayer(const Header *poHeader, const char* pszFilename, uint64_t offset)
{
    CPLDebug("FlatGeobuf", "offset: %zu", offset);

    if (poHeader == nullptr)
        CPLError(CE_Fatal, CPLE_OpenFailed, "poHeader is null.\n");

    if (pszFilename == nullptr)
        CPLError(CE_Fatal, CPLE_OpenFailed, "pszFilename is null.\n");

    m_poHeader = poHeader;
    // TODO: free
    m_pszFilename = CPLStrdup(pszFilename);
    m_offsetInit = offset;
    m_offset = offset;
    m_create = false;

    m_featuresCount = m_poHeader->features_count();
    m_geometryType = m_poHeader->geometry_type();
    m_hasM = m_poHeader->hasM();
    m_hasZ = m_poHeader->hasZ();
    m_hasT = m_poHeader->hasT();

    auto crs = m_poHeader->crs();
    if (crs != nullptr) {
        auto code = crs->code();
        if (code != 0) {
            m_poSRS = new OGRSpatialReference();
            m_poSRS->importFromEPSG(code);
        }
    }

    auto eGType = getOGRwkbGeometryType();

    CPLDebug("FlatGeobuf", "m_featuresCount: %zu", m_featuresCount);

    //CPLError(CE_Failure, CPLE_OpenFailed, "m_featuresCount is %d.\n", m_featuresCount);
    //CPLError(CE_Failure, CPLE_OpenFailed, "m_poHeader->name()->c_str() is %s.\n", m_poHeader->name()->c_str());

    const char *pszName = nullptr;
    if (m_poHeader->name()) {
        pszName = m_poHeader->name()->c_str();
        CPLDebug("FlatGeobuf", "header name: %s", pszName);
    }
    m_poFeatureDefn = new OGRFeatureDefn(pszName);
    SetDescription(m_poFeatureDefn->GetName());
    m_poFeatureDefn->SetGeomType(wkbNone);
    OGRGeomFieldDefn *poGeomFieldDefn = new OGRGeomFieldDefn(nullptr, eGType);
    if (m_poSRS != nullptr)
        poGeomFieldDefn->SetSpatialRef(m_poSRS);
    m_poFeatureDefn->AddGeomFieldDefn(poGeomFieldDefn, false);
    readColumns();
    m_poFeatureDefn->Reference();
}

OGRFlatGeobufLayer::OGRFlatGeobufLayer(
    const char *pszLayerName,
    const char *pszFilename,
    OGRSpatialReference *poSpatialRef,
    OGRwkbGeometryType eGType)
{
    CPLDebug("FlatGeobuf", "Request to create layer %s", pszLayerName);

    // TODO: free
    m_pszLayerName = CPLStrdup(pszLayerName);
    m_pszFilename = CPLStrdup(pszFilename);
    m_create = true;
    m_eGType = eGType;
    translateOGRwkbGeometryType();
    m_poSRS = poSpatialRef;

    CPLDebug("FlatGeobuf", "eGType: %d", (int) eGType);
    CPLDebug("FlatGeobuf", "m_geometryType: %d", (int) m_geometryType);

    SetMetadataItem(OLMD_FID64, "YES");

    m_poFeatureDefn = new OGRFeatureDefn(pszLayerName);
    SetDescription(m_poFeatureDefn->GetName());
    m_poFeatureDefn->SetGeomType(eGType);
    m_poFeatureDefn->Reference();
}

void OGRFlatGeobufLayer::translateOGRwkbGeometryType()
{
    switch (wkbFlatten(m_eGType)) {
        case OGRwkbGeometryType::wkbPoint: m_geometryType = GeometryType::Point; break;
        case OGRwkbGeometryType::wkbMultiPoint: m_geometryType = GeometryType::MultiPoint; break;
        case OGRwkbGeometryType::wkbLineString: m_geometryType = GeometryType::LineString; break;
        case OGRwkbGeometryType::wkbMultiLineString: m_geometryType = GeometryType::MultiLineString; break;
        case OGRwkbGeometryType::wkbPolygon: m_geometryType = GeometryType::Polygon; break;
        case OGRwkbGeometryType::wkbMultiPolygon: m_geometryType = GeometryType::MultiPolygon; break;
        default:
            CPLError(CE_Fatal, CPLE_NotSupported, "toGeometryType: Unknown OGRwkbGeometryType %d", (int) m_eGType);
    }
    if wkbHasZ(m_eGType)
        m_hasZ = true;
    if wkbHasM(m_eGType)
        m_hasM = true;
    return;
}

OGRwkbGeometryType OGRFlatGeobufLayer::getOGRwkbGeometryType()
{
    switch (m_geometryType) {
        case GeometryType::Point: return OGRwkbGeometryType::wkbPoint;
        case GeometryType::MultiPoint: return OGRwkbGeometryType::wkbMultiPoint;
        case GeometryType::LineString: return OGRwkbGeometryType::wkbLineString;
        case GeometryType::MultiLineString: return OGRwkbGeometryType::wkbMultiLineString;
        case GeometryType::Polygon: return OGRwkbGeometryType::wkbPolygon;
        case GeometryType::MultiPolygon: return OGRwkbGeometryType::wkbMultiPolygon;
        default:
            CPLError(CE_Fatal, CPLE_NotSupported, "toOGRwkbGeometryType: Unknown FlatGeobuf::GeometryType %d", (int) m_geometryType);
    }
    return OGRwkbGeometryType::wkbUnknown;
}

ColumnType OGRFlatGeobufLayer::toColumnType(OGRFieldType type, OGRFieldSubType /* subType */)
{
    switch (type) {
        case OGRFieldType::OFTInteger: return ColumnType::Int;
        case OGRFieldType::OFTInteger64: return ColumnType::Long;
        case OGRFieldType::OFTReal: return ColumnType::Double;
        case OGRFieldType::OFTString: return ColumnType::String;
        default: CPLError(CE_Fatal, CPLE_AppDefined, "toColumnType: Unknown OGRFieldType %d", type);
    }
    return ColumnType::String;
}

OGRFieldType OGRFlatGeobufLayer::toOGRFieldType(ColumnType type)
{
    switch (type) {
        case ColumnType::Int: return OGRFieldType::OFTInteger;
        case ColumnType::Long: return OGRFieldType::OFTInteger64;
        case ColumnType::Double: return OGRFieldType::OFTReal;
        case ColumnType::String: return OGRFieldType::OFTString;
        default: CPLError(CE_Fatal, CPLE_AppDefined, "toOGRFieldType: Unknown ColumnType %d", (int) type);
    }
    return OGRFieldType::OFTString;
}

const std::vector<Offset<Column>> OGRFlatGeobufLayer::writeColumns(FlatBufferBuilder &fbb)
{
    std::vector<Offset<Column>> columns;
    for (int i = 0; i < m_poFeatureDefn->GetFieldCount(); i++) {
        auto field = m_poFeatureDefn->GetFieldDefn(i);
        auto name = field->GetNameRef();
        auto columnType = toColumnType(field->GetType(), field->GetSubType());
        auto column = CreateColumnDirect(fbb, name, columnType);
        columns.push_back(column);
    }
    return columns;
}

void OGRFlatGeobufLayer::readColumns()
{
    auto columns = m_poHeader->columns();
    if (columns == nullptr)
        return;
    for (size_t i = 0; i < columns->size(); i++) {
        auto column = columns->Get(i);
        auto name = column->name()->c_str();
        auto type = toOGRFieldType(column->type());
        OGRFieldDefn field(name, type);
        m_poFeatureDefn->AddFieldDefn(&field);
    }
}

OGRFlatGeobufLayer::~OGRFlatGeobufLayer()
{
    if (m_create) {
        CPLDebug("FlatGeobuf", "Request to create %zu features", m_featuresCount);
        size_t c;

        //const char *filename = CPLFormFilename("", m_pszLayerName, "fgb");
        CPLDebug("FlatGeobuf", "m_pszFilename: %s", m_pszFilename);
        VSILFILE *fp = VSIFOpenL(m_pszFilename, "wb");

        if(fp == nullptr) {
            CPLError(CE_Fatal, CPLE_OpenFailed,
                        "Failed to create %s:\n%s",
                        m_pszFilename, VSIStrerror(errno));
            return;
        }

        c = VSIFWriteL(&magicbytes, sizeof(magicbytes), 1, fp);
        CPLDebug("FlatGeobuf", "Wrote magicbytes (%zu bytes)", c * sizeof(magicbytes));

        Rect extent = calcExtent(m_featureItems);
        const auto extentVector = extent.toVector();

        FlatBufferBuilder fbb;
        auto columns = writeColumns(fbb);

        uint16_t indexNodeSize = bCreateSpatialIndexAtClose ? 16 : 0;

        auto crs = 0;
        // TODO: this can crash for some inputs
        /*if (m_poSRS != nullptr) {
            auto code = m_poSRS->GetEPSGGeogCS();
            if (code != -1) {
                CPLDebug("FlatGeobuf", "Creating SRS with EPSG code %d", code);
                crs = CreateCrsDirect(fbb, code);
            }
        }*/

        auto header = CreateHeaderDirect(
            fbb, m_pszLayerName, &extentVector, m_geometryType, m_hasZ, m_hasM, m_hasT, &columns, m_featuresCount, true, indexNodeSize, 0);
        fbb.FinishSizePrefixed(header);
        c = VSIFWriteL(fbb.GetBufferPointer(), 1, fbb.GetSize(), fp);
        CPLDebug("FlatGeobuf", "Wrote header (%zu bytes)", c);

        if (bCreateSpatialIndexAtClose) {
            CPLDebug("FlatGeobuf", "Sorting items for Packed R-tree");
            hilbertSort(m_featureItems);
            CPLDebug("FlatGeobuf", "Creating Packed R-tree");
            PackedRTree tree(m_featureItems, extent);
            CPLDebug("FlatGeobuf", "PackedRTree extent %f, %f, %f, %f", extentVector[0], extentVector[1], extentVector[2], extentVector[3]);
            c = VSIFWriteL(tree.toData(), 1, tree.size(), fp);
            CPLDebug("FlatGeobuf", "Wrote tree (%zu bytes)", c);
        }

        c = 0;
        for (uint64_t i = 0, offset = 0; i < m_featuresCount; i++) {
            c += VSIFWriteL(&offset, 8, 1, fp);
            offset += static_cast<FeatureItem *>(m_featureItems[i])->size;
        }
        CPLDebug("FlatGeobuf", "Wrote feature offsets (%zu bytes)", c * 8);

        c = 0;
        for (uint64_t i = 0; i < m_featuresCount; i++)
            c += VSIFWriteL(static_cast<FeatureItem *>(m_featureItems[i])->data, 1, static_cast<FeatureItem *>(m_featureItems[i])->size, fp);
        CPLDebug("FlatGeobuf", "Wrote feature buffers (%zu bytes)", c);

        VSIFCloseL(fp);
    }

    if (m_poFeatureDefn != nullptr)
        m_poFeatureDefn->Release();

    if (m_padfX != nullptr)
        CPLFree(m_padfX);
    if (m_padfY != nullptr)
        CPLFree(m_padfY);
    if (m_padfZ != nullptr)
        CPLFree(m_padfZ);
    if (m_padfM != nullptr)
        CPLFree(m_padfM);

    if (m_featureBuf != nullptr)
        VSIFree(m_featureBuf);

    if (m_featureOffsets != nullptr)
        VSIFree(m_featureOffsets);

    m_processedSpatialIndex = false;
}

OGRFeature *OGRFlatGeobufLayer::GetFeature(GIntBig /* nFeatureId */)
{
    CPLError(CE_Fatal, CPLE_AppDefined, "GetFeature: Not implemented");
    return nullptr;
}

void OGRFlatGeobufLayer::processSpatialIndex() {
    if (m_poFilterGeom != nullptr && !m_processedSpatialIndex) {
        CPLDebug("FlatGeobuf", "Will do spatial index search");
        m_processedSpatialIndex = true;
        if (m_poFp == nullptr) {
            CPLDebug("FlatGeobuf", "processSpatialIndex (will attempt to open file %s)", m_pszFilename);
            m_poFp = VSIFOpenL(m_pszFilename, "rb");
            //m_poFp = (VSILFILE*) VSICreateCachedFile ( (VSIVirtualHandle*) m_poFp);
        }
        VSIFSeekL(m_poFp, sizeof(magicbytes), SEEK_SET); // skip magic bytes
        uoffset_t headerSize;
        VSIFReadL(&headerSize, sizeof(uoffset_t), 1, m_poFp);
        auto featuresCount = m_poHeader->features_count();
        auto treeSize = PackedRTree::size(featuresCount);
        OGREnvelope env;
        m_poFilterGeom->getEnvelope(&env);
        Rect r { env.MinX, env.MinY, env.MaxX, env.MaxY };
        CPLDebug("FlatGeobuf", "Spatial index search on %f,%f,%f,%f", env.MinX, env.MinY, env.MaxX, env.MaxY);
        auto readNode = [this, headerSize] (uint8_t *buf, uint32_t i, uint32_t s) {
            VSIFSeekL(m_poFp, sizeof(magicbytes) + sizeof(uoffset_t) + headerSize + i, SEEK_SET);
            VSIFReadL(buf, 1, s, m_poFp);
        };
        m_foundFeatureIndices = PackedRTree::streamSearch(featuresCount, 16, r, readNode);
        m_featuresCount = m_foundFeatureIndices.size();
        if (m_featuresCount == 0) {
            CPLDebug("FlatGeobuf", "No found features in spatial index search");
            return;
        }
        CPLDebug("FlatGeobuf", "%zu features found in spatial index search", m_featuresCount);
        m_featureOffsets = static_cast<uint64_t *>(VSI_MALLOC_VERBOSE(featuresCount * 8));
        VSIFSeekL(m_poFp, sizeof(magicbytes) + sizeof(uoffset_t) + headerSize + treeSize, SEEK_SET);
        VSIFReadL(m_featureOffsets, 8, featuresCount, m_poFp);
    } else {
        CPLDebug("FlatGeobuf", "processSpatialIndex noop");
    }
}

GIntBig OGRFlatGeobufLayer::GetFeatureCount(int /* bForce */) {
    processSpatialIndex();
    return m_featuresCount;
}

OGRFeature *OGRFlatGeobufLayer::GetNextFeature()
{
    if (m_featuresPos >= m_featuresCount) {
        CPLDebug("FlatGeobuf", "Iteration end");
        if (m_poFp != nullptr) {
            VSIFCloseL(m_poFp);
            m_poFp = nullptr;
        }
        return nullptr;
    }

    if (m_poFp == nullptr) {
        CPLDebug("FlatGeobuf", "Iteration start (will attempt to open file %s)", m_pszFilename);
        m_poFp = VSIFOpenL(m_pszFilename, "rb");
        //m_poFp = (VSILFILE*) VSICreateCachedFile ( (VSIVirtualHandle*) m_poFp);

        processSpatialIndex();
        if (m_featuresCount == 0) {
            CPLDebug("FlatGeobuf", "No features found");
            if (m_poFp != nullptr) {
                VSIFCloseL(m_poFp);
                m_poFp = nullptr;
            }
            return nullptr;
        }
    }

    if (m_processedSpatialIndex)
        m_offset = m_offsetInit + m_featureOffsets[m_foundFeatureIndices[m_featuresPos]];
    else if (m_featuresPos > 0)
        m_offset += m_featureSize + sizeof(uoffset_t);

    OGRFeature* poFeature = new OGRFeature(m_poFeatureDefn);

    VSIFSeekL(m_poFp, m_offset, SEEK_SET);
    VSIFReadL(&m_featureSize, sizeof(uoffset_t), 1, m_poFp);
    if (m_featureBufSize == 0) {
        m_featureBufSize = std::max(1024U * 32U, m_featureSize);
        CPLDebug("FlatGeobuf", "m_featureBufSize: %d", m_featureBufSize);
        m_featureBuf = static_cast<GByte *>(VSI_MALLOC_VERBOSE(m_featureBufSize));
    } else if (m_featureBufSize < m_featureSize) {
        m_featureBufSize = std::max(m_featureBufSize * 2, m_featureSize);
        CPLDebug("FlatGeobuf", "m_featureBufSize: %d", m_featureBufSize);
        m_featureBuf = static_cast<GByte *>(VSI_REALLOC_VERBOSE(m_featureBuf, m_featureBufSize));
    }
    VSIFReadL(m_featureBuf, 1, m_featureSize, m_poFp);

#ifdef DEBUG
    const uint8_t * vBuf = const_cast<const uint8_t *>(reinterpret_cast<uint8_t *>(m_featureBuf));
    Verifier v(vBuf, m_featureSize);
    auto ok = VerifyFeatureBuffer(v);
    if (!ok) {
        CPLDebug("FlatGeobuf", "VerifyFeatureBuffer says not ok");
        CPLDebug("FlatGeobuf", "m_offset: %zu", m_offset);
        CPLDebug("FlatGeobuf", "m_featuresPos: %zu", m_featuresPos);
        CPLDebug("FlatGeobuf", "featureSize: %d", m_featureSize);
    }
#endif

    auto feature = GetRoot<Feature>(m_featureBuf);
    auto fid = feature->fid();
    poFeature->SetFID(fid);
    //CPLDebug("FlatGeobuf", "fid: %zu", fid);
    auto ogrGeometry = readGeometry(feature);
#ifdef DEBUG
    //char *wkt;
    //ogrGeometry->exportToWkt(&wkt);
    //CPLDebug("FlatGeobuf", "readGeometry as wkt: %s", wkt);
#endif
    // TODO: find out why this is done in other drivers
    //if (poSRS != nullptr)
    //    ogrGeometry->assignSpatialReference(poSRS);
    poFeature->SetGeometry(ogrGeometry);

    auto properties = feature->properties();

    if (properties != nullptr) {
        auto data = properties->data();
        auto size = properties->size();
        //CPLDebug("FlatGeobuf", "properties->size: %d", size);
        uoffset_t offset = 0;
        while (offset < (size-1)) {
            uint16_t i = *((uint16_t *)(data + offset));
            offset += sizeof(uint16_t);
            //CPLDebug("FlatGeobuf", "i: %d", i);
            auto column = m_poHeader->columns()->Get(i);
            auto type = column->type();
            auto ogrField = poFeature->GetRawFieldRef(i);
            switch (type) {
                case ColumnType::Int:
                    ogrField->Integer = *((int32_t *)(data + offset));
                    offset += sizeof(int32_t);
                    break;
                case ColumnType::Long:
                    ogrField->Integer64 = *((int64_t *)(data + offset));
                    offset += sizeof(int64_t);
                    break;
                case ColumnType::Double:
                    ogrField->Real = *((double *)(data + offset));
                    offset += sizeof(double);
                    break;
                case ColumnType::String: {
                    uint32_t len = *((uint32_t *)(data + offset));
                    offset += sizeof(uint32_t);
                    uint8_t *str = new uint8_t[len + 1];
                    memcpy(str, data + offset, len);
                    offset += len;
                    str[len] = '\0';
                    ogrField->String = (char *) str;
                    break;
                }
                default:
                    CPLError(CE_Fatal, CPLE_AppDefined, "GetNextFeature: Unknown column->type: %d", (int) type);
            }
        }
    }

    m_featuresPos++;
    return poFeature;
}

void OGRFlatGeobufLayer::ensurePadfBuffers(size_t count)
{
    size_t requiredSize = count * sizeof(double);
    if (m_padfSize == 0) {
        m_padfSize = std::max(1024 * sizeof(double), requiredSize);
        CPLDebug("FlatGeobuf", "m_padfSize: %zu", m_padfSize);
        m_padfX = static_cast<double *>(CPLMalloc(m_padfSize));
        m_padfY = static_cast<double *>(CPLMalloc(m_padfSize));
        if (m_hasZ)
            m_padfZ = static_cast<double *>(CPLMalloc(m_padfSize));
        if (m_hasM)
            m_padfM = static_cast<double *>(CPLMalloc(m_padfSize));
    } else if (m_padfSize < requiredSize) {
        m_padfSize = std::max(m_padfSize * 2, requiredSize);
        CPLDebug("FlatGeobuf", "m_padfSize: %zu", m_padfSize);
        m_padfX = static_cast<double *>(CPLRealloc(m_padfX, m_padfSize));
        m_padfY = static_cast<double *>(CPLRealloc(m_padfY, m_padfSize));
        if (m_hasZ)
            m_padfZ = static_cast<double *>(CPLRealloc(m_padfZ, m_padfSize));
        if (m_hasM)
            m_padfM = static_cast<double *>(CPLRealloc(m_padfM, m_padfSize));
    }
}

OGRPoint *OGRFlatGeobufLayer::readPoint(const double *coords, uint32_t offset)
{
    return new OGRPoint { coords[offset + 0], coords[offset + 1] };
    // TODO: Z M ??
    /*
    else if (dimensions == 3)
        return new OGRPoint { coords[offset + 0], coords[offset + 1], coords[offset + 2] };
    else if (dimensions == 4)
        return new OGRPoint { coords[offset + 0], coords[offset + 1], coords[offset + 2], coords[offset + 3] };
    */
    CPLError(CE_Fatal, CPLE_AppDefined, "readPoint: Unsupported number of dimensions");
    return nullptr;
}

OGRMultiPoint *OGRFlatGeobufLayer::readMultiPoint(const double *coords, uint32_t coordsLength)
{
    auto mp = new OGRMultiPoint();
    for (size_t i = 0; i < coordsLength; i = i + 2)
        mp->addGeometryDirectly(readPoint(coords, i));
    return mp;
}

OGRLineString *OGRFlatGeobufLayer::readLineString(const double *coords, uint32_t coordsLength, uint32_t offset)
{
    auto ls = new OGRLineString();

    ensurePadfBuffers(coordsLength >> 1);

    unsigned int c = 0;
    for (size_t i = offset; i < offset + coordsLength; i = i + 2) {
        m_padfX[c] = coords[i];
        m_padfY[c] = coords[i+1];
        if (m_hasZ)
            m_padfZ[c] = coords[i+2];
        if (m_hasM)
            m_padfM[c] = coords[i+3];
        c++;
    }
    ls->setNumPoints(coordsLength >> 1, 0);
    ls->setPoints(coordsLength >> 1, m_padfX, m_padfY);
    // TODO: handle hasZ hasM
    /*else if (dimensions == 3)
        ls->setPoints(dimLength, m_padfX, m_padfY, m_padfZ);
    else if (dimensions == 4)
        ls->setPoints(dimLength, m_padfX, m_padfY, m_padfZ, m_padfM);
    */
    return ls;
}

OGRMultiLineString *OGRFlatGeobufLayer::readMultiLineString(const double *coords, const flatbuffers::Vector<uint32_t> *ends)
{
    auto mls = new OGRMultiLineString();
    uint32_t offset = 0;
    for (size_t i = 0; i < ends->size(); i++) {
        auto end = ends->Get(i);
        auto ls = readLineString(coords, end, offset);
        mls->addGeometryDirectly(ls);
        offset = end;
    }
    return mls;
}

OGRLinearRing *OGRFlatGeobufLayer::readLinearRing(const double *coords, uint32_t coordsLength, uint32_t offset)
{
    auto ls = new OGRLinearRing();

    ensurePadfBuffers(coordsLength);

    unsigned int c = 0;
    for (size_t i = offset; i < offset + coordsLength; i = i + 2) {
        m_padfX[c] = coords[i];
        m_padfY[c] = coords[i+1];
        if (m_hasZ)
            m_padfZ[c] = coords[i+2];
        if (m_hasM)
            m_padfM[c] = coords[i+3];
        c++;
    }
    ls->setNumPoints(coordsLength >> 1, 0);
    ls->setPoints(coordsLength >> 1, m_padfX, m_padfY);
    // TODO: handle Z / M
    /*
    else if (dimensions == 3)
        ls->setPoints(dimLength, m_padfX, m_padfY, m_padfZ);
    else if (dimensions == 4)
        ls->setPoints(dimLength, m_padfX, m_padfY, m_padfZ, m_padfM);
    */
    return ls;
}

OGRPolygon *OGRFlatGeobufLayer::readPolygon(const double *coords, uint32_t coordsLength, const Vector<uint32_t> *ends, uint32_t offset)
{
    auto p = new OGRPolygon();
    if (ends == nullptr || ends->size() < 2) {
        p->addRingDirectly(readLinearRing(coords, coordsLength));
    } else {
        for (size_t i = 0; i < ends->size(); i++) {
            auto end = ends->Get(i);
            p->addRingDirectly(readLinearRing(coords, end - offset, offset));
            offset = end;
        }
    }
    return p;
}

OGRMultiPolygon *OGRFlatGeobufLayer::readMultiPolygon(
    const double *coords,
    uint32_t coordsLength,
    const Vector<uint32_t> *ends,
    const Vector<uint32_t> *endss)
{
    auto mp = new OGRMultiPolygon();
    if (endss == nullptr || endss->size() < 2) {
        mp->addGeometryDirectly(readPolygon(coords, coordsLength, ends));
    } else {
        uint32_t offset = 0;
        size_t roffset = 0;
        for (size_t i = 0; i < endss->size(); i++) {
            auto p = new OGRPolygon();
            uint32_t ringCount = endss->Get(i);
            for (size_t j = 0; j < ringCount; j++) {
                uint32_t end = ends->Get(roffset++);
                p->addRingDirectly(readLinearRing(coords, end - offset, offset));
                offset = end;
            }
            mp->addGeometryDirectly(p);
        }
    }
    return mp;
}

OGRGeometry *OGRFlatGeobufLayer::readGeometry(const Feature *feature)
{
    auto pXy = feature->xy();
    if (pXy == nullptr)
        CPLError(CE_Fatal, CPLE_AppDefined, "readGeometry: Geometry has no coordinates");
    auto xy = pXy->data();
    auto xySize = pXy->size();
    switch (m_geometryType) {
        case GeometryType::Point:
            return readPoint(xy);
        case GeometryType::MultiPoint:
            return readMultiPoint(xy, xySize);
        case GeometryType::LineString:
            return readLineString(xy, xySize);
        case GeometryType::MultiLineString:
            return readMultiLineString(xy, feature->ends());
        case GeometryType::Polygon:
            return readPolygon(xy, xySize, feature->ends());
        case GeometryType::MultiPolygon:
            return readMultiPolygon(xy, xySize, feature->ends(), feature->endss());
        default:
            CPLError(CE_Fatal, CPLE_AppDefined, "readGeometry: Unknown FlatGeobuf::GeometryType %d", (int) m_geometryType);
    }
    return nullptr;
}

OGRErr OGRFlatGeobufLayer::CreateField(OGRFieldDefn *poField, int /* bApproxOK */)
{
    CPLDebug("FlatGeobuf", "CreateField %s %s", poField->GetNameRef(), poField->GetFieldTypeName(poField->GetType()));
    if(!TestCapability(OLCCreateField))
    {
        CPLError(CE_Failure, CPLE_AppDefined, "Unable to create new fields after first feature written.");
        return OGRERR_FAILURE;
    }

    m_poFeatureDefn->AddFieldDefn(poField);

    return OGRERR_NONE;
}

OGRErr OGRFlatGeobufLayer::ICreateFeature(OGRFeature *poNewFeature)
{
    auto fid = poNewFeature->GetFID();
    if (fid == OGRNullFID)
        fid = m_featuresCount;

    uint8_t *propertiesBuffer = new uint8_t[1000000];
    uint32_t propertiesOffset = 0;
    FlatBufferBuilder fbb;

    for (int i = 0; i < m_poFeatureDefn->GetFieldCount(); i++) {
        auto fieldDef = m_poFeatureDefn->GetFieldDefn(i);
        if (poNewFeature->IsFieldNull(i))
            continue;

        uint16_t column_index = i;
        memcpy(propertiesBuffer + propertiesOffset, &column_index, sizeof(uint16_t));
        propertiesOffset += sizeof(uint16_t);

        auto fieldType = fieldDef->GetType();
        auto field = poNewFeature->GetRawFieldRef(i);
        switch (fieldType) {
            case OGRFieldType::OFTInteger: {
                memcpy(propertiesBuffer + propertiesOffset, &field->Integer, sizeof(int32_t));
                propertiesOffset += sizeof(int32_t);
                break;
            }
            case OGRFieldType::OFTInteger64: {
                memcpy(propertiesBuffer + propertiesOffset, &field->Integer64, sizeof(int64_t));
                propertiesOffset += sizeof(int64_t);
                break;
            }
            case OGRFieldType::OFTReal: {
                memcpy(propertiesBuffer + propertiesOffset, &field->Real, sizeof(double));
                propertiesOffset += sizeof(double);
                break;
            }
            case OGRFieldType::OFTString: {
                uint32_t len = strlen(field->String);
                memcpy(propertiesBuffer + propertiesOffset, &len, sizeof(uint32_t));
                propertiesOffset += sizeof(len);
                memcpy(propertiesBuffer + propertiesOffset, field->String, len);
                propertiesOffset += len;
                break;
            }
            default:
                CPLError(CE_Failure, CPLE_AppDefined, "ICreateFeature: Missing implementation for OGRFieldType %d", fieldType);
                return OGRERR_FAILURE;
        }
    }

    auto ogrGeometry = poNewFeature->GetGeometryRef();
#ifdef DEBUG
    //char *wkt;
    //ogrGeometry->exportToWkt(&wkt);
    //CPLDebug("FlatGeobuf", "poNewFeature as wkt: %s", wkt);
#endif
    if (ogrGeometry == nullptr)
        return 0;
    if (ogrGeometry->getGeometryType() != m_eGType) {
        CPLError(CE_Failure, CPLE_AppDefined, "ICreateFeature: Mismatched geometry type");
        return OGRERR_FAILURE;
    }

    std::vector<double> xy;
    std::vector<uint32_t> ends;
    std::vector<uint32_t> endss;
    switch (m_geometryType) {
        case GeometryType::Point:
            writePoint(ogrGeometry->toPoint(), xy);
            break;
        case GeometryType::MultiPoint:
            writeMultiPoint(ogrGeometry->toMultiPoint(), xy);
            break;
        case GeometryType::LineString:
            writeLineString(ogrGeometry->toLineString(), xy);
            break;
        case GeometryType::MultiLineString:
            writeMultiLineString(ogrGeometry->toMultiLineString(), xy, ends);
            break;
        case GeometryType::Polygon:
            writePolygon(ogrGeometry->toPolygon(), xy, ends, false, 0);
            break;
        case GeometryType::MultiPolygon:
            writeMultiPolygon(ogrGeometry->toMultiPolygon(), xy, ends, endss);
            break;
        default:
            CPLError(CE_Failure, CPLE_AppDefined, "ICreateFeature: Unknown FlatGeobuf::GeometryType %d", (int) m_geometryType);
            return OGRERR_FAILURE;
    }
    CPLDebug("FlatGeobuf", "geom encoded");
    auto pEnds = ends.size() == 0 ? nullptr : &ends;
    auto pEndss = endss.size() == 0 ? nullptr : &endss;
    std::vector<uint8_t> properties ( propertiesBuffer, propertiesBuffer + propertiesOffset );
    auto feature = CreateFeatureDirect(fbb, fid, pEnds, pEndss, &xy, nullptr, nullptr, nullptr, &properties);
    //auto feature = CreateFeatureDirect(fbb, fid, pEnds, pEndss, &coords, nullptr);
    fbb.FinishSizePrefixed(feature);
    delete propertiesBuffer;

    OGREnvelope psEnvelope;
    ogrGeometry->getEnvelope(&psEnvelope);

    auto item = new FeatureItem();
    item->buf = fbb.Release();
    item->data = item->buf.data();
    item->size = item->buf.size();
    item->rect = {
        psEnvelope.MinX,
        psEnvelope.MinY,
        psEnvelope.MaxX,
        psEnvelope.MaxY
    };

    m_featureItems.push_back(item);

    m_featuresCount++;

    return OGRERR_NONE;
}

void OGRFlatGeobufLayer::writePoint(OGRPoint *p, std::vector<double> &coords)
{
    CPLDebug("FlatGeobuf", "writePoint %f %f", p->getX(), p->getY());
    coords.push_back(p->getX());
    coords.push_back(p->getY());
}

void OGRFlatGeobufLayer::writeMultiPoint(OGRMultiPoint *mp, std::vector<double> &coords)
{
    for (int i = 0; i < mp->getNumGeometries(); i++)
        writePoint(mp->getGeometryRef(i)->toPoint(), coords);
}

uint32_t OGRFlatGeobufLayer::writeLineString(OGRLineString *ls, std::vector<double> &coords)
{
    OGRPoint p;
    uint32_t length = 0;
    for (int i = 0; i < ls->getNumPoints(); i++) {
        ls->getPoint(i, &p);
        writePoint(&p, coords);
        length += 2;
    }
    return length;
}

void OGRFlatGeobufLayer::writeMultiLineString(OGRMultiLineString *mls, std::vector<double> &coords, std::vector<uint32_t> &ends)
{
    auto end = 0;
    if (mls->getNumGeometries() > 1)
        for (int i = 0; i < mls->getNumGeometries(); i++)
            ends.push_back(end += writeLineString(mls->getGeometryRef(i)->toLineString(), coords));
    else
        ends.push_back(writeLineString(mls->getGeometryRef(0)->toLineString(), coords));
}

uint32_t OGRFlatGeobufLayer::writePolygon(OGRPolygon *p, std::vector<double> &coords, std::vector<uint32_t> &ends, bool isMulti, uint32_t end)
{
    auto exteriorRing = p->getExteriorRing();
    auto numInteriorRings = p->getNumInteriorRings();
    end += writeLineString(exteriorRing, coords);
    if (numInteriorRings > 0 || isMulti) {
        ends.push_back(end);
        for (int i = 0; i < numInteriorRings; i++)
            ends.push_back(end += writeLineString(p->getInteriorRing(i), coords));
    }
    return end;
}

void OGRFlatGeobufLayer::writeMultiPolygon(
    OGRMultiPolygon *mp,
    std::vector<double> &coords,
    std::vector<uint32_t> &ends,
    std::vector<uint32_t> &endss)
{
    uint32_t end = 0;
    auto isMulti = mp->getNumGeometries() > 1;
    for (int i = 0; i < mp->getNumGeometries(); i++) {
        auto p = mp->getGeometryRef(i)->toPolygon();
        end = writePolygon(p, coords, ends, isMulti, end);
        if (isMulti)
            endss.push_back(p->getNumInteriorRings() + 1);
    }
}

int OGRFlatGeobufLayer::TestCapability(const char *pszCap)
{
    if (EQUAL(pszCap, ODrCCreateDataSource))
        return m_create;
    else if (EQUAL(pszCap, ODsCCreateLayer))
        return m_create;
    else if (EQUAL(pszCap, OLCCreateField))
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


void OGRFlatGeobufLayer::ResetReading()
{
    CPLDebug("FlatGeobuf", "ResetReading");
    m_offset = m_offsetInit;
    m_featuresPos = 0;
    m_processedSpatialIndex = false;
    return;
}
