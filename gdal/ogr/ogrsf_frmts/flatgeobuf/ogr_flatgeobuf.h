#ifndef OGR_FLATGEOBUF_H_INCLUDED
#define OGR_FLATGEOBUF_H_INCLUDED

#include "ogrsf_frmts.h"
#include "ogr_p.h"

#include "flatgeobuf_generated.h"
#include "packedrtree.h"

using namespace FlatGeobuf;

class OGRFlatGeobufDataset;

struct FeatureItem : Item {
    flatbuffers::DetachedBuffer buf;
    uint8_t *data;
    size_t size;
};

class OGRFlatGeobufLayer : public OGRLayer
{
    private:
        const char *m_pszFilename = nullptr;
        VSILFILE *m_poFp = nullptr;
        const char *m_pszLayerName = nullptr;

        const Header *m_poHeader = nullptr;
        GeometryType m_geometryType;
        uint8_t m_dimensions = 2;
        uint64_t m_featuresCount = 0;

        OGRFeatureDefn *m_poFeatureDefn = nullptr;
        OGRSpatialReference *m_poSRS = nullptr;

        // iterator variables
        uint64_t m_featuresPos = 0;
        uint64_t m_featuresSize = 0;
        uint64_t m_offset = 0;
        uint64_t m_offsetInit = 0;
        uint64_t *m_featureOffsets = nullptr;
        std::vector<uint64_t> m_foundFeatureIndices;
        bool m_processedSpatialIndex = false;

        bool m_create = false;

        // creation buffers
        std::vector<Item *> m_featureItems;
        GByte *m_featureBuf = nullptr;
        uint32_t m_featureSize;
        uint32_t m_featureBufSize = 0;
        size_t m_padfSize = 0;
        double *m_padfX = nullptr;
        double *m_padfY = nullptr;
        double *m_padfZ = nullptr;
        double *m_padfM = nullptr;

        // deserialize
        void ensurePadfBuffers(size_t count, uint8_t dimensions);
        OGRLineString *readLineString(const double *coords, uint32_t coordsLength, uint8_t dimensions, uint32_t offset = 0);
        OGRLinearRing *readLinearRing(const double *coords, uint32_t coordsLength, uint8_t dimensions, uint32_t offset = 0);
        OGRPolygon *readPolygon(const double *coords, uint32_t coordsLength, const flatbuffers::Vector<uint32_t> *ringLengths, uint8_t dimensions);
        OGRGeometry *readGeometry(const Geometry* geometry, uint8_t dimensions);
        ColumnType toColumnType(OGRFieldType fieldType, OGRFieldSubType subType);
        OGRFieldType toOGRFieldType(ColumnType type);
        const std::vector<flatbuffers::Offset<Column>> writeColumns(flatbuffers::FlatBufferBuilder &fbb);
        void readColumns();
        void processSpatialIndex();

        // serialize
        void writePoint(OGRPoint *p, std::vector<double> &coords);
        uint32_t writeLineString(OGRLineString *ls, std::vector<double> &coords);
        void writePolygon(OGRPolygon *p, std::vector<double> &coords, std::vector<uint32_t> &ringLengths);
        flatbuffers::Offset<Geometry> writeGeometry(flatbuffers::FlatBufferBuilder& fbb, OGRGeometry* ogrGeometry);
    public:
        OGRFlatGeobufLayer(const Header*, const char* pszFilename, uint64_t offset);
        OGRFlatGeobufLayer(const char *pszLayerName, const char *pszFilename, OGRSpatialReference *poSpatialRef, OGRwkbGeometryType eGType);
        virtual ~OGRFlatGeobufLayer();

        virtual OGRFeature *GetFeature(GIntBig nFeatureId) override;
        virtual OGRFeature *GetNextFeature() override;
        virtual OGRErr CreateField(OGRFieldDefn *poField, int bApproxOK = true) override;
        virtual OGRErr ICreateFeature(OGRFeature *poFeature) override;
        virtual int TestCapability(const char *) override;

        virtual void ResetReading() override;
        virtual OGRFeatureDefn *GetLayerDefn() override { return m_poFeatureDefn; }
        virtual GIntBig GetFeatureCount(int bForce) override;
};

class OGRFlatGeobufDataset final: public GDALDataset
{
    private:
        const char *m_pszName = nullptr;
        const char *m_pszFilename = nullptr;
        std::vector<std::unique_ptr<OGRLayer>> m_apoLayers;
        bool m_create = false;
    public:
        explicit OGRFlatGeobufDataset();
        explicit OGRFlatGeobufDataset(const char *pszName);
        ~OGRFlatGeobufDataset();

        static GDALDataset *Open(GDALOpenInfo*);
        static GDALDataset *Create( const char *pszName,
                                        CPL_UNUSED int nBands,
                                        CPL_UNUSED int nXSize,
                                        CPL_UNUSED int nYSize,
                                        CPL_UNUSED GDALDataType eDT,
                                        char **papszOptions );
        virtual OGRLayer *GetLayer( int ) override;
        int TestCapability( const char * pszCap ) override;
        virtual OGRLayer *ICreateLayer( const char *pszName,
                                     OGRSpatialReference *poSpatialRef = nullptr,
                                     OGRwkbGeometryType eGType = wkbUnknown,
                                     char ** papszOptions = nullptr ) override;

        virtual int GetLayerCount() override { return static_cast<int>(m_apoLayers.size()); }

        static GeometryType toGeometryType(OGRwkbGeometryType eGType);
        static OGRwkbGeometryType toOGRwkbGeometryType(GeometryType geometryType);
};

#endif /* ndef OGR_FLATGEOBUF_H_INCLUDED */
