// automatically generated by the FlatBuffers compiler, do not modify


#ifndef FLATBUFFERS_GENERATED_HEADER_FLATGEOBUF_H_
#define FLATBUFFERS_GENERATED_HEADER_FLATGEOBUF_H_

#include "flatbuffers/flatbuffers.h"

namespace FlatGeobuf {

struct Column;

struct Header;

enum class GeometryType : uint8_t {
  Point = 0,
  MultiPoint = 1,
  LineString = 2,
  MultiLineString = 3,
  Polygon = 4,
  MultiPolygon = 5,
  MIN = Point,
  MAX = MultiPolygon
};

inline const GeometryType (&EnumValuesGeometryType())[6] {
  static const GeometryType values[] = {
    GeometryType::Point,
    GeometryType::MultiPoint,
    GeometryType::LineString,
    GeometryType::MultiLineString,
    GeometryType::Polygon,
    GeometryType::MultiPolygon
  };
  return values;
}

inline const char * const *EnumNamesGeometryType() {
  static const char * const names[] = {
    "Point",
    "MultiPoint",
    "LineString",
    "MultiLineString",
    "Polygon",
    "MultiPolygon",
    nullptr
  };
  return names;
}

inline const char *EnumNameGeometryType(GeometryType e) {
  if (e < GeometryType::Point || e > GeometryType::MultiPolygon) return "";
  const size_t index = static_cast<size_t>(e);
  return EnumNamesGeometryType()[index];
}

enum class ColumnType : uint8_t {
  Byte = 0,
  UByte = 1,
  Bool = 2,
  Short = 3,
  UShort = 4,
  Int = 5,
  UInt = 6,
  Long = 7,
  ULong = 8,
  Float = 9,
  Double = 10,
  String = 11,
  Json = 12,
  DateTime = 13,
  MIN = Byte,
  MAX = DateTime
};

inline const ColumnType (&EnumValuesColumnType())[14] {
  static const ColumnType values[] = {
    ColumnType::Byte,
    ColumnType::UByte,
    ColumnType::Bool,
    ColumnType::Short,
    ColumnType::UShort,
    ColumnType::Int,
    ColumnType::UInt,
    ColumnType::Long,
    ColumnType::ULong,
    ColumnType::Float,
    ColumnType::Double,
    ColumnType::String,
    ColumnType::Json,
    ColumnType::DateTime
  };
  return values;
}

inline const char * const *EnumNamesColumnType() {
  static const char * const names[] = {
    "Byte",
    "UByte",
    "Bool",
    "Short",
    "UShort",
    "Int",
    "UInt",
    "Long",
    "ULong",
    "Float",
    "Double",
    "String",
    "Json",
    "DateTime",
    nullptr
  };
  return names;
}

inline const char *EnumNameColumnType(ColumnType e) {
  if (e < ColumnType::Byte || e > ColumnType::DateTime) return "";
  const size_t index = static_cast<size_t>(e);
  return EnumNamesColumnType()[index];
}

struct Column FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
  enum FlatBuffersVTableOffset FLATBUFFERS_VTABLE_UNDERLYING_TYPE {
    VT_NAME = 4,
    VT_TYPE = 6
  };
  const flatbuffers::String *name() const {
    return GetPointer<const flatbuffers::String *>(VT_NAME);
  }
  bool KeyCompareLessThan(const Column *o) const {
    return *name() < *o->name();
  }
  int KeyCompareWithValue(const char *val) const {
    return strcmp(name()->c_str(), val);
  }
  ColumnType type() const {
    return static_cast<ColumnType>(GetField<uint8_t>(VT_TYPE, 0));
  }
  bool Verify(flatbuffers::Verifier &verifier) const {
    return VerifyTableStart(verifier) &&
           VerifyOffsetRequired(verifier, VT_NAME) &&
           verifier.VerifyString(name()) &&
           VerifyField<uint8_t>(verifier, VT_TYPE) &&
           verifier.EndTable();
  }
};

struct ColumnBuilder {
  flatbuffers::FlatBufferBuilder &fbb_;
  flatbuffers::uoffset_t start_;
  void add_name(flatbuffers::Offset<flatbuffers::String> name) {
    fbb_.AddOffset(Column::VT_NAME, name);
  }
  void add_type(ColumnType type) {
    fbb_.AddElement<uint8_t>(Column::VT_TYPE, static_cast<uint8_t>(type), 0);
  }
  explicit ColumnBuilder(flatbuffers::FlatBufferBuilder &_fbb)
        : fbb_(_fbb) {
    start_ = fbb_.StartTable();
  }
  ColumnBuilder &operator=(const ColumnBuilder &);
  flatbuffers::Offset<Column> Finish() {
    const auto end = fbb_.EndTable(start_);
    auto o = flatbuffers::Offset<Column>(end);
    fbb_.Required(o, Column::VT_NAME);
    return o;
  }
};

inline flatbuffers::Offset<Column> CreateColumn(
    flatbuffers::FlatBufferBuilder &_fbb,
    flatbuffers::Offset<flatbuffers::String> name = 0,
    ColumnType type = ColumnType::Byte) {
  ColumnBuilder builder_(_fbb);
  builder_.add_name(name);
  builder_.add_type(type);
  return builder_.Finish();
}

inline flatbuffers::Offset<Column> CreateColumnDirect(
    flatbuffers::FlatBufferBuilder &_fbb,
    const char *name = nullptr,
    ColumnType type = ColumnType::Byte) {
  auto name__ = name ? _fbb.CreateString(name) : 0;
  return FlatGeobuf::CreateColumn(
      _fbb,
      name__,
      type);
}

struct Header FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
  enum FlatBuffersVTableOffset FLATBUFFERS_VTABLE_UNDERLYING_TYPE {
    VT_NAME = 4,
    VT_ENVELOPE = 6,
    VT_GEOMETRY_TYPE = 8,
    VT_DIMENSIONS = 10,
    VT_COLUMNS = 12,
    VT_FEATURES_COUNT = 14,
    VT_FIDS = 16,
    VT_INDEX_NODE_SIZE = 18,
    VT_SRS_CODE = 20,
    VT_SRS_ORG = 22
  };
  const flatbuffers::String *name() const {
    return GetPointer<const flatbuffers::String *>(VT_NAME);
  }
  const flatbuffers::Vector<double> *envelope() const {
    return GetPointer<const flatbuffers::Vector<double> *>(VT_ENVELOPE);
  }
  GeometryType geometry_type() const {
    return static_cast<GeometryType>(GetField<uint8_t>(VT_GEOMETRY_TYPE, 0));
  }
  uint8_t dimensions() const {
    return GetField<uint8_t>(VT_DIMENSIONS, 2);
  }
  const flatbuffers::Vector<flatbuffers::Offset<Column>> *columns() const {
    return GetPointer<const flatbuffers::Vector<flatbuffers::Offset<Column>> *>(VT_COLUMNS);
  }
  uint64_t features_count() const {
    return GetField<uint64_t>(VT_FEATURES_COUNT, 0);
  }
  bool fids() const {
    return GetField<uint8_t>(VT_FIDS, 1) != 0;
  }
  uint16_t index_node_size() const {
    return GetField<uint16_t>(VT_INDEX_NODE_SIZE, 16);
  }
  int32_t srs_code() const {
    return GetField<int32_t>(VT_SRS_CODE, 0);
  }
  const flatbuffers::String *srs_org() const {
    return GetPointer<const flatbuffers::String *>(VT_SRS_ORG);
  }
  bool Verify(flatbuffers::Verifier &verifier) const {
    return VerifyTableStart(verifier) &&
           VerifyOffset(verifier, VT_NAME) &&
           verifier.VerifyString(name()) &&
           VerifyOffset(verifier, VT_ENVELOPE) &&
           verifier.VerifyVector(envelope()) &&
           VerifyField<uint8_t>(verifier, VT_GEOMETRY_TYPE) &&
           VerifyField<uint8_t>(verifier, VT_DIMENSIONS) &&
           VerifyOffset(verifier, VT_COLUMNS) &&
           verifier.VerifyVector(columns()) &&
           verifier.VerifyVectorOfTables(columns()) &&
           VerifyField<uint64_t>(verifier, VT_FEATURES_COUNT) &&
           VerifyField<uint8_t>(verifier, VT_FIDS) &&
           VerifyField<uint16_t>(verifier, VT_INDEX_NODE_SIZE) &&
           VerifyField<int32_t>(verifier, VT_SRS_CODE) &&
           VerifyOffset(verifier, VT_SRS_ORG) &&
           verifier.VerifyString(srs_org()) &&
           verifier.EndTable();
  }
};

struct HeaderBuilder {
  flatbuffers::FlatBufferBuilder &fbb_;
  flatbuffers::uoffset_t start_;
  void add_name(flatbuffers::Offset<flatbuffers::String> name) {
    fbb_.AddOffset(Header::VT_NAME, name);
  }
  void add_envelope(flatbuffers::Offset<flatbuffers::Vector<double>> envelope) {
    fbb_.AddOffset(Header::VT_ENVELOPE, envelope);
  }
  void add_geometry_type(GeometryType geometry_type) {
    fbb_.AddElement<uint8_t>(Header::VT_GEOMETRY_TYPE, static_cast<uint8_t>(geometry_type), 0);
  }
  void add_dimensions(uint8_t dimensions) {
    fbb_.AddElement<uint8_t>(Header::VT_DIMENSIONS, dimensions, 2);
  }
  void add_columns(flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Column>>> columns) {
    fbb_.AddOffset(Header::VT_COLUMNS, columns);
  }
  void add_features_count(uint64_t features_count) {
    fbb_.AddElement<uint64_t>(Header::VT_FEATURES_COUNT, features_count, 0);
  }
  void add_fids(bool fids) {
    fbb_.AddElement<uint8_t>(Header::VT_FIDS, static_cast<uint8_t>(fids), 1);
  }
  void add_index_node_size(uint16_t index_node_size) {
    fbb_.AddElement<uint16_t>(Header::VT_INDEX_NODE_SIZE, index_node_size, 16);
  }
  void add_srs_code(int32_t srs_code) {
    fbb_.AddElement<int32_t>(Header::VT_SRS_CODE, srs_code, 0);
  }
  void add_srs_org(flatbuffers::Offset<flatbuffers::String> srs_org) {
    fbb_.AddOffset(Header::VT_SRS_ORG, srs_org);
  }
  explicit HeaderBuilder(flatbuffers::FlatBufferBuilder &_fbb)
        : fbb_(_fbb) {
    start_ = fbb_.StartTable();
  }
  HeaderBuilder &operator=(const HeaderBuilder &);
  flatbuffers::Offset<Header> Finish() {
    const auto end = fbb_.EndTable(start_);
    auto o = flatbuffers::Offset<Header>(end);
    return o;
  }
};

inline flatbuffers::Offset<Header> CreateHeader(
    flatbuffers::FlatBufferBuilder &_fbb,
    flatbuffers::Offset<flatbuffers::String> name = 0,
    flatbuffers::Offset<flatbuffers::Vector<double>> envelope = 0,
    GeometryType geometry_type = GeometryType::Point,
    uint8_t dimensions = 2,
    flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<Column>>> columns = 0,
    uint64_t features_count = 0,
    bool fids = true,
    uint16_t index_node_size = 16,
    int32_t srs_code = 0,
    flatbuffers::Offset<flatbuffers::String> srs_org = 0) {
  HeaderBuilder builder_(_fbb);
  builder_.add_features_count(features_count);
  builder_.add_srs_org(srs_org);
  builder_.add_srs_code(srs_code);
  builder_.add_columns(columns);
  builder_.add_envelope(envelope);
  builder_.add_name(name);
  builder_.add_index_node_size(index_node_size);
  builder_.add_fids(fids);
  builder_.add_dimensions(dimensions);
  builder_.add_geometry_type(geometry_type);
  return builder_.Finish();
}

inline flatbuffers::Offset<Header> CreateHeaderDirect(
    flatbuffers::FlatBufferBuilder &_fbb,
    const char *name = nullptr,
    const std::vector<double> *envelope = nullptr,
    GeometryType geometry_type = GeometryType::Point,
    uint8_t dimensions = 2,
    const std::vector<flatbuffers::Offset<Column>> *columns = nullptr,
    uint64_t features_count = 0,
    bool fids = true,
    uint16_t index_node_size = 16,
    int32_t srs_code = 0,
    const char *srs_org = nullptr) {
  auto name__ = name ? _fbb.CreateString(name) : 0;
  auto envelope__ = envelope ? _fbb.CreateVector<double>(*envelope) : 0;
  auto columns__ = columns ? _fbb.CreateVector<flatbuffers::Offset<Column>>(*columns) : 0;
  auto srs_org__ = srs_org ? _fbb.CreateString(srs_org) : 0;
  return FlatGeobuf::CreateHeader(
      _fbb,
      name__,
      envelope__,
      geometry_type,
      dimensions,
      columns__,
      features_count,
      fids,
      index_node_size,
      srs_code,
      srs_org__);
}

inline const FlatGeobuf::Header *GetHeader(const void *buf) {
  return flatbuffers::GetRoot<FlatGeobuf::Header>(buf);
}

inline const FlatGeobuf::Header *GetSizePrefixedHeader(const void *buf) {
  return flatbuffers::GetSizePrefixedRoot<FlatGeobuf::Header>(buf);
}

inline bool VerifyHeaderBuffer(
    flatbuffers::Verifier &verifier) {
  return verifier.VerifyBuffer<FlatGeobuf::Header>(nullptr);
}

inline bool VerifySizePrefixedHeaderBuffer(
    flatbuffers::Verifier &verifier) {
  return verifier.VerifySizePrefixedBuffer<FlatGeobuf::Header>(nullptr);
}

inline void FinishHeaderBuffer(
    flatbuffers::FlatBufferBuilder &fbb,
    flatbuffers::Offset<FlatGeobuf::Header> root) {
  fbb.Finish(root);
}

inline void FinishSizePrefixedHeaderBuffer(
    flatbuffers::FlatBufferBuilder &fbb,
    flatbuffers::Offset<FlatGeobuf::Header> root) {
  fbb.FinishSizePrefixed(root);
}

}  // namespace FlatGeobuf

#endif  // FLATBUFFERS_GENERATED_HEADER_FLATGEOBUF_H_
