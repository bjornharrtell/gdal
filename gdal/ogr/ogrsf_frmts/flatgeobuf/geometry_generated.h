// automatically generated by the FlatBuffers compiler, do not modify


#ifndef FLATBUFFERS_GENERATED_GEOMETRY_FLATGEOBUF_H_
#define FLATBUFFERS_GENERATED_GEOMETRY_FLATGEOBUF_H_

#include "flatbuffers/flatbuffers.h"

namespace FlatGeobuf {

struct Geometry;

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
  const size_t index = static_cast<int>(e);
  return EnumNamesGeometryType()[index];
}

struct Geometry FLATBUFFERS_FINAL_CLASS : private flatbuffers::Table {
  enum {
    VT_RING_COUNTS = 4,
    VT_RING_LENGTHS = 6,
    VT_LENGTHS = 8,
    VT_COORDS = 10
  };
  const flatbuffers::Vector<uint32_t> *ring_counts() const {
    return GetPointer<const flatbuffers::Vector<uint32_t> *>(VT_RING_COUNTS);
  }
  const flatbuffers::Vector<uint32_t> *ring_lengths() const {
    return GetPointer<const flatbuffers::Vector<uint32_t> *>(VT_RING_LENGTHS);
  }
  const flatbuffers::Vector<uint32_t> *lengths() const {
    return GetPointer<const flatbuffers::Vector<uint32_t> *>(VT_LENGTHS);
  }
  const flatbuffers::Vector<double> *coords() const {
    return GetPointer<const flatbuffers::Vector<double> *>(VT_COORDS);
  }
  bool Verify(flatbuffers::Verifier &verifier) const {
    return VerifyTableStart(verifier) &&
           VerifyOffset(verifier, VT_RING_COUNTS) &&
           verifier.VerifyVector(ring_counts()) &&
           VerifyOffset(verifier, VT_RING_LENGTHS) &&
           verifier.VerifyVector(ring_lengths()) &&
           VerifyOffset(verifier, VT_LENGTHS) &&
           verifier.VerifyVector(lengths()) &&
           VerifyOffsetRequired(verifier, VT_COORDS) &&
           verifier.VerifyVector(coords()) &&
           verifier.EndTable();
  }
};

struct GeometryBuilder {
  flatbuffers::FlatBufferBuilder &fbb_;
  flatbuffers::uoffset_t start_;
  void add_ring_counts(flatbuffers::Offset<flatbuffers::Vector<uint32_t>> ring_counts) {
    fbb_.AddOffset(Geometry::VT_RING_COUNTS, ring_counts);
  }
  void add_ring_lengths(flatbuffers::Offset<flatbuffers::Vector<uint32_t>> ring_lengths) {
    fbb_.AddOffset(Geometry::VT_RING_LENGTHS, ring_lengths);
  }
  void add_lengths(flatbuffers::Offset<flatbuffers::Vector<uint32_t>> lengths) {
    fbb_.AddOffset(Geometry::VT_LENGTHS, lengths);
  }
  void add_coords(flatbuffers::Offset<flatbuffers::Vector<double>> coords) {
    fbb_.AddOffset(Geometry::VT_COORDS, coords);
  }
  explicit GeometryBuilder(flatbuffers::FlatBufferBuilder &_fbb)
        : fbb_(_fbb) {
    start_ = fbb_.StartTable();
  }
  GeometryBuilder &operator=(const GeometryBuilder &);
  flatbuffers::Offset<Geometry> Finish() {
    const auto end = fbb_.EndTable(start_);
    auto o = flatbuffers::Offset<Geometry>(end);
    fbb_.Required(o, Geometry::VT_COORDS);
    return o;
  }
};

inline flatbuffers::Offset<Geometry> CreateGeometry(
    flatbuffers::FlatBufferBuilder &_fbb,
    flatbuffers::Offset<flatbuffers::Vector<uint32_t>> ring_counts = 0,
    flatbuffers::Offset<flatbuffers::Vector<uint32_t>> ring_lengths = 0,
    flatbuffers::Offset<flatbuffers::Vector<uint32_t>> lengths = 0,
    flatbuffers::Offset<flatbuffers::Vector<double>> coords = 0) {
  GeometryBuilder builder_(_fbb);
  builder_.add_coords(coords);
  builder_.add_lengths(lengths);
  builder_.add_ring_lengths(ring_lengths);
  builder_.add_ring_counts(ring_counts);
  return builder_.Finish();
}

inline flatbuffers::Offset<Geometry> CreateGeometryDirect(
    flatbuffers::FlatBufferBuilder &_fbb,
    const std::vector<uint32_t> *ring_counts = nullptr,
    const std::vector<uint32_t> *ring_lengths = nullptr,
    const std::vector<uint32_t> *lengths = nullptr,
    const std::vector<double> *coords = nullptr) {
  return FlatGeobuf::CreateGeometry(
      _fbb,
      ring_counts ? _fbb.CreateVector<uint32_t>(*ring_counts) : 0,
      ring_lengths ? _fbb.CreateVector<uint32_t>(*ring_lengths) : 0,
      lengths ? _fbb.CreateVector<uint32_t>(*lengths) : 0,
      coords ? _fbb.CreateVector<double>(*coords) : 0);
}

inline const FlatGeobuf::Geometry *GetGeometry(const void *buf) {
  return flatbuffers::GetRoot<FlatGeobuf::Geometry>(buf);
}

inline const FlatGeobuf::Geometry *GetSizePrefixedGeometry(const void *buf) {
  return flatbuffers::GetSizePrefixedRoot<FlatGeobuf::Geometry>(buf);
}

inline bool VerifyGeometryBuffer(
    flatbuffers::Verifier &verifier) {
  return verifier.VerifyBuffer<FlatGeobuf::Geometry>(nullptr);
}

inline bool VerifySizePrefixedGeometryBuffer(
    flatbuffers::Verifier &verifier) {
  return verifier.VerifySizePrefixedBuffer<FlatGeobuf::Geometry>(nullptr);
}

inline void FinishGeometryBuffer(
    flatbuffers::FlatBufferBuilder &fbb,
    flatbuffers::Offset<FlatGeobuf::Geometry> root) {
  fbb.Finish(root);
}

inline void FinishSizePrefixedGeometryBuffer(
    flatbuffers::FlatBufferBuilder &fbb,
    flatbuffers::Offset<FlatGeobuf::Geometry> root) {
  fbb.FinishSizePrefixed(root);
}

}  // namespace FlatGeobuf

#endif  // FLATBUFFERS_GENERATED_GEOMETRY_FLATGEOBUF_H_
