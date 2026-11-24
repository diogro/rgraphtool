#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "graph_types.h"

namespace {

constexpr uint8_t GT_MAGIC[6] = {0xe2, 0x9b, 0xbe, 0x20, 0x67, 0x74};
constexpr uint8_t GT_VERSION = 0x01;

enum GtKeyType : uint8_t {
    GT_KEY_GRAPH = 0x00,
    GT_KEY_VERTEX = 0x01,
    GT_KEY_EDGE = 0x02
};

enum GtValType : uint8_t {
    GT_VAL_BOOL = 0x00,
    GT_VAL_INT16 = 0x01,
    GT_VAL_INT32 = 0x02,
    GT_VAL_INT64 = 0x03,
    GT_VAL_DOUBLE = 0x04,
    GT_VAL_LONG_DOUBLE = 0x05,
    GT_VAL_STRING = 0x06,
    GT_VAL_VEC_BOOL = 0x07,
    GT_VAL_VEC_INT16 = 0x08,
    GT_VAL_VEC_INT32 = 0x09,
    GT_VAL_VEC_INT64 = 0x0a,
    GT_VAL_VEC_DOUBLE = 0x0b,
    GT_VAL_VEC_LONG_DOUBLE = 0x0c,
    GT_VAL_VEC_STRING = 0x0d,
    GT_VAL_PYTHON_OBJECT = 0x0e
};

inline bool host_is_big_endian() {
    uint16_t x = 1;
    return (*reinterpret_cast<const uint8_t*>(&x) == 0);
}

inline size_t vertex_index_width(uint64_t n) {
    if (n < (uint64_t(1) << 8)) return 1;
    if (n < (uint64_t(1) << 16)) return 2;
    if (n < (uint64_t(1) << 32)) return 4;
    return 8;
}

class GtReader {
public:
    GtReader(const uint8_t* data, size_t size)
        : data_(data), size_(size), pos_(0), swap_(false) {}

    void set_swap(bool swap) { swap_ = swap; }
    size_t remaining() const { return size_ - pos_; }

    void read_bytes(void* dst, size_t n) {
        if (n > remaining()) {
            Rcpp::stop("Unexpected end of .gt file (needed %llu bytes, %llu remaining)",
                       static_cast<unsigned long long>(n),
                       static_cast<unsigned long long>(remaining()));
        }
        if (n > 0) {
            std::memcpy(dst, data_ + pos_, n);
            pos_ += n;
        }
    }

    template <typename T>
    T read_pod() {
        T val{};
        read_bytes(&val, sizeof(T));
        if (swap_ && sizeof(T) > 1) {
            uint8_t* p = reinterpret_cast<uint8_t*>(&val);
            std::reverse(p, p + sizeof(T));
        }
        return val;
    }

    std::string read_string() {
        uint64_t len = read_pod<uint64_t>();
        if (len > remaining()) {
            Rcpp::stop("Invalid string length (%llu) in .gt file",
                       static_cast<unsigned long long>(len));
        }
        std::string s(static_cast<size_t>(len), '\0');
        if (len > 0) {
            read_bytes(&s[0], static_cast<size_t>(len));
        }
        return s;
    }

    Rcpp::RawVector read_raw() {
        uint64_t len = read_pod<uint64_t>();
        if (len > remaining()) {
            Rcpp::stop("Invalid raw/object length (%llu) in .gt file",
                       static_cast<unsigned long long>(len));
        }
        Rcpp::RawVector out(static_cast<R_xlen_t>(len));
        if (len > 0) {
            read_bytes(RAW(out), static_cast<size_t>(len));
        }
        return out;
    }

    uint64_t read_vertex_index(size_t width) {
        switch (width) {
            case 1: return static_cast<uint64_t>(read_pod<uint8_t>());
            case 2: return static_cast<uint64_t>(read_pod<uint16_t>());
            case 4: return static_cast<uint64_t>(read_pod<uint32_t>());
            case 8: return read_pod<uint64_t>();
            default: Rcpp::stop("Internal error: invalid vertex index width");
        }
    }

private:
    const uint8_t* data_;
    size_t size_;
    size_t pos_;
    bool swap_;
};

class GtWriter {
public:
    explicit GtWriter(size_t reserve_bytes = 1024) {
        buf_.reserve(reserve_bytes);
    }

    void write_bytes(const void* src, size_t n) {
        if (n == 0) return;
        const uint8_t* p = reinterpret_cast<const uint8_t*>(src);
        buf_.insert(buf_.end(), p, p + n);
    }

    template <typename T>
    void write_pod(const T& val) {
        write_bytes(&val, sizeof(T));
    }

    void write_string(const std::string& s) {
        uint64_t len = static_cast<uint64_t>(s.size());
        write_pod<uint64_t>(len);
        write_bytes(s.data(), s.size());
    }

    void write_vertex_index(uint64_t v, size_t width) {
        switch (width) {
            case 1: write_pod<uint8_t>(static_cast<uint8_t>(v)); break;
            case 2: write_pod<uint16_t>(static_cast<uint16_t>(v)); break;
            case 4: write_pod<uint32_t>(static_cast<uint32_t>(v)); break;
            case 8: write_pod<uint64_t>(v); break;
            default: Rcpp::stop("Internal error: invalid vertex index width");
        }
    }

    const std::vector<uint8_t>& buffer() const { return buf_; }

private:
    std::vector<uint8_t> buf_;
};

SEXP read_int64_vector(GtReader& r, size_t m) {
    std::vector<int64_t> tmp(m);
    bool fits_int32 = true;
    const int64_t r_int_min = static_cast<int64_t>(std::numeric_limits<int>::min()) + 1;
    const int64_t r_int_max = static_cast<int64_t>(std::numeric_limits<int>::max());
    for (size_t i = 0; i < m; ++i) {
        tmp[i] = r.read_pod<int64_t>();
        if (tmp[i] < r_int_min || tmp[i] > r_int_max) {
            fits_int32 = false;
        }
    }
    if (fits_int32) {
        Rcpp::IntegerVector out(static_cast<R_xlen_t>(m));
        for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = static_cast<int>(tmp[i]);
        return out;
    }
    Rcpp::NumericVector out(static_cast<R_xlen_t>(m));
    for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = static_cast<double>(tmp[i]);
    return out;
}

SEXP read_property_values(GtReader& r, uint8_t vtype, size_t m, bool is_graph_prop) {
    switch (vtype) {
        case GT_VAL_BOOL: {
            Rcpp::LogicalVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = (r.read_pod<uint8_t>() != 0);
            return out;
        }
        case GT_VAL_INT16: {
            Rcpp::IntegerVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = static_cast<int>(r.read_pod<int16_t>());
            return out;
        }
        case GT_VAL_INT32: {
            Rcpp::IntegerVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = static_cast<int>(r.read_pod<int32_t>());
            return out;
        }
        case GT_VAL_INT64: {
            return read_int64_vector(r, m);
        }
        case GT_VAL_DOUBLE: {
            Rcpp::NumericVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = r.read_pod<double>();
            return out;
        }
        case GT_VAL_LONG_DOUBLE: {
            Rcpp::NumericVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = static_cast<double>(r.read_pod<long double>());
            return out;
        }
        case GT_VAL_STRING: {
            Rcpp::CharacterVector out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) out[static_cast<R_xlen_t>(i)] = r.read_string();
            return out;
        }
        case GT_VAL_VEC_BOOL: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining()) Rcpp::stop("Invalid vector<bool> length in .gt file");
                Rcpp::LogicalVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = (r.read_pod<uint8_t>() != 0);
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_INT16: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(int16_t)) Rcpp::stop("Invalid vector<int16_t> length in .gt file");
                Rcpp::IntegerVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = static_cast<int>(r.read_pod<int16_t>());
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_INT32: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(int32_t)) Rcpp::stop("Invalid vector<int32_t> length in .gt file");
                Rcpp::IntegerVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = static_cast<int>(r.read_pod<int32_t>());
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_INT64: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(int64_t)) Rcpp::stop("Invalid vector<int64_t> length in .gt file");
                out[static_cast<R_xlen_t>(i)] = read_int64_vector(r, static_cast<size_t>(len));
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_DOUBLE: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(double)) Rcpp::stop("Invalid vector<double> length in .gt file");
                Rcpp::NumericVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = r.read_pod<double>();
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_LONG_DOUBLE: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(long double)) Rcpp::stop("Invalid vector<long double> length in .gt file");
                Rcpp::NumericVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = static_cast<double>(r.read_pod<long double>());
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_VEC_STRING: {
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                uint64_t len = r.read_pod<uint64_t>();
                if (len > r.remaining() / sizeof(uint64_t)) Rcpp::stop("Invalid vector<string> length in .gt file");
                Rcpp::CharacterVector v(static_cast<R_xlen_t>(len));
                for (uint64_t j = 0; j < len; ++j) v[static_cast<R_xlen_t>(j)] = r.read_string();
                out[static_cast<R_xlen_t>(i)] = v;
            }
            return is_graph_prop ? SEXP(out[0]) : SEXP(out);
        }
        case GT_VAL_PYTHON_OBJECT: {
            if (is_graph_prop) {
                Rcpp::RawVector raw = r.read_raw();
                raw.attr("class") = "gt_python_object";
                return raw;
            }
            Rcpp::List out(static_cast<R_xlen_t>(m));
            for (size_t i = 0; i < m; ++i) {
                Rcpp::RawVector raw = r.read_raw();
                raw.attr("class") = "gt_python_object";
                out[static_cast<R_xlen_t>(i)] = raw;
            }
            out.attr("class") = "gt_python_object_list";
            return out;
        }
        default:
            Rcpp::stop("Unsupported property value type 0x%02x in .gt file",
                       static_cast<unsigned int>(vtype));
    }
}

Rcpp::List parse_gt_buffer(const uint8_t* data, size_t size) {
    GtReader r(data, size);

    uint8_t magic[6];
    r.read_bytes(magic, 6);
    if (std::memcmp(magic, GT_MAGIC, 6) != 0) {
        Rcpp::stop("Invalid .gt file: magic header does not match graph-tool binary format");
    }

    uint8_t version = r.read_pod<uint8_t>();
    if (version != GT_VERSION) {
        Rcpp::stop("Unsupported .gt format version %u (expected %u)",
                   static_cast<unsigned int>(version),
                   static_cast<unsigned int>(GT_VERSION));
    }

    uint8_t endian = r.read_pod<uint8_t>();
    if (endian != 0 && endian != 1) {
        Rcpp::stop("Invalid endianness byte %u in .gt header", static_cast<unsigned int>(endian));
    }
    bool file_is_be = (endian == 1);
    r.set_swap(file_is_be != host_is_big_endian());

    std::string comment = r.read_string();

    uint8_t directed_byte = r.read_pod<uint8_t>();
    if (directed_byte != 0 && directed_byte != 1) {
        Rcpp::stop("Invalid directedness flag %u in .gt file", static_cast<unsigned int>(directed_byte));
    }
    bool directed = (directed_byte == 1);

    uint64_t n_vertices = r.read_pod<uint64_t>();
    if (n_vertices > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
        Rcpp::stop("Graph has too many vertices (%llu) for R's integer vertex indexing",
                   static_cast<unsigned long long>(n_vertices));
    }

    Rcpp::XPtr<Graph> graph(new Graph(directed), true);
    for (uint64_t v = 0; v < n_vertices; ++v) {
        boost::add_vertex(graph->adjacency);
    }

    size_t idx_width = vertex_index_width(n_vertices);
    for (uint64_t v = 0; v < n_vertices; ++v) {
        uint64_t k = r.read_pod<uint64_t>();
        if (k > r.remaining() / idx_width) {
            Rcpp::stop("Invalid neighbor count (%llu) for vertex %llu in .gt file",
                       static_cast<unsigned long long>(k),
                       static_cast<unsigned long long>(v));
        }
        for (uint64_t j = 0; j < k; ++j) {
            uint64_t nbr = r.read_vertex_index(idx_width);
            if (nbr >= n_vertices) {
                Rcpp::stop("Out-of-range neighbor index %llu (num_vertices = %llu) in .gt file",
                           static_cast<unsigned long long>(nbr),
                           static_cast<unsigned long long>(n_vertices));
            }
            boost::add_edge(static_cast<size_t>(v), static_cast<size_t>(nbr), graph->adjacency);
        }
    }

    size_t n_edges = boost::num_edges(graph->adjacency);
    uint64_t n_props = r.read_pod<uint64_t>();

    Rcpp::List graph_props;
    Rcpp::CharacterVector graph_prop_names;
    Rcpp::List vertex_props;
    Rcpp::CharacterVector vertex_prop_names;
    Rcpp::List edge_props_extra;
    Rcpp::CharacterVector edge_prop_extra_names;

    for (uint64_t p = 0; p < n_props; ++p) {
        uint8_t ktype = r.read_pod<uint8_t>();
        std::string name = r.read_string();
        uint8_t vtype = r.read_pod<uint8_t>();

        size_t m = 0;
        if (ktype == GT_KEY_GRAPH) {
            m = 1;
        } else if (ktype == GT_KEY_VERTEX) {
            m = static_cast<size_t>(n_vertices);
        } else if (ktype == GT_KEY_EDGE) {
            m = n_edges;
        } else {
            Rcpp::stop("Invalid property key type 0x%02x in .gt file",
                       static_cast<unsigned int>(ktype));
        }

        SEXP vals = read_property_values(r, vtype, m, ktype == GT_KEY_GRAPH);

        if (ktype == GT_KEY_GRAPH) {
            graph_props.push_back(vals);
            graph_prop_names.push_back(name);
        } else if (ktype == GT_KEY_VERTEX) {
            vertex_props.push_back(vals);
            vertex_prop_names.push_back(name);
        } else {
            // Numeric scalar edge properties go directly into C++ Graph::edge_props
            // so SBM inference and C++ extensions can access them natively.
            if (vtype == GT_VAL_DOUBLE || vtype == GT_VAL_LONG_DOUBLE ||
                vtype == GT_VAL_INT16 || vtype == GT_VAL_INT32 || vtype == GT_VAL_INT64) {
                std::vector<double>& vec = graph->edge_props[name];
                vec.resize(n_edges);
                if (TYPEOF(vals) == REALSXP) {
                    Rcpp::NumericVector nv(vals);
                    for (size_t i = 0; i < n_edges; ++i) vec[i] = nv[static_cast<R_xlen_t>(i)];
                } else if (TYPEOF(vals) == INTSXP) {
                    Rcpp::IntegerVector iv(vals);
                    for (size_t i = 0; i < n_edges; ++i) {
                        int val = iv[static_cast<R_xlen_t>(i)];
                        vec[i] = (val == NA_INTEGER) ? NA_REAL : static_cast<double>(val);
                    }
                }
            } else {
                edge_props_extra.push_back(vals);
                edge_prop_extra_names.push_back(name);
            }
        }
    }

    if (graph_props.size() > 0) graph_props.names() = graph_prop_names;
    if (vertex_props.size() > 0) vertex_props.names() = vertex_prop_names;
    if (edge_props_extra.size() > 0) edge_props_extra.names() = edge_prop_extra_names;

    return Rcpp::List::create(
        Rcpp::Named("pointer") = graph,
        Rcpp::Named("directed") = directed,
        Rcpp::Named("comment") = comment,
        Rcpp::Named("graph_props") = graph_props,
        Rcpp::Named("vertex_props") = vertex_props,
        Rcpp::Named("edge_props_extra") = edge_props_extra
    );
}

void write_graph_property(GtWriter& w, const std::string& name, SEXP val) {
    w.write_pod<uint8_t>(GT_KEY_GRAPH);
    w.write_string(name);

    if (TYPEOF(val) == RAWSXP) {
        w.write_pod<uint8_t>(GT_VAL_PYTHON_OBJECT);
        Rcpp::RawVector rv(val);
        w.write_pod<uint64_t>(static_cast<uint64_t>(rv.size()));
        if (rv.size() > 0) w.write_bytes(RAW(rv), static_cast<size_t>(rv.size()));
        return;
    }

    R_xlen_t n = Rf_xlength(val);
    if (TYPEOF(val) == LGLSXP) {
        Rcpp::LogicalVector lv(val);
        if (n == 1) {
            w.write_pod<uint8_t>(GT_VAL_BOOL);
            w.write_pod<uint8_t>(lv[0] == TRUE ? 1 : 0);
        } else {
            w.write_pod<uint8_t>(GT_VAL_VEC_BOOL);
            w.write_pod<uint64_t>(static_cast<uint64_t>(n));
            for (R_xlen_t i = 0; i < n; ++i) w.write_pod<uint8_t>(lv[i] == TRUE ? 1 : 0);
        }
    } else if (TYPEOF(val) == INTSXP) {
        Rcpp::IntegerVector iv(val);
        if (n == 1) {
            w.write_pod<uint8_t>(GT_VAL_INT32);
            w.write_pod<int32_t>(static_cast<int32_t>(iv[0]));
        } else {
            w.write_pod<uint8_t>(GT_VAL_VEC_INT32);
            w.write_pod<uint64_t>(static_cast<uint64_t>(n));
            for (R_xlen_t i = 0; i < n; ++i) w.write_pod<int32_t>(static_cast<int32_t>(iv[i]));
        }
    } else if (TYPEOF(val) == REALSXP) {
        Rcpp::NumericVector nv(val);
        if (n == 1) {
            w.write_pod<uint8_t>(GT_VAL_DOUBLE);
            w.write_pod<double>(nv[0]);
        } else {
            w.write_pod<uint8_t>(GT_VAL_VEC_DOUBLE);
            w.write_pod<uint64_t>(static_cast<uint64_t>(n));
            for (R_xlen_t i = 0; i < n; ++i) w.write_pod<double>(nv[i]);
        }
    } else if (TYPEOF(val) == STRSXP) {
        Rcpp::CharacterVector cv(val);
        if (n == 1) {
            w.write_pod<uint8_t>(GT_VAL_STRING);
            w.write_string(cv[0] == NA_STRING ? "" : std::string(CHAR(cv[0])));
        } else {
            w.write_pod<uint8_t>(GT_VAL_VEC_STRING);
            w.write_pod<uint64_t>(static_cast<uint64_t>(n));
            for (R_xlen_t i = 0; i < n; ++i) {
                w.write_string(cv[i] == NA_STRING ? "" : std::string(CHAR(cv[i])));
            }
        }
    } else {
        Rcpp::stop("Unsupported R type for graph property '%s' when writing .gt", name.c_str());
    }
}

void write_seq_property(GtWriter& w, uint8_t ktype, const std::string& name,
                        SEXP val, size_t expected_len) {
    w.write_pod<uint8_t>(ktype);
    w.write_string(name);

    // Handle matrices (nrow == expected_len) as vector<T> per vertex/edge
    if (Rf_isMatrix(val)) {
        int nr = Rf_nrows(val);
        int nc = Rf_ncols(val);
        if (static_cast<size_t>(nr) != expected_len) {
            Rcpp::stop("Matrix property '%s' must have %llu rows (got %d)",
                       name.c_str(), static_cast<unsigned long long>(expected_len), nr);
        }
        if (TYPEOF(val) == REALSXP) {
            Rcpp::NumericMatrix mat(val);
            w.write_pod<uint8_t>(GT_VAL_VEC_DOUBLE);
            for (int i = 0; i < nr; ++i) {
                w.write_pod<uint64_t>(static_cast<uint64_t>(nc));
                for (int j = 0; j < nc; ++j) w.write_pod<double>(mat(i, j));
            }
            return;
        } else if (TYPEOF(val) == INTSXP) {
            Rcpp::IntegerMatrix mat(val);
            w.write_pod<uint8_t>(GT_VAL_VEC_INT32);
            for (int i = 0; i < nr; ++i) {
                w.write_pod<uint64_t>(static_cast<uint64_t>(nc));
                for (int j = 0; j < nc; ++j) w.write_pod<int32_t>(static_cast<int32_t>(mat(i, j)));
            }
            return;
        } else if (TYPEOF(val) == LGLSXP) {
            Rcpp::LogicalMatrix mat(val);
            w.write_pod<uint8_t>(GT_VAL_VEC_BOOL);
            for (int i = 0; i < nr; ++i) {
                w.write_pod<uint64_t>(static_cast<uint64_t>(nc));
                for (int j = 0; j < nc; ++j) w.write_pod<uint8_t>(mat(i, j) == TRUE ? 1 : 0);
            }
            return;
        } else if (TYPEOF(val) == STRSXP) {
            Rcpp::CharacterMatrix mat(val);
            w.write_pod<uint8_t>(GT_VAL_VEC_STRING);
            for (int i = 0; i < nr; ++i) {
                w.write_pod<uint64_t>(static_cast<uint64_t>(nc));
                for (int j = 0; j < nc; ++j) {
                    SEXP el = mat(i, j);
                    w.write_string(el == NA_STRING ? "" : std::string(CHAR(el)));
                }
            }
            return;
        }
    }

    R_xlen_t len = Rf_xlength(val);
    if (static_cast<size_t>(len) != expected_len) {
        Rcpp::stop("Property '%s' length (%lld) does not match expected length (%llu)",
                   name.c_str(), static_cast<long long>(len),
                   static_cast<unsigned long long>(expected_len));
    }

    if (TYPEOF(val) == LGLSXP) {
        Rcpp::LogicalVector lv(val);
        w.write_pod<uint8_t>(GT_VAL_BOOL);
        for (R_xlen_t i = 0; i < len; ++i) w.write_pod<uint8_t>(lv[i] == TRUE ? 1 : 0);
    } else if (TYPEOF(val) == INTSXP) {
        Rcpp::IntegerVector iv(val);
        w.write_pod<uint8_t>(GT_VAL_INT32);
        for (R_xlen_t i = 0; i < len; ++i) w.write_pod<int32_t>(static_cast<int32_t>(iv[i]));
    } else if (TYPEOF(val) == REALSXP) {
        Rcpp::NumericVector nv(val);
        w.write_pod<uint8_t>(GT_VAL_DOUBLE);
        for (R_xlen_t i = 0; i < len; ++i) w.write_pod<double>(nv[i]);
    } else if (TYPEOF(val) == STRSXP) {
        Rcpp::CharacterVector cv(val);
        w.write_pod<uint8_t>(GT_VAL_STRING);
        for (R_xlen_t i = 0; i < len; ++i) {
            w.write_string(cv[i] == NA_STRING ? "" : std::string(CHAR(cv[i])));
        }
    } else if (TYPEOF(val) == VECSXP) {
        Rcpp::List lst(val);
        // Determine vector element type from non-empty entries
        bool has_raw = false, has_str = false, has_real = false, has_int = false, has_lgl = false;
        for (R_xlen_t i = 0; i < len; ++i) {
            SEXP el = lst[i];
            if (el == R_NilValue) continue;
            switch (TYPEOF(el)) {
                case RAWSXP: has_raw = true; break;
                case STRSXP: has_str = true; break;
                case REALSXP: has_real = true; break;
                case INTSXP: has_int = true; break;
                case LGLSXP: has_lgl = true; break;
                default:
                    Rcpp::stop("Unsupported element type in list property '%s'", name.c_str());
            }
        }
        if (has_raw && !has_str && !has_real && !has_int && !has_lgl) {
            w.write_pod<uint8_t>(GT_VAL_PYTHON_OBJECT);
            for (R_xlen_t i = 0; i < len; ++i) {
                SEXP el = lst[i];
                if (el == R_NilValue) {
                    w.write_pod<uint64_t>(0);
                } else {
                    Rcpp::RawVector rv(el);
                    w.write_pod<uint64_t>(static_cast<uint64_t>(rv.size()));
                    if (rv.size() > 0) w.write_bytes(RAW(rv), static_cast<size_t>(rv.size()));
                }
            }
        } else if (has_str) {
            w.write_pod<uint8_t>(GT_VAL_VEC_STRING);
            for (R_xlen_t i = 0; i < len; ++i) {
                SEXP el = lst[i];
                if (el == R_NilValue) {
                    w.write_pod<uint64_t>(0);
                } else {
                    Rcpp::CharacterVector cv = Rcpp::as<Rcpp::CharacterVector>(el);
                    w.write_pod<uint64_t>(static_cast<uint64_t>(cv.size()));
                    for (R_xlen_t j = 0; j < cv.size(); ++j) {
                        w.write_string(cv[j] == NA_STRING ? "" : std::string(CHAR(cv[j])));
                    }
                }
            }
        } else if (has_real || (!has_int && !has_lgl)) {
            w.write_pod<uint8_t>(GT_VAL_VEC_DOUBLE);
            for (R_xlen_t i = 0; i < len; ++i) {
                SEXP el = lst[i];
                if (el == R_NilValue) {
                    w.write_pod<uint64_t>(0);
                } else {
                    Rcpp::NumericVector nv = Rcpp::as<Rcpp::NumericVector>(el);
                    w.write_pod<uint64_t>(static_cast<uint64_t>(nv.size()));
                    for (R_xlen_t j = 0; j < nv.size(); ++j) w.write_pod<double>(nv[j]);
                }
            }
        } else if (has_int) {
            w.write_pod<uint8_t>(GT_VAL_VEC_INT32);
            for (R_xlen_t i = 0; i < len; ++i) {
                SEXP el = lst[i];
                if (el == R_NilValue) {
                    w.write_pod<uint64_t>(0);
                } else {
                    Rcpp::IntegerVector iv = Rcpp::as<Rcpp::IntegerVector>(el);
                    w.write_pod<uint64_t>(static_cast<uint64_t>(iv.size()));
                    for (R_xlen_t j = 0; j < iv.size(); ++j) w.write_pod<int32_t>(static_cast<int32_t>(iv[j]));
                }
            }
        } else {
            w.write_pod<uint8_t>(GT_VAL_VEC_BOOL);
            for (R_xlen_t i = 0; i < len; ++i) {
                SEXP el = lst[i];
                if (el == R_NilValue) {
                    w.write_pod<uint64_t>(0);
                } else {
                    Rcpp::LogicalVector lv = Rcpp::as<Rcpp::LogicalVector>(el);
                    w.write_pod<uint64_t>(static_cast<uint64_t>(lv.size()));
                    for (R_xlen_t j = 0; j < lv.size(); ++j) w.write_pod<uint8_t>(lv[j] == TRUE ? 1 : 0);
                }
            }
        }
    } else {
        Rcpp::stop("Unsupported R type for property '%s' when writing .gt", name.c_str());
    }
}

std::vector<uint8_t> serialize_gt_buffer(const Graph& graph,
                                         const std::string& comment,
                                         const Rcpp::List& graph_props,
                                         const Rcpp::List& vertex_props,
                                         const Rcpp::List& edge_props_extra) {
    const size_t n_vertices = boost::num_vertices(graph.adjacency);
    const size_t n_edges = boost::num_edges(graph.adjacency);
    const size_t idx_width = vertex_index_width(static_cast<uint64_t>(n_vertices));

    GtWriter w(256 + n_vertices * 16 + n_edges * idx_width);

    // 1. Header
    w.write_bytes(GT_MAGIC, 6);
    w.write_pod<uint8_t>(GT_VERSION);
    w.write_pod<uint8_t>(host_is_big_endian() ? 1 : 0);
    w.write_string(comment);

    // 2. Adjacency
    w.write_pod<uint8_t>(graph.directed ? 1 : 0);
    w.write_pod<uint64_t>(static_cast<uint64_t>(n_vertices));

    // Record edge insertion indices in exact (v, out_edges(v)) emission order
    std::vector<size_t> edge_idx_order;
    edge_idx_order.reserve(n_edges);
    for (size_t v = 0; v < n_vertices; ++v) {
        uint64_t k = static_cast<uint64_t>(boost::out_degree(v, graph.adjacency));
        w.write_pod<uint64_t>(k);
        auto range = boost::out_edges(v, graph.adjacency);
        for (auto it = range.first; it != range.second; ++it) {
            w.write_vertex_index(static_cast<uint64_t>(boost::target(*it, graph.adjacency)), idx_width);
            edge_idx_order.push_back(it->idx);
        }
    }

    // 3. Property maps
    const uint64_t n_total_props =
        static_cast<uint64_t>(graph_props.size()) +
        static_cast<uint64_t>(vertex_props.size()) +
        static_cast<uint64_t>(graph.edge_props.size()) +
        static_cast<uint64_t>(edge_props_extra.size());
    w.write_pod<uint64_t>(n_total_props);

    if (graph_props.size() > 0) {
        Rcpp::CharacterVector gnames = graph_props.names();
        for (R_xlen_t i = 0; i < graph_props.size(); ++i) {
            write_graph_property(w, std::string(CHAR(gnames[i])), graph_props[i]);
        }
    }

    if (vertex_props.size() > 0) {
        Rcpp::CharacterVector vnames = vertex_props.names();
        for (R_xlen_t i = 0; i < vertex_props.size(); ++i) {
            write_seq_property(w, GT_KEY_VERTEX, std::string(CHAR(vnames[i])),
                               vertex_props[i], n_vertices);
        }
    }

    // C++ numeric edge properties
    for (const auto& kv : graph.edge_props) {
        w.write_pod<uint8_t>(GT_KEY_EDGE);
        w.write_string(kv.first);
        w.write_pod<uint8_t>(GT_VAL_DOUBLE);
        const std::vector<double>& vec = kv.second;
        for (size_t i = 0; i < edge_idx_order.size(); ++i) {
            size_t eidx = edge_idx_order[i];
            double val = (eidx < vec.size()) ? vec[eidx] : NA_REAL;
            w.write_pod<double>(val);
        }
    }

    // Extra (non-numeric) R edge properties (aligned with boost::edges order)
    if (edge_props_extra.size() > 0) {
        Rcpp::CharacterVector enames = edge_props_extra.names();
        for (R_xlen_t i = 0; i < edge_props_extra.size(); ++i) {
            write_seq_property(w, GT_KEY_EDGE, std::string(CHAR(enames[i])),
                               edge_props_extra[i], n_edges);
        }
    }

    return w.buffer();
}

} // namespace

// [[Rcpp::export]]
Rcpp::List cpp_read_gt_raw(Rcpp::RawVector bytes) {
    const uint8_t* ptr = bytes.size() > 0 ? reinterpret_cast<const uint8_t*>(RAW(bytes)) : nullptr;
    return parse_gt_buffer(ptr, static_cast<size_t>(bytes.size()));
}

// [[Rcpp::export]]
Rcpp::List cpp_read_gt_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) {
        Rcpp::stop("Could not open .gt file for reading: %s", path.c_str());
    }
    std::streamsize sz = in.tellg();
    if (sz < 0) {
        Rcpp::stop("Could not determine size of .gt file: %s", path.c_str());
    }
    in.seekg(0, std::ios::beg);
    std::vector<uint8_t> buf(static_cast<size_t>(sz));
    if (sz > 0 && !in.read(reinterpret_cast<char*>(buf.data()), sz)) {
        Rcpp::stop("Error reading .gt file: %s", path.c_str());
    }
    return parse_gt_buffer(buf.data(), buf.size());
}

// [[Rcpp::export]]
Rcpp::RawVector cpp_write_gt_raw(SEXP pointer,
                                 const std::string& comment,
                                 Rcpp::List graph_props,
                                 Rcpp::List vertex_props,
                                 Rcpp::List edge_props_extra) {
    const Graph& graph = get_graph(pointer);
    std::vector<uint8_t> buf = serialize_gt_buffer(graph, comment, graph_props,
                                                   vertex_props, edge_props_extra);
    Rcpp::RawVector out(static_cast<R_xlen_t>(buf.size()));
    if (!buf.empty()) {
        std::memcpy(RAW(out), buf.data(), buf.size());
    }
    return out;
}

// [[Rcpp::export]]
void cpp_write_gt_file(SEXP pointer,
                       const std::string& path,
                       const std::string& comment,
                       Rcpp::List graph_props,
                       Rcpp::List vertex_props,
                       Rcpp::List edge_props_extra) {
    const Graph& graph = get_graph(pointer);
    std::vector<uint8_t> buf = serialize_gt_buffer(graph, comment, graph_props,
                                                   vertex_props, edge_props_extra);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        Rcpp::stop("Could not open .gt file for writing: %s", path.c_str());
    }
    if (!buf.empty()) {
        out.write(reinterpret_cast<const char*>(buf.data()),
                  static_cast<std::streamsize>(buf.size()));
    }
    if (!out) {
        Rcpp::stop("Error writing .gt file: %s", path.c_str());
    }
}
