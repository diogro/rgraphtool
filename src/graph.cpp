#include <Rcpp.h>
#include <limits>
#include <map>
#include <cmath>
#include "graph_types.h"

// [[Rcpp::export]]
SEXP cpp_graph_create(int n, bool directed) {
    Rcpp::XPtr<Graph> graph(new Graph(directed), true);
    for (int i = 0; i < n; ++i)
        boost::add_vertex(graph->adjacency);
    return graph;
}

// [[Rcpp::export]]
void cpp_add_vertices(SEXP pointer, int n) {
    Graph& graph = get_graph(pointer);
    for (int i = 0; i < n; ++i)
        boost::add_vertex(graph.adjacency);
}

// [[Rcpp::export]]
void cpp_add_edges(SEXP pointer, Rcpp::IntegerMatrix endpoints) {
    Graph& graph = get_graph(pointer);
    size_t old_e = boost::num_edges(graph.adjacency);
    for (int i = 0; i < endpoints.nrow(); ++i)
        boost::add_edge(static_cast<size_t>(endpoints(i, 0) - 1),
                        static_cast<size_t>(endpoints(i, 1) - 1),
                        graph.adjacency);
    size_t new_e = boost::num_edges(graph.adjacency);
    if (new_e > old_e && !graph.edge_props.empty()) {
        for (auto& kv : graph.edge_props)
            kv.second.resize(new_e, NA_REAL);
    }
}

// [[Rcpp::export]]
void cpp_set_edge_property(SEXP pointer, const std::string& name,
                           Rcpp::NumericVector values, bool by_insertion_order,
                           int offset) {
    Graph& graph = get_graph(pointer);
    size_t count = boost::num_edges(graph.adjacency);
    if (by_insertion_order) {
        if (offset < 0 || static_cast<size_t>(offset + values.size()) > count)
            Rcpp::stop("Edge property length does not match number of edges");
        auto& vec = graph.edge_props[name];
        if (vec.size() != count)
            vec.resize(count, NA_REAL);
        for (R_xlen_t i = 0; i < values.size(); ++i)
            vec[static_cast<size_t>(offset + i)] = values[i];
    } else {
        if (static_cast<size_t>(values.size()) != count)
            Rcpp::stop("Edge property length must equal graph_num_edges(graph)");
        std::vector<double>& vec = graph.edge_props[name];
        vec.assign(count, NA_REAL);
        auto range = boost::edges(graph.adjacency);
        R_xlen_t pos = 0;
        for (auto it = range.first; it != range.second; ++it, ++pos) {
            if (it->idx < count)
                vec[it->idx] = values[pos];
        }
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_get_edge_property(SEXP pointer, const std::string& name) {
    Graph& graph = get_graph(pointer);
    auto found = graph.edge_props.find(name);
    if (found == graph.edge_props.end())
        Rcpp::stop("Edge property '%s' does not exist on graph", name.c_str());
    const std::vector<double>& vec = found->second;
    size_t count = boost::num_edges(graph.adjacency);
    Rcpp::NumericVector out(static_cast<R_xlen_t>(count));
    auto range = boost::edges(graph.adjacency);
    R_xlen_t pos = 0;
    for (auto it = range.first; it != range.second; ++it, ++pos)
        out[pos] = (it->idx < vec.size()) ? vec[it->idx] : NA_REAL;
    return out;
}

// [[Rcpp::export]]
Rcpp::CharacterVector cpp_edge_property_names(SEXP pointer) {
    Graph& graph = get_graph(pointer);
    Rcpp::CharacterVector out(static_cast<R_xlen_t>(graph.edge_props.size()));
    R_xlen_t i = 0;
    for (const auto& kv : graph.edge_props)
        out[i++] = kv.first;
    return out;
}

// [[Rcpp::export]]
double cpp_num_vertices(SEXP pointer) {
    return static_cast<double>(boost::num_vertices(get_graph(pointer).adjacency));
}

// [[Rcpp::export]]
double cpp_num_edges(SEXP pointer) {
    return static_cast<double>(boost::num_edges(get_graph(pointer).adjacency));
}

// [[Rcpp::export]]
Rcpp::IntegerVector cpp_edge_order(SEXP pointer) {
    Graph& graph = get_graph(pointer);
    size_t count = boost::num_edges(graph.adjacency);
    if (count > static_cast<size_t>(std::numeric_limits<int>::max()))
        Rcpp::stop("Too many edges to return edge order");
    Rcpp::IntegerVector out(static_cast<R_xlen_t>(count));
    auto range = boost::edges(graph.adjacency);
    R_xlen_t pos = 0;
    for (auto it = range.first; it != range.second; ++it, ++pos) {
        out[pos] = static_cast<int>(it->idx) + 1;
    }
    return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_edges(SEXP pointer) {
    Graph& graph = get_graph(pointer);
    size_t count = boost::num_edges(graph.adjacency);
    if (count > static_cast<size_t>(std::numeric_limits<int>::max()))
        Rcpp::stop("Too many edges to return as a matrix");
    Rcpp::IntegerMatrix result(static_cast<int>(count), 2);
    auto range = boost::edges(graph.adjacency);
    int i = 0;
    for (auto it = range.first; it != range.second; ++it, ++i) {
        result(i, 0) = static_cast<int>(it->s + 1);
        result(i, 1) = static_cast<int>(it->t + 1);
    }
    Rcpp::colnames(result) = Rcpp::CharacterVector::create("source", "target");
    return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector cpp_degree(SEXP pointer, const std::string& mode) {
    Graph& graph = get_graph(pointer);
    if (!graph.directed && mode != "total")
        Rcpp::stop("Undirected graphs support only total degree");
    size_t n = boost::num_vertices(graph.adjacency);
    Rcpp::NumericVector result(n);
    for (size_t i = 0; i < n; ++i) {
        size_t in = boost::in_degree(i, graph.adjacency);
        size_t out = boost::out_degree(i, graph.adjacency);
        result[i] = static_cast<double>(mode == "in" ? in : mode == "out" ? out : in + out);
    }
    return result;
}

// [[Rcpp::export]]
Rcpp::List cpp_block_statistics(SEXP pointer, Rcpp::IntegerVector blocks, int nblocks) {
    Graph& graph = get_graph(pointer);
    const Adjacency& adjacency = graph.adjacency;
    Rcpp::NumericMatrix counts(nblocks, nblocks);
    Rcpp::NumericVector sizes(nblocks);
    for (R_xlen_t i = 0; i < blocks.size(); ++i)
        ++sizes[blocks[i] - 1];

    std::map<std::pair<size_t, size_t>, size_t> multiplicities;
    auto range = boost::edges(adjacency);
    for (auto it = range.first; it != range.second; ++it) {
        size_t source = it->s, target = it->t;
        int r = blocks[source] - 1, s = blocks[target] - 1;
        ++counts(r, s);
        if (!graph.directed)
            ++counts(s, r);
        if (!graph.directed && source > target)
            std::swap(source, target);
        ++multiplicities[std::make_pair(source, target)];
    }

    double edge_factor = 0;
    for (const auto& entry : multiplicities) {
        edge_factor += std::lgamma(static_cast<double>(entry.second) + 1);
        if (!graph.directed && entry.first.first == entry.first.second)
            edge_factor += entry.second * std::log(2.0);
    }
    double degree_factor = 0;
    Rcpp::NumericVector in_degrees(boost::num_vertices(adjacency));
    Rcpp::NumericVector out_degrees(boost::num_vertices(adjacency));
    for (size_t i = 0; i < boost::num_vertices(adjacency); ++i) {
        size_t in = boost::in_degree(i, adjacency);
        size_t out = boost::out_degree(i, adjacency);
        in_degrees[i] = graph.directed ? in : in + out;
        out_degrees[i] = graph.directed ? out : in + out;
        if (graph.directed)
            degree_factor += std::lgamma(static_cast<double>(in) + 1) +
                             std::lgamma(static_cast<double>(out) + 1);
        else
            degree_factor += std::lgamma(static_cast<double>(in + out) + 1);
    }
    return Rcpp::List::create(Rcpp::Named("counts") = counts,
                              Rcpp::Named("sizes") = sizes,
                              Rcpp::Named("directed") = graph.directed,
                              Rcpp::Named("edge_factor") = edge_factor,
                              Rcpp::Named("degree_factor") = degree_factor,
                              Rcpp::Named("in_degrees") = in_degrees,
                              Rcpp::Named("out_degrees") = out_degrees);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_log_partition_table(int max_degree, int max_vertices) {
    const double cells = (static_cast<double>(max_degree) + 1) *
                         (static_cast<double>(max_vertices) + 1);
    if (cells > 10000000)
        Rcpp::stop("Degree partition table exceeds 10 million entries");
    Rcpp::NumericMatrix table(max_degree + 1, max_vertices + 1);
    std::fill(table.begin(), table.end(), R_NegInf);
    for (int k = 0; k <= max_vertices; ++k)
        table(0, k) = 0;
    for (int k = 1; k <= max_vertices; ++k) {
        Rcpp::checkUserInterrupt();
        for (int m = 1; m <= max_degree; ++m) {
            double a = table(m, k - 1);
            double b = m >= k ? table(m - k, k) : R_NegInf;
            if (a == R_NegInf)
                table(m, k) = b;
            else if (b == R_NegInf)
                table(m, k) = a;
            else {
                double high = std::max(a, b);
                table(m, k) = high + std::log1p(std::exp(std::min(a, b) - high));
            }
        }
    }
    return table;
}
