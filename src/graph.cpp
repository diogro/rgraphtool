#include <Rcpp.h>
#include <limits>
#include "graph_adjacency.h"

typedef boost::adj_list<size_t> Adjacency;

struct Graph {
    Adjacency adjacency;
    bool directed;
    explicit Graph(bool directed_) : directed(directed_) {}
};

static Graph& get_graph(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL)
        Rcpp::stop("Graph pointer is invalid; graphs cannot be restored from serialized R objects");
    return *Rcpp::XPtr<Graph>(pointer);
}

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
    for (int i = 0; i < endpoints.nrow(); ++i)
        boost::add_edge(static_cast<size_t>(endpoints(i, 0) - 1),
                        static_cast<size_t>(endpoints(i, 1) - 1),
                        graph.adjacency);
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
