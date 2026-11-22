// graph-tool-R: shared graph type used by the native interface.
// Copyright (C) 2026 Diogo Melo. LGPL-3.0-or-later.
#ifndef GRAPHTOOL_R_GRAPH_TYPES_H
#define GRAPHTOOL_R_GRAPH_TYPES_H

#include <Rcpp.h>
#include "vendor/graph_adjacency.hh"

typedef boost::adj_list<size_t> Adjacency;

struct Graph {
    Adjacency adjacency;
    bool directed;
    explicit Graph(bool directed_) : directed(directed_) {}
};

inline Graph& get_graph(SEXP pointer) {
    if (TYPEOF(pointer) != EXTPTRSXP || R_ExternalPtrAddr(pointer) == NULL)
        Rcpp::stop("Graph pointer is invalid; graphs cannot be restored from serialized R objects");
    return *Rcpp::XPtr<Graph>(pointer);
}

#endif
