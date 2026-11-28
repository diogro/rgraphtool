# Common accessors for every block-state family, so callers need no class
# switches of their own.

#' @noRd
#' @export
description_length <- function(state, ...) UseMethod("description_length")

#' @export
description_length.default <- function(state, ...) {
  stop("no description length for an object of class ",
       paste(class(state), collapse = "/"), call. = FALSE)
}
#' @export
description_length.rgraphtool_block_state <- function(state, ...) block_description_length(state, ...)
#' @export
description_length.rgraphtool_nested_block_state <- function(state, ...) nested_description_length(state, ...)
#' @export
description_length.rgraphtool_overlap_block_state <- function(state, ...) overlap_description_length(state, ...)
#' @export
description_length.rgraphtool_nested_overlap_block_state <- function(state, ...) {
  nested_overlap_description_length(state, ...)
}
#' @export
description_length.rgraphtool_weighted_block_state <- function(state, ...) weighted_description_length(state, ...)
#' @export
description_length.rgraphtool_nested_weighted_block_state <- function(state, ...) {
  nested_weighted_description_length(state, ...)
}
#' @export
description_length.rgraphtool_layered_block_state <- function(state, ...) layered_description_length(state, ...)
#' @export
description_length.rgraphtool_nested_layered_block_state <- function(state, ...) {
  nested_layered_description_length(state, ...)
}
#' @export
description_length.rgraphtool_pp_block_state <- function(state, ...) pp_description_length(state, ...)
#' @export
description_length.rgraphtool_ranked_block_state <- function(state, ...) ranked_description_length(state, ...)
#' @export
description_length.rgraphtool_nested_ranked_block_state <- function(state, ...) {
  nested_ranked_description_length(state, ...)
}
#' @export
description_length.rgraphtool_bipartite_block_state <- function(state, ...) bipartite_description_length(state)
#' @export
description_length.rgraphtool_nested_bipartite_block_state <- function(state, ...) {
  nested_bipartite_description_length(state)
}
#' @export
description_length.rgraphtool_dynamics_block_state <- function(state, ...) dynamics_description_length(state, ...)
#' @export
description_length.rgraphtool_latent_closure_block_state <- function(state, ...) {
  closure_description_length(state, ...)
}
#' @export
description_length.rgraphtool_uncertain_block_state <- function(state, ...) {
  uncertain_description_length(state, ...)
}

.nested_state_classes <- c("rgraphtool_nested_block_state", "rgraphtool_nested_weighted_block_state",
                           "rgraphtool_nested_layered_block_state",
                           "rgraphtool_nested_ranked_block_state")
.flat_state_classes <- c("rgraphtool_block_state", "rgraphtool_weighted_block_state",
                         "rgraphtool_layered_block_state", "rgraphtool_pp_block_state",
                         "rgraphtool_ranked_block_state", "rgraphtool_bipartite_block_state",
                         "rgraphtool_modularity_state")

#' @noRd
#' @export
state_blocks <- function(state, level = 0L) {
  if (inherits(state, .nested_state_classes)) return(nested_project(state, level))
  if (inherits(state, "rgraphtool_dynamics_block_state")) {
    # the SBM prior over the reconstructed graph: a hierarchy in state$bs
    bs <- state$bs
    depth <- length(bs)
    if (!is.numeric(level) || length(level) != 1L || is.na(level) || level < 0 ||
        level >= depth || level != floor(level)) {
      stop(sprintf("level must be an integer between 0 and %d", depth - 1L), call. = FALSE)
    }
    b <- as.integer(bs[[1L]])
    for (l in seq_len(level)) b <- as.integer(bs[[l + 1L]][b])
    return(b)
  }
  if (inherits(state, .flat_state_classes)) {
    if (!identical(as.integer(level), 0L)) {
      stop("a flat state has only level 0", call. = FALSE)
    }
    return(block_membership(state))
  }
  stop("state_blocks() does not support objects of class ",
       paste(class(state), collapse = "/"), call. = FALSE)
}

#' @noRd
#' @export
is_nested_state <- function(state) inherits(state, .nested_state_classes)

#' @noRd
#' @export
is_block_state <- function(x) {
  inherits(x, c(.nested_state_classes, .flat_state_classes,
                "rgraphtool_overlap_block_state", "rgraphtool_nested_overlap_block_state",
                "rgraphtool_nested_bipartite_block_state", "rgraphtool_dynamics_block_state",
                "rgraphtool_latent_closure_block_state", "rgraphtool_uncertain_block_state"))
}
