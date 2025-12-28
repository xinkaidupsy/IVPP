#' Centrality (degree/EI) and bridge centrality for an network matrix
#'
#' Computes node centrality measures from a weighted network matrix.
#' Degree and expected-influence measures are obtained via
#' \code{\link[qgraph]{centrality}} and bridge centrality is obtained via
#' \code{\link[networktools]{bridge}} when communities are supplied.
#'
#' @param graph A square numeric network matrix (`p x p`). Can be weighted and/or signed.
#'   If directed, provide a non-symmetric matrix. Row and column names should be identical
#'   and in the same order (node names).
#' @param communities A character vector indicating community membership of the nodes
#' (e.g., c("Comm1", "Comm1", "Comm2", "Comm2)).
#' @param alpha The tuning parameter. Defaults to 1.
#' @param posfun 	A function that converts positive and negative values to only positive.
#' Defaults to the absolute value.
#' @param weighted Logical, set to FALSE to set all edge weights to 1 or -1
#' @param signed Logical, set to FALSE to make all edge weights absolute.
#' @param directed Logical. Whether the input network is directed.
#' Automatically detected if set to "NULL" (the default).
#' Symmetric network matrices will be undirected, asymmetric matrices will be directed
#' @param bridge_normalize logical. Bridge centralities are divided by their highest possible value
#' (assuming max edge strength=1) in order to normalize by different community sizes
#' @param useCommunities Character string passed to \code{networktools::bridge(useCommunities = ...)}.
#'
#' @return A list with components:
#' \describe{
#'   \item{degree}{Named list of numeric vectors (length `p`):
#'     `OutDegree`, `InDegree`, `OutExpectedInfluence`, `InExpectedInfluence`.}
#'   \item{bridge}{If `communities` is provided, the object returned by
#'     \code{networktools::bridge()}; otherwise `NULL`.}
#'   \item{meta}{List with metadata: `directed`, `alpha`, `weighted`, `signed`,
#'     and the expression name of `posfun`.}
#' }
#'
#' @seealso \code{\link[qgraph]{centrality}}, \code{\link[networktools]{bridge}}
#'
#' @examples
#' # Undirected signed weighted example:
#' W <- matrix(c(
#'   0,  0.4, -0.2,
#'   0.4, 0,   0.1,
#'  -0.2, 0.1, 0
#' ), nrow = 3, byrow = TRUE)
#' colnames(W) <- rownames(W) <- c("A","B","C")
#'
#' cent <- centrality(W)
#' cent$degree$OutDegree
#' cent$degree$OutExpectedInfluence
#'
#' # Bridge centrality:
#' comm <- c(A = "X", B = "X", C = "Y")
#' cent_b <- centrality(W, communities = comm)
#' @importFrom qgraph qgraph centrality
#' @importFrom networktools bridge
#' @export centrality

centrality <- function(
    graph,
    communities = NULL,
    alpha = 1,
    posfun = abs,
    weighted = TRUE,
    signed = TRUE,
    directed = NULL,
    bridge_normalize = FALSE,
    useCommunities = "all"
) {
  # ---- minimal matrix sanitation ----
  if (!is.matrix(graph)) stop("`graph` must be a network matrix.")
  if (nrow(graph) != ncol(graph)) stop("`graph` must be square (p x p).")

  p <- ncol(graph)

  # Node names (needed for consistent outputs)
  nodes <- rownames(graph)
  if (is.null(nodes)) nodes <- colnames(graph)
  if (is.null(nodes)) nodes <- paste0("V", seq_len(p))
  rownames(graph) <- colnames(graph) <- nodes

  g <- qgraph::qgraph(graph, DoNotPlot = TRUE, directed = directed, labels = nodes)

  # Directedness (used for meta + passed to bridge)
  if (is.null(directed)) {
    directed <- !isTRUE(all.equal(graph, t(graph)))
  }

  # ---- Degree + Expected Influence ----

  cent_q <- qgraph::centrality(
    g,
    alpha = alpha,
    posfun = posfun,
    weighted = weighted,
    signed = signed
  )

  keep <- c("OutDegree","InDegree","OutExpectedInfluence","InExpectedInfluence")
  degree_list <- setNames(lapply(keep, function(nm) {
    x <- cent_q[[nm]]
    if (!is.null(names(x))) x[nodes] else setNames(as.numeric(x), nodes)
  }), keep)


  # ---- Bridge centrality (networktools) ----
  bridge_list <- NULL
  if (!is.null(communities)) {
    if (!requireNamespace("networktools", quietly = TRUE)) {
      stop("`communities` was provided, but the networktools package is not installed.")
    }

    # If communities is a named vector, align it to node order
    if (is.atomic(communities) && !is.null(names(communities))) {
      communities <- as.character(communities[nodes])
    }

    bridge_list <- networktools::bridge(
      network = g,
      communities = communities,
      useCommunities = useCommunities,
      directed = directed,
      nodes = nodes,
      normalize = bridge_normalize
    )
    # Note: bridge strength uses abs(weights); bridge betweenness/closeness only for positive weights;
    # negative edges are dropped for those metrics. :contentReference[oaicite:2]{index=2}
  }

  list(
    degree = degree_list,
    bridge = bridge_list,
    meta = list(
      directed = directed,
      alpha = alpha,
      weighted = weighted,
      signed = signed,
      posfun = deparse(substitute(posfun))
    )
  )
}


#' Radar plot for selected centrality measures
#'
#' Creates a radar chart (via \code{\link[fmsb]{radarchart}}) from the output of
#' \code{\link{centrality}}.
#'
#' @param cent A list returned by \code{\link{centrality}}.
#' @param measures Character vector of measures to plot. Allowed values:
#'   \itemize{
#'     \item `outdegree`, `indegree`, `outEI`, `inEI`
#'     \item `bridgestrength`, `bridgeoutdegree`, `bridgeindegree`,
#'           `bridgebetweenness`, `bridgecloseness`, `bridgeEI_1`, `bridgeEI_2`
#'   }
#' @param title Optional character title. If provided and `...` does not already set `title`,
#'   it is passed to \code{fmsb::radarchart(title = ...)}.
#' @param axis_range Numeric length-2 vector `c(min, max)` for the radial axis. If `NULL`,
#'   computed from the selected measures.
#' @param add_legend Logical. If `TRUE`, adds a base-R legend with one entry per selected measure.
#' @param legend_pos Character position, e.g. `"bottomright"`, `"topright"`, etc.
#' @param ... Additional arguments forwarded to \code{\link[fmsb]{radarchart}}.
#'
#' @return Invisibly returns the data frame used for plotting (in \pkg{fmsb} format):
#' the first row is `max`, the second row is `min`, and subsequent rows correspond
#' to the requested centrality measures.
#'
#' @seealso \code{\link{centrality}}, \code{\link[fmsb]{radarchart}}
#'
#' @examples
# Minimal example (requires fmsb installed to actually draw the plot)
#' W <- matrix(c(
#'   0,  0.4, -0.2, 0,
#'   0.4, 0,   0.1, 0.3,
#'   -0.2, 0.1, 0,   0.2,
#'   0,   0.3, 0.2, 0
#' ), nrow = 4, byrow = TRUE)
#' colnames(W) <- rownames(W) <- c("A", "B", "C", "D")
#'
#' cent <- centrality(W)
#'
#' # Plot degree + EI
#' # pdf('radar_plot.pdf')
#' centrality_radar(cent, measures = c("outdegree", "outEI"),
#'                  axis_range = c(0, 1), plwd = 2,
#'                  seg = 5, caxislabel = seq(0,1,0.2))
#' # dev.off()
#'
#' # If bridge metrics were computed:
#'
#' comm <- c(A="X", B="X", C="Y", D="Y")
#' cent_b <- centrality(W, communities = comm)
#' centrality_radar(cent_b, measures = c("bridgestrength", "bridgebetweenness"))
#'
#' @importFrom fmsb radarchart
#' @importFrom graphics legend
#' @export centrality_radar

centrality_radar <- function(
    cent,
    measures = c(
      "outdegree",
      "indegree",
      "outEI",
      "inEI",
      "bridgestrength",
      "bridgeoutdegree",
      "bridgeindegree",
      "bridgebetweenness",
      "bridgecloseness",
      "bridgeEI_1",
      "bridgeEI_2"
    ),
    title = NULL,
    axis_range = NULL,          # c(min, max). If NULL, computed from selected measures
    add_legend = TRUE,
    legend_pos = "bottomright",
    ...
) {

  if(missing(cent)) {
    stop('Must input centrality list computed with centrality()')
  }

  # measures
  measures <- match.arg(measures, several.ok = TRUE)

  # --- map measures with keys in the centrality() output ---
  dict <- list(
    outdegree = c("OutDegree"),
    indegree = c("InDegree"),
    outEI = c("OutExpectedInfluence"),
    inEI = c("InExpectedInfluence"),

    bridgestrength = c("Bridge Strength"),
    bridgebetweenness = c("Bridge Betweenness"),
    bridgecloseness = c("Bridge Closeness"),
    bridgeEI_1 = c("Bridge Expected Influence (1-step)"),
    bridgeEI_2 = c("Bridge Expected Influence (2-step)"),

    bridgeoutdegree = c("Bridge Outdegree"),
    bridgeindegree = c("Bridge Indegree")
  )

  # --- helper: fetch one measure vector (named by nodes) from cent$degree or cent$bridge ---
  get_measure <- function(k) {
    for (key in dict[[k]]) {
      if (!is.null(cent$degree[[key]])) return(cent$degree[[key]])
      if (!is.null(cent$bridge[[key]])) return(cent$bridge[[key]])
    }
    NULL
  }

  # --- collect requested measures ---
  c_ls <- lapply(measures, get_measure)

  # If any requested measures are missing, stop with a clear message
  missing_meas <- measures[vapply(c_ls, is.null, FALSE)]
  if (length(missing_meas) > 0) {
    stop(
      "These centrality measures were requested but not found: ",
      paste(missing_meas, collapse = ", "),
      "\nNote: 'bridgeoutdegree' / 'bridgeindegree' only work if your network is directed"
    )
  }

  # --- build a matrix: rows = measures (series), columns = nodes (axes) ---
  # nodes
  nodes <- names(c_ls[[1]])

  # matrix of centrality
  c_mat <- do.call(rbind, lapply(c_ls, function(v) as.numeric(v[nodes])))
  rownames(c_mat) <- measures
  colnames(c_mat) <- nodes

  if (ncol(c_mat) < 3) stop("Radar plot needs at least 3 nodes (axes).")

  # --- choose radial axis range ---
  if (is.null(axis_range)) {
    axis_range <- range(c_mat, finite = TRUE, na.rm = TRUE)
  } else {
    if (length(axis_range) != 2) stop("`axis_range` must be length-2: c(min, max).")
    axis_range <- sort(axis_range)
  }

  # --- fmsb formated df ---
  df_plot <- as.data.frame(rbind(
    max = rep(axis_range[2], ncol(c_mat)),
    min = rep(axis_range[1], ncol(c_mat)),
    c_mat
  ))
  colnames(df_plot) <- nodes

  # --- pass ... to radarchart(), with a couple of defaults ---
  args <- list(df = df_plot)
  dots <- list(...)
  args[names(dots)] <- dots
  if (is.null(args$axistype)) args$axistype <- 1
  if (!is.null(title) && is.null(args$title)) args$title <- title

  do.call(fmsb::radarchart, args)

  # --- legend: one line per selected measure ---
  if (add_legend) {
    # Use user-supplied line colors/types/widths if provided; otherwise make simple defaults.
    pcol <- if (!is.null(args$pcol)) rep(args$pcol, length.out = length(measures)) else seq_along(measures)
    plty <- if (!is.null(args$plty)) rep(args$plty, length.out = length(measures)) else rep(1, length(measures))
    plwd <- if (!is.null(args$plwd)) rep(args$plwd, length.out = length(measures)) else rep(1, length(measures))

    legend(legend_pos, legend = measures, col = pcol, lty = plty, lwd = plwd, bty = "n")
  }

  invisible(df_plot)
}

#' Bridge centrality by target community
#'
#' Computes bridge strength by community pairs
#' @param graph p x p adjacency matrix (can be directed or undirected)
#' @param community character vector of length p in the SAME order as
#' rows/cols of graph e.g., rep(c("dep","anx","mech"), each = 3)
#' @export bridge_by_community
#'
bridge_by_community <- function(graph, community) {


  # --- sanity checks ---
  #
  if (!is.matrix(graph)) stop("graph must be a matrix.")
  if (nrow(graph) != ncol(graph)) stop("graph must be square (p x p).")
  p <- nrow(graph)

  # check community input
  if (length(community) != p) stop("community must have length nrow(graph).")
  community <- as.character(community)

  # set node names
  nodes <- rownames(graph)
  if (is.null(nodes)) nodes <- colnames(graph)
  if (is.null(nodes)) nodes <- paste0("V", seq_len(p))
  rownames(graph) <- colnames(graph) <- nodes

  # NA to 0 & ignore self-loops
  graph[is.na(graph)] <- 0
  diag(graph) <- 0

  absW <- abs(graph)
  comms <- unique(community)

  # --- calculate strength and EI by target community ---
  out_strength <- matrix(NA, p, length(comms), dimnames = list(nodes, comms))
  in_strength <- matrix(NA, p, length(comms), dimnames = list(nodes, comms))
  out_ei <- matrix(NA, p, length(comms), dimnames = list(nodes, comms))
  in_ei <- matrix(NA, p, length(comms), dimnames = list(nodes, comms))

  for (C in comms) {
    idx <- which(community == C)

    # "out" uses row i -> columns in target community C
    out_strength[, C] <- rowSums(absW[, idx, drop = FALSE])
    out_ei[, C] <- rowSums(graph[, idx, drop = FALSE])

    # "in" uses rows in target community C -> column i
    in_strength[, C] <- colSums(absW[idx, , drop = FALSE])
    in_ei[, C] <- colSums(graph[idx, , drop = FALSE])

    # within-community entries are not "bridge" edges, set to NA
    out_strength[idx, C] <- NA
    out_ei[idx, C] <- NA
    in_strength[idx, C] <- NA
    in_ei[idx, C] <- NA
  }

  colnames(out_strength) <- paste0("bridge_strength_out_to_", comms)
  colnames(in_strength) <- paste0("bridge_strength_in_from_", comms)
  colnames(out_ei) <- paste0("bridge_EI_out_to_", comms)
  colnames(in_ei) <- paste0("bridge_EI_in_from_", comms)

  data.frame(
    node = nodes,
    community = community,
    out_strength,
    total_bridge_strength_out = rowSums(out_strength, na.rm = TRUE),
    in_strength,
    total_bridge_strength_in = rowSums(in_strength, na.rm = TRUE),
    out_ei,
    total_bridge_EI_out = rowSums(out_ei, na.rm = TRUE),
    in_ei,
    total_bridge_EI_in = rowSums(in_ei, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE
  )
}
