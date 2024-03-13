


rewire <- function(net, p = 0, directed){
  net <- as.matrix(net)
  rownames(net) <- colnames(net) <- NULL
  if (missing(directed)){
    directed <- !all(net == t(net))
  }

  if (directed){
    edges <- c(net)
  } else {
    edges <- net[lower.tri(net,diag=FALSE)]
  }

  n_edge <- length(edges)

  for (i in which(edges!=0 & runif(n_edge) < p)){
    rewire_to <- sample(seq_len(n_edge)[-i],1)
    edge1 <- edges[i]
    edge2 <- edges[rewire_to]
    edges[i] <- edge2
    edges[rewire_to] <- edge1
  }

  if (directed){
    net[] <- edges
  } else {
    net[lower.tri(net,diag=FALSE)] <- edges
    net[upper.tri(net,diag=FALSE)] <-  t(net)[upper.tri(net,diag=FALSE)]
  }

  return(net)
}
