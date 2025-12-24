
# capture fit measures without printing

quiet_fit <- function(model, ...) {
  utils::capture.output( res <- fit(model, ...) )
  res
}

# calculate Frobenius norm based similarity measure for any pair
# This serves as an effect size measure according to Ulitzsch et al. (2023)
sF <- function(mat) {

  # focus on beta and cont
  net <- c("beta", "contemporaneous")
  # numbero of groups
  n_g <- length(mat$beta)
  # group names
  g <- names(mat$beta)
  # all possible comparison pairs
  comb <- utils::combn(n_g, 2, simplify = FALSE)

  # calculate norm
  norms <- lapply(net, function(n) {
    lapply(comb, function(k) {
      i <- k[1]; j <- k[2]; p = nrow(mat[[n]][[i]])
      data.frame(
        i = g[i],
        j = g[j],
        # calculate the similarity measure in Ulitzsch et al. (2023)
        sf = (1 / (1 + norm(mat[[n]][[i]] - mat[[n]][[j]], type = "F") / sqrt(p / 2))) %>% round(3)
      )
    }) %>% bind_rows
  }) %>% setNames(net)

  return(norms)

}
