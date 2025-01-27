shortr <- function(adj.mat, k, method = base::c("brute.force", "simulated.annealing"), absolute = TRUE, start.temp = 1, cool.fact = 0.999, stop.temp = 0.001, max.iter = 1000, n.runs = 1000, seed = 5107) {

  ##### adj.mat, the adjacency matrix, must be a symmetric adjacency matrix

  if (!base::is.matrix(adj.mat) || !base::isSymmetric(adj.mat)) {

    stop("adj.mat, the adjacency matrix, must be a symmetric adjacency matrix\n")

  }

  ##### k, the size of the subset, must be a non-null positive integer

  if (!base::is.numeric(k) || k < 0 || k != base::round(k)) {

    stop("k, the size of the subset, must be a non-null positive integer\n")

  }

  ##### k, the size of the subset, must not be greater than n, the size of the set

  if (k > base::nrow(adj.mat)) {

    stop("k, the size of the subset, must not be greater than n, the size of the set\n")

  }

  ##### method, the combinatorial search algorithm, must match either "brute.force" or "simulated.annealing"

  if (base::length(method) != 1 || !(method %in% base::c("brute.force", "simulated.annealing"))) {

    stop('method, the combinatorial search algorithm, must match either "brute.force" or "simulated.annealing"\n')

  }

  ##### absolute, whether absolute values of adj.mat, the symmetric adjacency matrix, should be computed, must match either FALSE or TRUE

  if (!base::is.logical(absolute) || base::length(absolute) != 1) {

    stop("absolute, whether absolute values of adj.mat, the symmetric adjacency matrix, should be computed, must match either FALSE or TRUE\n")

  }

  ##### start.temp, the starting temperature, must be a non-null positive numeric lying within the interval (0, 1]

  if (!base::is.numeric(start.temp) || start.temp <= 0 || start.temp > 1) {

    stop("start.temp, the starting temperature, must be a non-null positive numeric lying within the interval (0, 1]\n")

  }

  ##### cool.fact, the cooling factor, must be a non-null positive numeric lying within the interval (0, 1)

  if (!base::is.numeric(cool.fact) || cool.fact <= 0 || cool.fact >= 1) {

    stop("cool.fact, the cooling factor, must be a non-null positive numeric lying within the interval (0, 1)\n")

  }

  ##### stop.temp, the stopping temperature, must be a non-null positive numeric lying within the interval (0, 1)

  if (!base::is.numeric(stop.temp) || stop.temp <= 0 || stop.temp >= 1) {

    stop("stop.temp, the stopping temperature, must be a non-null positive numeric lying within the interval (0, 1)\n")

  }

  ##### start.temp, the starting temperature, must be greater than stop.temp, the stopping temperature

  if (start.temp <= stop.temp) {

    stop("start.temp, the starting temperature, must be greater than stop.temp, the stopping temperature\n")

  }

  ##### max.iter, the maximal number of iterations in the simulated annealing search, must be a non-null positive integer

  if (!base::is.numeric(max.iter) || max.iter < 1 || max.iter != base::round(max.iter)) {

    stop("max.iter, the maximal number of iterations in the simulated annealing search, must be a non-null positive integer\n")

  }

  ##### n.runs, the number of runs in the simulated annealing search, must be a non-null positive integer

  if (!base::is.numeric(n.runs) || n.runs < 1 || n.runs != base::round(n.runs)) {
    stop("n.runs, the number of runs in the simulated annealing search, must be a non-null positive integer\n")
  }

  ##### seed, the optional random number generator state for random number generation in the simulated annealing search, must be a non-null positive integer

  if (!base::is.null(seed)) {

    if (!base::is.numeric(seed) || seed < 1 || seed != base::round(seed)) {

      stop("seed, the optional random number generator state for random number generation in the simulated annealing search, must be a non-null positive integer\n")

    }

    base::set.seed(seed)

  }

  ##### if applicable, set not available values in the symmetric adjacency matrix to null values

  if (base::any(base::is.na(adj.mat))) {

    adj.mat[base::is.na(adj.mat)] <- 0

  }

  ##### if applicable, set not available dimnames in the symmetric adjacency to a regular sequence of non-null positive integers

  if (base::is.null(base::dimnames(adj.mat))) {

    base::dimnames(adj.mat) <- base::list(base::seq_len(base::nrow(adj.mat)), base::seq_len(base::nrow(adj.mat)))

  }

  ##### consider an undirected weighted network model G = (V, E), where V denotes the set of nodes, and E denotes the set of edges
  ##### each edge eij is associated with a positive or negative weight wij
  ##### let S be a subset of nodes from V with a fixed size k, and S- denote the complement of S in V
  ##### function.S returns the sum of the (absolute) values of the edge weights connecting S with its complement S- of size n-k

  if (base::isFALSE(absolute)) {

    function.S <- function(S) {

      base::sum(adj.mat[S, base::setdiff(base::seq_len(base::nrow(adj.mat)), S)])

    }

  } else {

    function.S <- function(S) {

      base::sum(base::abs(adj.mat)[S, base::setdiff(base::seq_len(base::nrow(adj.mat)), S)])

    }

  }

  ##### brute force search

  if (method == "brute.force") {

    binomial.coefficient <- base::choose(base::nrow(adj.mat), k)

    if (binomial.coefficient <= 1) {

      if (k == 0) {

        base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to zero?\n"))
        return(base::invisible(base::list(optimal.S = base::character(0), optimal.function.S = 0)))

      } else if (k == base::nrow(adj.mat)) {

        base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to n, the size of the set?\n"))
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat), optimal.function.S = function.S(base::seq_len(base::nrow(adj.mat))))))

      } else {

        base::cat(base::sprintf("\nExamining the only candidate subset...\n\n"))
        best.function.S <- function.S(base::seq_len(k))
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[base::seq_len(k)], optimal.function.S = best.function.S)))

      }

    }

    best.function.S <- -Inf
    best.S <- NULL

    base::cat(base::sprintf("\nExamining all %.0f subsets...\n\n", binomial.coefficient))

    progress.bar <- utils::txtProgressBar(min = 0, max = binomial.coefficient, style = 3)

    candidate.S <- utils::combn(base::nrow(adj.mat), k)

    for (i in base::seq_len(base::ncol(candidate.S))) {

      utils::setTxtProgressBar(progress.bar, i)

      candidate.S.i <- candidate.S[, i]
      candidate.function.S.i <- function.S(candidate.S.i)

      if (candidate.function.S.i > best.function.S) {

        best.function.S <- candidate.function.S.i
        best.S <- candidate.S.i

      }

    }

    utils::setTxtProgressBar(progress.bar, binomial.coefficient)
    base::close(progress.bar)

    best.S <- base::sort(best.S)

    base::cat("\nNumber of candidate subsets of nodes S from V with a fixed size k examined:", binomial.coefficient, "\n")
    base::cat("\nOptimal subset identified!\n")
    base::cat("\nOptimal subset:", base::paste(base::colnames(adj.mat)[best.S], collapse = " -- "), "\n")
    base::cat("\nOptimal value:", base::round(best.function.S, 3), "\n\n")

    return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[best.S], optimal.function.S = best.function.S)))

  } else if (method == "simulated.annealing") {

    ##### simulated annealing search

    binomial.coefficient <- base::choose(base::nrow(adj.mat), k)

    if (binomial.coefficient <= 1) {

      if (k == 0) {

        base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to zero?\n"))
        return(base::invisible(base::list(optimal.S = base::character(0), optimal.function.S = 0)))

      } else if (k == base::nrow(adj.mat)) {

        base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to n, the size of the set?\n"))
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat), optimal.function.S = function.S(base::seq_len(base::nrow(adj.mat))))))

      } else {

        base::cat(base::sprintf("\nExamining the only candidate subset...\n\n"))
        best.function.S <- function.S(base::seq_len(k))
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[base::seq_len(k)], optimal.function.S = best.function.S)))

      }

    }

    best.function.S <- -Inf
    best.S <- NULL
    iteration.i <- 0

    candidate.i <- base::new.env(hash = TRUE, parent = base::emptyenv())

    base::cat(base::sprintf("\nExamining part of all %.0f subsets...\n\n", binomial.coefficient))

    progress.bar <- utils::txtProgressBar(min = 0, max = n.runs * base::min(base::ceiling(base::log(stop.temp / start.temp) / base::log(cool.fact)), max.iter), style = 3)

    for (run.i in base::seq_len(n.runs)) {

      candidate.S.i <- base::sample.int(base::nrow(adj.mat), k)
      candidate.function.S.i <- function.S(candidate.S.i)

      best.S.run.i <- candidate.S.i
      best.function.S.run.i <- candidate.function.S.i

      temp.i <- start.temp

      for (i in base::seq_len(max.iter)) {

        if (temp.i < stop.temp || k < 2)

          break

        candidate.S.j <- candidate.S.i
        candidate.S.j[base::sample.int(k, 1)] <- base::sample(base::setdiff(base::seq_len(base::nrow(adj.mat)), candidate.S.i), 1)
        candidate.function.S.j <- function.S(candidate.S.j)

        candidate.j <- base::paste(base::sort(candidate.S.j), collapse = "-")

        if (!base::exists(candidate.j, envir = candidate.i)) {

          base::assign(candidate.j, TRUE, envir = candidate.i)

        }

        if (candidate.function.S.j > candidate.function.S.i || stats::runif(1) < base::exp((candidate.function.S.j - candidate.function.S.i) / temp.i)) {

          candidate.S.i <- candidate.S.j
          candidate.function.S.i <- candidate.function.S.j

          if (candidate.function.S.j > best.function.S.run.i) {

            best.function.S.run.i <- candidate.function.S.j
            best.S.run.i <- candidate.S.j

          }

        }

        temp.i <- temp.i * cool.fact
        iteration.i <- iteration.i + 1
        utils::setTxtProgressBar(progress.bar, iteration.i)

      }

      if (best.function.S.run.i > best.function.S) {

        best.function.S <- best.function.S.run.i
        best.S <- best.S.run.i

      }

    }

    utils::setTxtProgressBar(progress.bar, n.runs * base::min(base::ceiling(base::log(stop.temp / start.temp) / base::log(cool.fact)), max.iter))
    base::close(progress.bar)

    best.S <- base::sort(best.S)

    candidate.k <- base::length(base::ls(candidate.i))

    base::cat("\nNumber of candidate subsets of nodes S from V with a fixed size k examined:", candidate.k, "\n")
    base::cat("\nOptimal subset identified!\n")
    base::cat("\nOptimal subset:", base::paste(base::colnames(adj.mat)[best.S], collapse = " -- "), "\n")
    base::cat("\nOptimal value:", base::round(best.function.S, 3), "\n\n")

    return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[best.S], optimal.function.S = best.function.S)))

  }

}
