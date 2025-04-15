#' @title
#' Optimal Subset Identification in Undirected Weighted Network Models
#'
#' @description
#' Identify the optimal subset such that the sum of the (absolute) values of the edge weights connecting the optimal subset with its complement is maximized.
#'
#' @param adj.mat
#' The adjacency matrix. Must be a symmetric adjacency matrix.
#'
#' @param k
#' The size of the subset. Must be a non-null positive integer.
#'
#' @param method
#' The combinatorial search algorithm. Must match either `"brute.force"` or `"simulated.annealing"`. Default is `"brute.force"`.
#'
#' @param absolute
#' Whether absolute values of adj.mat, the symmetric adjacency matrix, should be computed. Must match either `FALSE` or `TRUE`. Default is `TRUE`.
#'
#' @param start.temp
#' The starting temperature in the simulated annealing search. Must be a non-null positive numeric lying within the interval (0, 1], and must be greater than stop.temp, the stopping temperature. Default is `1`.
#'
#' @param cool.fact
#' The cooling factor in the simulated annealing search. Must be a non-null positive numeric lying within the interval (0, 1). Default is `0.999`.
#'
#' @param stop.temp
#' The stopping temperature in the simulated annealing search. Must be a non-null positive numeric lying within the interval (0, 1), and must be less than start.temp, the starting temperature. Default is `0.001`.
#'
#' @param max.iter
#' The maximal number of iterations in the simulated annealing search. Must be a non-null positive integer. Default is `1000`.
#'
#' @param n.runs
#' The number of runs in the simulated annealing search. Must be a non-null positive integer. Default is `1000`.
#'
#' @param seed
#' The optional random number generator state for random number generation in the simulated annealing search. Must be a non-null positive integer. Optional. Default is `5107`.
#'
#' @param verbose
#' Whether information messages should be printed to the console. Must match either `FALSE` or `TRUE`. Default is `TRUE`.
#'
#' @return
#' A list of two named objects:
#' \describe{
#'  \item{\strong{optimal.S}}{A character vector denoting the optimal subset \eqn{S} of size \eqn{k}.}
#'  \item{\strong{optimal.function.S}}{A numeric value denoting the sum of the (absolute) values of the edge weights connecting the optimal subset \eqn{S} of size \eqn{k} with its complement \eqn{\bar{S}} of size \eqn{n - k}.}
#' }
#'
#' @details
#' In psychometrics, the network approach refers to a set of methods employed to model and analyze the relationships among psychological variables. Unlike traditional psychometric approaches, such as the structural equation approach focusing on latent variables, the network approach emphasizes the interconnections among observed variables. Due to the latter emphasis, modeling and analyzing network models to complement structural equation models when developing and evaluating psychometric instruments offers several advantages. Most notably, in undirected weighted network models, a subtype of network models, observed variables are represented by nodes, and associations between observed variables, each assigned a weight that represents the magnitude of associations, are represented by edges. In this perspective, undirected weighted network models provide estimates of the magnitude of associations (i.e., the shared variance) among items of psychometric instruments that structural equation models cannot, providing critical insight into the construct-level content validity of subsets of items of psychometric instruments. To illustrate, if an undirected weighted network model suggests that a subset of items of a psychometric instrument presents a high magnitude of associations with another subset of items, the shared variance of the said subsets of items is therefore high: the content they assess and the information they provide is highly similar. From the standpoint of the latter illustration, undirected weighted network modeling allows for the estimation of whether a subset of items assesses a “narrow” or a “broad” proportion of the construct-level content of a psychometric instrument. Hence, identifying an optimal subset of a desired number of items that assesses the “broadest” proportion of the construct-level content of a psychometric instrument consists of a combinatorial optimization problem.\cr\cr
#' Consider an undirected weighted network model \eqn{G = (V, E)}, where \eqn{V} denotes the set of nodes, and \eqn{E} denotes the set of edges. Each edge \eqn{e_{ij}} is associated with a positive or negative weight \eqn{w_{ij}}. Let \eqn{S} be a subset of nodes from \eqn{V} with a fixed size \eqn{k}, and \eqn{\bar{S}} denote the complement of \eqn{S} in \eqn{V}. The objective is to identify the optimal subset \eqn{S} of size \eqn{k} such that the sum of the (absolute) values of the edge weights connecting \eqn{S} with its complement \eqn{\bar{S}} of size \eqn{n - k} is maximized. Formally, the combinatorial optimization problem can be expressed as:\cr\cr
#' \deqn{\max_{S \subset V, |S| = k} \left( \sum_{i \in S, j \in \bar{S}} |w_{ij}| \right)}\cr
#' Solving the combinatorial optimization problem allows identifying what optimal subset of a desired number of items presents the highest magnitude of associations (i.e., the highest shared variance) within the set of items. In this light, combinatorial search algorithms (e.g., brute force search algorithm, simulated annealing search algorithm) allow to identify what optimal subset of a desired number of items should be retained in a short version of a psychometric instrument to assess the “broadest” proportion of the construct-level content of the set of items included in the original version of the said psychometric instrument.
#'
#' @examples
#' adj.mat <- (
#'   stats::runif(n = 25^2, min = -1, max = 1) |>
#'   base::matrix(nrow = 25, ncol = 25) |>
#'   (\(m) (m + base::t(m)) / 2)() |>
#'   (\(m) {base::diag(m) <- 0; m})()
#' )
#'
#' shortr::shortr(
#'   adj.mat = adj.mat,
#'   k = 5,
#'   method = base::c("brute.force"),
#'   absolute = TRUE
#' )
#'
#' @references
#' Fournier, L., Heeren, A., Baggio, S., Clark, L., Verdejo-García, A., Perales, J. C., & Billieux, J. (2025). *shortr: Optimal Subset Identification in Undirected Weighted Network Models* (Version 1.0.1) \[Computer software\]. \doi{doi:10.32614/CRAN.package.shortr}
#'
#' @encoding
#' UTF-8
#'
#' @importFrom
#' stats
#' runif
#'
#' @importFrom
#' utils
#' combn
#' setTxtProgressBar
#' txtProgressBar
#'
#' @export

shortr <- function(adj.mat, k, method = base::c("brute.force", "simulated.annealing"), absolute = TRUE, start.temp = 1, cool.fact = 0.999, stop.temp = 0.001, max.iter = 1000, n.runs = 1000, seed = 5107, verbose = TRUE) {
  
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
  
  ##### verbose, whether information messages should be printed to the console, must match either FALSE or TRUE
  
  if (!base::is.logical(verbose) || base::length(verbose) != 1) {
    
    stop("verbose, whether information messages should be printed to the console, must match either FALSE or TRUE\n")
    
  }
  
  ##### if applicable, set not available values in the symmetric adjacency matrix to null values
  
  if (base::any(base::is.na(adj.mat))) {
    
    adj.mat[base::is.na(adj.mat)] <- 0
    
  }
  
  ##### if applicable, set not available dimnames in the symmetric adjacency to a regular sequence of non-null positive integers
  
  if (base::is.null(base::dimnames(adj.mat))) {
    
    base::dimnames(adj.mat) <- base::list(base::seq_len(base::nrow(adj.mat)), base::seq_len(base::nrow(adj.mat)))
    
  }
  
  ##### sum of the (absolute) values of the edge weights connecting a subset with its complement
  
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
        
        if (verbose) {
          
          base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to zero?\n"))
          
        }
        
        return(base::invisible(base::list(optimal.S = base::character(0), optimal.function.S = 0)))
        
      } else if (k == base::nrow(adj.mat)) {
        
        if (verbose) {
          
          base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to n, the size of the set?\n"))
          
        }
        
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat), optimal.function.S = function.S(base::seq_len(base::nrow(adj.mat))))))
        
      } else {
        
        if (verbose) {
          
          base::cat(base::sprintf("\nExamining the only candidate subset...\n\n"))
          
        }
        
        best.function.S <- function.S(base::seq_len(k))
        
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[base::seq_len(k)], optimal.function.S = best.function.S)))
        
      }
      
    }
    
    best.function.S <- -Inf
    
    best.S <- NULL
    
    if (verbose) {
      
      base::cat(base::sprintf("\nExamining all %.0f subsets...\n\n", binomial.coefficient))
      
      progress.bar <- utils::txtProgressBar(min = 0, max = binomial.coefficient, style = 3)
      
    }
    
    candidate.S <- utils::combn(base::nrow(adj.mat), k)
    
    for (i in base::seq_len(base::ncol(candidate.S))) {
      
      if (verbose) {
        
        utils::setTxtProgressBar(progress.bar, i)
        
      }
      
      candidate.S.i <- candidate.S[, i]
      
      candidate.function.S.i <- function.S(candidate.S.i)
      
      if (candidate.function.S.i > best.function.S) {
        
        best.function.S <- candidate.function.S.i
        
        best.S <- candidate.S.i
        
      }
      
    }
    
    if (verbose) {
      
      utils::setTxtProgressBar(progress.bar, binomial.coefficient)
      
      base::close(progress.bar)
      
    }
    
    best.S <- base::sort(best.S)
    
    if (verbose) {
      
      base::cat("\nNumber of candidate subsets of nodes S from V with a fixed size k examined:", binomial.coefficient, "\n")
      
      base::cat("\nOptimal subset identified!\n")
      
      base::cat("\nOptimal subset:", base::paste(base::colnames(adj.mat)[best.S], collapse = " -- "), "\n")
      
      base::cat("\nOptimal value:", base::round(best.function.S, 3), "\n\n")
      
    }
    
    return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[best.S], optimal.function.S = best.function.S)))
    
  } else if (method == "simulated.annealing") {
    
    ##### simulated annealing search
    
    binomial.coefficient <- base::choose(base::nrow(adj.mat), k)
    
    if (binomial.coefficient <= 1) {
      
      if (k == 0) {
        
        if (verbose) {
          
          base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to zero?\n"))
          
        }
        
        return(base::invisible(base::list(optimal.S = base::character(0), optimal.function.S = 0)))
        
      } else if (k == base::nrow(adj.mat)) {
        
        if (verbose) {
          
          base::cat(base::sprintf("\nWhy would you set k, the size of the subset, to n, the size of the set?\n"))
          
        }
        
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat), optimal.function.S = function.S(base::seq_len(base::nrow(adj.mat))))))
        
      } else {
        
        if (verbose) {
          
          base::cat(base::sprintf("\nExamining the only candidate subset...\n\n"))
          
        }
        
        best.function.S <- function.S(base::seq_len(k))
        
        return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[base::seq_len(k)], optimal.function.S = best.function.S)))
        
      }
      
    }
    
    best.function.S <- -Inf
    
    best.S <- NULL
    
    iteration.i <- 0
    
    candidate.i <- base::new.env(hash = TRUE, parent = base::emptyenv())
    
    if (verbose) {
      
      base::cat(base::sprintf("\nExamining part of all %.0f subsets...\n\n", binomial.coefficient))
      
      progress.bar <- utils::txtProgressBar(min = 0, max = n.runs * base::min(base::ceiling(base::log(stop.temp / start.temp) / base::log(cool.fact)), max.iter), style = 3)
      
    }
    
    for (run.i in base::seq_len(n.runs)) {
      
      candidate.S.i <- base::sample.int(base::nrow(adj.mat), k)
      
      candidate.function.S.i <- function.S(candidate.S.i)
      
      best.S.run.i <- candidate.S.i
      
      best.function.S.run.i <- candidate.function.S.i
      
      temp.i <- start.temp
      
      for (i in base::seq_len(max.iter)) {
        
        if (temp.i < stop.temp || k < 2) {
          
          break
          
        }
        
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
        
        if (verbose) {
          
          utils::setTxtProgressBar(progress.bar, iteration.i)
          
        }
        
      }
      
      if (best.function.S.run.i > best.function.S) {
        
        best.function.S <- best.function.S.run.i
        
        best.S <- best.S.run.i
        
      }
      
    }
    
    if (verbose) {
      
      utils::setTxtProgressBar(progress.bar, n.runs * base::min(base::ceiling(base::log(stop.temp / start.temp) / base::log(cool.fact)), max.iter))
      
      base::close(progress.bar)
      
    }
    
    best.S <- base::sort(best.S)
    
    candidate.k <- base::length(base::ls(candidate.i))
    
    if (verbose) {
      
      base::cat("\nNumber of candidate subsets of nodes S from V with a fixed size k examined:", candidate.k, "\n")
      
      base::cat("\nOptimal subset identified!\n")
      
      base::cat("\nOptimal subset:", base::paste(base::colnames(adj.mat)[best.S], collapse = " -- "), "\n")
      
      base::cat("\nOptimal value:", base::round(best.function.S, 3), "\n\n")
      
    }
    
    return(base::invisible(base::list(optimal.S = base::colnames(adj.mat)[best.S], optimal.function.S = best.function.S)))
    
  }
  
}
