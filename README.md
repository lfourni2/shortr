
<!-- README.md is generated from README.Rmd. Please edit that file -->

# *shortr*: Optimal Subset Identification in Undirected Weighted Network Models

<!-- badges: start -->
<!-- badges: end -->

In psychometrics, the network approach refers to a set of methods
employed to model and analyze the relationships among psychological
variables. Unlike traditional psychometric approaches, such as the
structural equation approach focusing on latent variables, the network
approach emphasizes the interconnections among observed variables. Due
to the latter emphasis, modeling and analyzing network models to
complement structural equation models when developing and evaluating
psychometric instruments offers several advantages. Most notably, in
undirected weighted network models, a subtype of network models,
observed variables are represented by nodes, and associations between
observed variables, each assigned a weight that represents the magnitude
of associations, are represented by edges. In this perspective,
undirected weighted network models provide estimates of the magnitude of
associations (i.e., the shared variance) among items of psychometric
instruments that structural equation models cannot, providing critical
insight into the construct-level content validity of subsets of items of
psychometric instruments. To illustrate, if an undirected weighted
network model suggests that a subset of items of a psychometric
instrument presents a high magnitude of associations with another subset
of items, the shared variance of the said subsets of items is therefore
high: the content they assess and the information they provide is highly
similar. From the standpoint of the latter illustration, undirected
weighted network modeling allows for the estimation of whether a subset
of items assesses a “narrow” or a “broad” proportion of the
construct-level content of a psychometric instrument. Hence, identifying
an optimal subset of a desired number of items that assesses the
“broadest” proportion of the construct-level content of a psychometric
instrument consists of a combinatorial optimization problem.

Consider an undirected weighted network model $G = (V, E)$, where $V$
denotes the set of nodes, and $E$ denotes the set of edges. Each edge
$e_{ij}$ is associated with a positive or negative weight $w_{ij}$. Let
$S$ be a subset of nodes from $V$ with a fixed size $k$, and $\bar{S}$
denote the complement of $S$ in $V$. The objective is to identify the
optimal subset $S$ of size $k$ such that the sum of the (absolute)
values of the edge weights connecting $S$ with its complement $\bar{S}$
of size $n - k$ is maximized. Formally, the combinatorial optimization
problem can be expressed as:

$$
\max_{S \subset V, |S| = k} \left( \sum_{i \in S, j \in \bar{S}} |w_{ij}| \right)
$$

Solving the combinatorial optimization problem allows identifying what
optimal subset of a desired number of items presents the highest
magnitude of associations (i.e., the highest shared variance) within the
set of items. In this light, combinatorial search algorithms (e.g.,
brute force search algorithm, simulated annealing search algorithm)
allow to identify what optimal subset of a desired number of items
should be retained in a short version of a psychometric instrument to
assess the “broadest” proportion of the construct-level content of the
set of items included in the original version of the said psychometric
instrument.

## Installation

``` r

utils::install.packages("shortr")

base::library("shortr")

```
