# Data documentation

#' OMNI
#'
#' OMNI is a directed protein protein interaction (PPI) database from omnypath pooled with and immune models (in house).
#'
#' @name OMNI
#' @docType data
#' @references \url{http://omnipathdb.org/} \url{https://acsn.curie.fr/ACSN2/ACSN2.html}
#' @keywords data
NULL

#' GenesSelectImmuno
#'
#' GenesSelecImmuno are immuno related gene sets. In house collection.
#'
#' @name GenesSelectImmuno
#' @docType data
NULL

#' GenesSelectHall
#'
#' GenesSelectHall are cancer hallmarks related gene sets from MSigDB hallmarks.
#'
#' @references \url{http://software.broadinstitute.org/gsea/msigdb/}
#' @name GenesSelectHall
#' @docType data
NULL

#' names_MECA
#'
#' names_MECA are the names of gene sets from GenesSelectImmuno and GenesSelectHall
#' @name names_MECA
#' @docType data
NULL


#' Default
#'
#' Default the default hyper-parameters used for building, simulation and training of *AMoNet* object.
#' @details
#' \code{print(Default)} to visualize default hyper-parameters.
#' @name Default
#' @docType data
NULL

#' Boundaries
#'
#' Boundaries are limits of the hyper-parameter search.
#' @details
#' \code{print(Boundaries)} to visualize default hyper-parameters.
#' @name Boundaries
#' @docType data
NULL

#' CGS
#'
#' Cancer Gene Census from Cosmic database.
#' @references \url{https://cancer.sanger.ac.uk/census}
#' @details
#' Used in the \code{build()} function if \code{FilterCGS=TRUE} to reduce *AMoNet* network buiding to cancer specific genes.
#' Used in the \code{LoadCleanTCGA()} function if number of \code{Species} selected is > 800. It selects only CGS genes among user-selected \code{Species}. Motivated by problems in downloading more \code{Species}.
#' @name CGS
#' @docType data
NULL
