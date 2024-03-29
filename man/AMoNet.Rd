% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3functions.R
\name{AMoNet}
\alias{AMoNet}
\title{Initiate an AMoNet object with gene selection.}
\usage{
AMoNet(GENESman=NULL,treatmt="")
}
\arguments{
\item{GENESman}{character vector. A vector of gene(s) selected to build the *AMoNet* object.}

\item{treatmt}{character vector. A vector of gene(s) targeted by treatment(s). Used in to build the *AMoNet* object.}
}
\description{
Initiate an AMoNet object with gene selection.
}
\details{
The default \code{AMoNet()} function initiate an *AMoNet* object.
It can be invoqued either directly by inputing GENESman, or in an interactive mode by just runing \code{AMoNet()}.
Initiation consist in storing in *AMoNet* object the gene(s) and treatment(s) queries, default hyper-parameters and boundaries (and import the latter to \code{globalenv()}).
\code{AMoNet()} is also included in the \code{buid.AMoNet()} function: the user will be asked to enter gene(s)' names interactively, if not performed before.
}
\examples{
\dontrun{
net<-AMoNet() # prints an interactive question to enter gene(s)' names.
}
}
