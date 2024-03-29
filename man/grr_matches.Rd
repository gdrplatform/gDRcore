% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{grr_matches}
\alias{grr_matches}
\title{Value Matching}
\usage{
grr_matches(
  x,
  y,
  all.x = TRUE,
  all.y = TRUE,
  list = FALSE,
  indexes = TRUE,
  nomatch = NA
)
}
\arguments{
\item{x}{vector.  The values to be matched.  Long vectors are not currently
supported.}

\item{y}{vector.  The values to be matched.  Long vectors are not currently
supported.}

\item{all.x}{logical; if \code{TRUE}, then each value in \code{x} will be
included even if it has no matching values in \code{y}}

\item{all.y}{logical; if \code{TRUE}, then each value in \code{y} will be
included even if it has no matching values in \code{x}}

\item{list}{logical.  If \code{TRUE}, the result will be returned as a list
of vectors, each vector being the matching values in y. If \code{FALSE},
result is returned as a data.table with repeated values for each match.}

\item{indexes}{logical.  Whether to return the indices of the matches or the
actual values.}

\item{nomatch}{the value to be returned in the case when no match is found.
If not provided and \code{indexes=TRUE}, items with no match will be
represented as \code{NA}.  If set to \code{NULL}, items with no match will
be set to an index value of \code{length+1}.  If \code{indexes=FALSE}, they will
default to \code{NA}.}
}
\value{
data.table
}
\description{
Returns a lookup table or list of the positions of ALL matches of its first
argument in its second and vice versa. Similar to \code{\link{match}}, though
that function only returns the first match.
}
\details{
This behavior can be imitated by using joins to create lookup tables, but
\code{matches} is simpler and faster: usually faster than the best joins in
other packages and thousands of times faster than the built in
\code{\link{merge}}.

\code{all.x/all.y} correspond to the four types of database joins in the
following way:

\describe{ \item{left}{\code{all.x=TRUE}, \code{all.y=FALSE}}
\item{right}{\code{all.x=FALSE}, \code{all.y=TRUE}}
\item{inner}{\code{all.x=FALSE}, \code{all.y=FALSE}}
\item{full}{\code{all.x=TRUE}, \code{all.y=TRUE}} }

Note that \code{NA} values will match other \code{NA} values.

Source of the function:
https://github.com/cran/grr/blob/master/R/grr.R
}
\examples{
mat_elem <- data.table::data.table(
  DrugName = rep(c("untreated", "drugA", "drugB", "untreated"), 2),
  DrugName_2 = rep(c("untreated", "vehicle", "drugA", "drugB"), 2),
  clid = rep(c("C1", "C2"), each = 4)
)
untreated_tag <- gDRutils::get_env_identifiers("untreated_tag")
ref_idx <- which(
  mat_elem$DrugName \%in\% untreated_tag |
   mat_elem$DrugName_2 \%in\% untreated_tag
)
ref <- mat_elem[ref_idx, ]
treated <- mat_elem[-ref_idx, ]
valid <- c("DrugName", "DrugName_2")
trt <- lapply(valid, function(x) {
  colnames <- c("clid", x) 
  treated[, colnames, with = FALSE]
})
trt <- do.call(paste, 
  do.call(rbind, lapply(trt, function(x) setNames(x, names(trt[[1]]))))
)
ref <- lapply(valid, function(x) {
  colnames <- c("clid", x) 
  ref[, colnames, with = FALSE]
})
ref <- do.call(paste, 
  do.call(rbind, lapply(ref, function(x) setNames(x, names(ref[[1]]))))
)
grr_matches(trt, ref, list = FALSE, all.y = FALSE)

}
\keyword{utils}
