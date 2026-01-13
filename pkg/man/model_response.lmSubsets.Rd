\name{model_response.lmSubsets}


\alias{model_response.lmSubsets}
\alias{model_response.lmSelect}


\title{Extract the model response from a subset regression}


\description{

  Return the model response.

}


\usage{
\method{model_response}{lmSubsets}(data, ...)

\method{model_response}{lmSelect}(data, ...)
}


\arguments{

  \item{data}{\code{"lmSubsets"}, \code{"lmSelect"}---a subset
    regression}

  \item{...}{ignored}

}


\value{

  \code{double[]}---the model response

}


\seealso{

  \itemize{

    \item{\code{\link[=lmSubsets]{lmSubsets()}} for all-subsets
      regression}

    \item{\code{\link[=lmSelect]{lmSelect()}} for best-subset
      regression}

    \item{\code{\link[=model_response]{model_response()}} for the S3
      generic}

  }

}
