\name{model_response}


\alias{model_response}
\alias{model_response.default}


\title{Model response}


\description{

  Extract the model response.

}


\usage{
model_response(data, ...)

\method{model_response}{default}(data, type = "any", ...)
}


\arguments{

  \item{data}{an object}

  \item{type}{\code{character}---the return type}

  \item{...}{further arguments}

}


\details{

  The default method simply forwards the call to
  \code{model.response()}.

}


\value{

  \code{double[]}---the model response

}


\seealso{

  \itemize{

    \item{\code{\link[=model.response]{model.response()}} for the
      default implementation}

  }

}
