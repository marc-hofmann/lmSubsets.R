## Copyright 2018  Marc Hofmann and Achim Zeileis
##
## This file is part of the 'lmSubsets' R extension.
##
## 'lmSubsets' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'lmSubsets' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



###########################
##  WORKHORSE FUNCTIONS  ##
###########################


## workhorse function
##
## Arguments:
##   x         - (double[,]) model matrix
##   y         - (double[]) model response
##   weights   - (double[])
##   offset    - (double[])
##   include   - (integer[]) regressors to force in
##   exclude   - (integer[]) regressors to force out
##   nmin      - (integer)
##   nmax      - (integer)
##   tolerance - (double[])
##   nbest     - (integer) number of best subsets
##   ...       - ignored
##   pradius   - (integer) preordering radius
##
## Result: (list)
##   NOBS         - (integer)
##   nobs         - (integer)
##   nvar         - (integer)
##   weights      - (double[])
##   intercept    - (logical)
##   include      - (logical[nvar])
##   exclude      - (logical[nvar])
##   size         - (integer[])
##   tolerance    - (double[])
##   nbest        - (integer)
##   submodel     - (data.frame)
##   subset       - (data.frame)
##   .interrupted - (logical)
##   .nodes       - (integer)
##
lmSubsets_fit <- function (x, y, weights = NULL, offset = NULL,
                           include = NULL, exclude = NULL,
                           nmin = NULL, nmax = NULL, tolerance = 0,
                           nbest = 1, ..., pradius = NULL) {
    ## dims
    NOBS <- NROW(x)

    ## intercept
    has_intercept <- all(x[, 1] == 1)

    ## offset
    if (is.null(o <- offset))  o <- rep(0, NOBS)

    y <- y - o

    ## weights
    if (is.null(w <- weights))  w <- rep(1, NOBS)

    ok <- w != 0
    w <- sqrt(w[ok])

    x <- w * x[ok, , drop = FALSE]
    y <- w * y[ok,   drop = FALSE]

    ## dims
    nobs <- NROW(x)
    nvar <- NCOL(x)

    ## variables names
    x_names <- colnames(x)

    ## include
    if (is.null(include)) {
        include <- logical(nvar)
    } else {
        if (is.numeric(include)) {
            ## numeric --> numeric
            include <- match(include, seq_len(nvar))
        } else if (is.character(include)) {
            ## character --> numeric
            include <- match(include, x_names)
        } else if (is.logical(include)) {
            ## logical --> numeric
            include <- which(rep(include, length.out = nvar))
        }

        if (any(is.na(include))) {
            warning ("non-existing columns selected in 'include'")
        }

        ## canonicalize
        include <- is.element(seq_len(nvar), include)
    }

    names(include) <- x_names

    ## exclude
    if (is.null(exclude)) {
        exclude <- logical(nvar)
    } else {
        if (is.numeric(exclude)) {
            ## numeric --> numeric
            exclude <- match(exclude, seq_len(nvar))
        } else if (is.character(exclude)) {
            ## character --> numeric
            exclude <- match(exclude, x_names)
        } else if (is.logical(exclude)) {
            ## logical --> numeric
            exclude <- which(rep(exclude, length.out = nvar))
        }

        if (any(is.na(include))) {
            warning ("invalid columns selected in 'exclude'")
        }

        ## canonicalize
        exclude <- is.element(seq_len(nvar), exclude)
    }

    names(exclude) <- x_names

    ## include/exclude non-overlapping
    if (any(include & exclude)) {
        warning ("'include' and 'exclude' overlap: fixing 'exclude'")

        exclude <- exclude & !include
    }

    ## intercept
    if (has_intercept) {
        if (include[1]) {
            ## OK: included
        } else if (exclude[1]) {
            ## OK: excluded
        } else {
            ## default: include
            include[1] <- TRUE
        }
    }

    ## which
    which_incl <- as.vector(which(include))
    which_excl <- as.vector(which(exclude))

    which_free <- as.vector(which(!include & !exclude))

    if (length(which_free) < 2) {
        stop ("number of free variables must be > 1")
    }

    which_fixup <- c(which_incl, which_free, which_excl)

    ## size
    nincl <- length(which_incl)
    nexcl <- length(which_excl)

    size <- nvar - nexcl
    nmin <- max(nmin, nincl + 1)
    nmax <- min(nmax, size)

    ## preordering radius
    ## TODO: validate
    if (is.null(pradius)) {
        pradius <- floor((size - nincl) / 3)
    }

    ## tolerance
    ## TODO: validate
    tolerance <- rep(tolerance, length.out = size)
    tolerance[-seq.int(nmin, nmax)] <- +Inf

    ## call C
    ans <- local({
        algo <- list(...)$.algo
        xy <- cbind(x[, which_fixup[seq_len(size)]], y)
        mark <- nincl
        tau <- tolerance + 1.0

        .Call(C_lmSubsets, algo, xy, mark, tau, nbest, pradius)
    })

    ## value
    ans$NOBS <- NOBS
    ans$nobs <- nobs
    ans$nvar <- nvar
    ans$weights <- weights
    ans$intercept <- has_intercept
    ans$include <- include
    ans$exclude <- exclude
    ans$size <- seq.int(nmin, nmax)
    ans$tolerance <- array(tolerance, dim = nvar,
                           dimnames = list(seq_len(nvar)))
    ans$nbest <- nbest

    ans$subset <- local({
        z <- matrix(FALSE, nrow = nrow(ans$subset), ncol = nexcl)
        z[with(ans$submodel, is.na(RSS)), ] <- NA

        z <- cbind(ans$subset, z)
        z <- z[, order(which_fixup), drop = FALSE]

        colnames(z) <- x_names

        as.data.frame(z)
    })

    ## done
    ans
}



#########################
##  GENERATOR METHODS  ##
#########################


## matrix interface
##
## Arguments:
##   formula   - (double[,])
##   y         - (double[]|integer|character)
##   intercept - (logical)
##   ...       - forwarded
##
## Result: (lmSubsets) see also 'lmSubsets.default'
##
lmSubsets.matrix <- function (formula, y, intercept = TRUE, ...) {
    x_which <- seq_len(ncol(formula))
    x_names <- colnames(formula)

    ## extract expression for x
    x <- substitute(formula)

    ## extract expression for y
    if (length(y) != 1) {
        y <- substitute(y)
    } else {
        if (is.numeric(y)) {
            y <- match(y, x_which)
        } else if (is.character(y)) {
            y <- match(y, x_names)
        } else {
            stop ("'y' must be numeric or character")
        }

        if (is.na(y)) {
            stop ("non-existing column selected in 'y'")
        }

        x_which <- x_which[-y]

        y <- bquote(.(x)[, .(y)], list(x = x, y = x_names[y]))
    }

    ## build formula
    if (!is.null(x_names))  x_which <- x_names[x_which]

    f <- bquote(.(x)[, .(j)], list(x = x, j = x_which[1]))
    for (j in x_which[-1]) {
        f <- bquote(.(f) + .(x)[, .(j)], list(f = f, x = x, j = j))
    }

    if (!intercept) {
        f <- bquote(.(f) - 1, list(f = f))
    }

    f <- bquote(.(y) ~ .(f), list(y = y, f = f))
    environment(f) <- parent.frame()

    ## forward call
    cl <- match.call()
    cl[[1]] <- quote(lmSubsets)
    cl$formula <- f
    cl$y <- cl$intercept <- NULL

    eval(cl, parent.frame())
}


## standard formula interface
##
## Arguments:
##   formula   - (formula) an object that can be coerced to class 'formula'
##   data      - (data.frame)
##   subset    - (vector) subset of observations
##   weights   - (double[])
##   na.action - (function)
##   model     - (logical)
##   x         - (logical)
##   y         - (logical)
##   contrasts - (list)
##   offset    - (double[])
##   ...       - forwarded to 'lmSubsets_fit'
##
## Result: (lmSubsets) see also 'lmSubsets_fit'
##   call      - (call)
##   na.action - ()
##   offset    - (double[])
##   contrasts - (list)
##   xlevels   - (list)
##   terms     - (terms, formula)
##   model     - (data.frame) model frame
##   x         - (double[,])
##   y         - (double[])
##
lmSubsets.default <- function (formula, data, subset, weights, na.action,
                               model = TRUE, x = FALSE, y = FALSE,
                               contrasts = NULL, offset, ...) {
    ret_m <- model;  model <- NULL
    ret_x <- x;  x <- NULL
    ret_y <- y;  y <- NULL

    ## keep call (of generic)
    cl <- match.call()
    cl[[1]] <- as.name("lmSubsets")

    ## model frame
    mf <- match.call(expand.dots = FALSE)
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- quote(model.frame)
    mf$drop.unused.levels <- TRUE
    mf <- eval(mf, parent.frame())

    ## model terms
    mt <- attr(mf, "terms")

    ## model matrix and response
    x <- model.matrix(mt, mf, contrasts)
    y <- model_response(mf, "numeric")

    ## weights
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) {
        stop ("'weights' must be a numeric vector")
    }

    ## offset
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset) && (length(offset) != NROW(y))) {
        stop ("length of 'offset' must equal number of observations")
    }

    ## fit subsets
    ans <- lmSubsets_fit(x, y, weights = w, offset = offset, ...)

    ## result
    class(ans) <- "lmSubsets"

    ans$call <- cl
    ans$na.action <- attr(mf, "na.action")
    ans$offset <- offset
    ans$contrasts <- attr(x, "contrasts")
    ans$xlevels <- .getXlevels(mt, mf)
    ans$terms <- mt

    if (ret_m)  ans$model <- mf
    if (ret_x)  ans$x <- x
    if (ret_y)  ans$y <- y

    ## done
    ans
}



########################
##  STANDARD METHODS  ##
########################


## print method for 'lmSubsets' objects
##
## Arguments:
##   x   - (lmSubsets)
##   ... - forwarded
##
## Result: (lmSubsets) invisible
##
print.lmSubsets <- function (x, ...) {
    catln <- function (...)  base::cat(..., "\n", sep = "")
    print <- function (..., skip = 0, indent = 0) {
        output <- capture.output(base::print(...))
        if (skip > 0) output <- output[-seq_len(skip)]
        indent <- paste0(rep(" ", indent), collapse = "")
        cat(paste0(indent, output, "\n"), sep = "")
    }

    ## call
    catln("Call:")
    catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))

    catln()

    ## deviance
    rss <- lapply(x$size, function (n) {
        z <- with(x$submodel, RSS[SIZE == n])
        ifelse(is.na(z), "", format_default(z, ...))
    })
    rss <- do.call(cbind, rss)

    rownames(rss) <- format_ordinal(seq_len(x$nbest))
    colnames(rss) <- paste0(x$size, " (", {
        z <- x$tolerance[x$size]
        format_default(z, ...)
    }, ")")

    catln("Deviance:")
    catln("  [best, size (tolerance)] = RSS")
    print(rss, quote = FALSE, indent = 2)

    catln()

    ## subset
    subset <- lapply(names(x$subset), function (name) {
        sapply(seq_len(x$nbest), function (b) {
            z <- with(x$submodel, BEST == b)
            format_which(x$subset[z, name])
        })
    })
    subset <- do.call(rbind, subset)

    rownames(subset) <- names(x$subset)
    rownames(subset)[x$include] <- paste0("+", rownames(subset)[x$include])
    rownames(subset)[x$exclude] <- paste0("-", rownames(subset)[x$exclude])
    colnames(subset) <- format_ordinal(seq_len(x$nbest))

    catln("Subset:")
    catln("  [variable, best] = size")
    print(subset, quote = FALSE, indent = 2)

    catln()

    ## done
    invisible(x)
}


## plot method for 'lmSubsets' objects
##
## Arguments:
##   x       - (lmSubsets)
##   penalty - (character|numeric|function)
##   ...     - graphical parameters
##
## Result: (lmSubsets) invisible
##
plot.lmSubsets <- function (x, penalty = "BIC", xlim, ylim_rss, ylim_ic,
                            type_rss = "o", type_ic = "o", main, sub, xlab,
                            ylab_rss, ylab_ic, legend_rss, legend_ic,
                            ann = par("ann"), axes = TRUE, lty_rss = c(1, 3),
                            pch_rss = c(16, 21), col_rss = "black",
                            bg_rss = "white", lty_ic = c(1, 3),
                            pch_ic = c(16, 21), col_ic = "red", bg_ic = "white",
                            ...) {
    object <- x;  x <- NULL

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    par(mar = c(5, 4, 4, 4) + 0.1)

    ## RSS
    if (missing(xlim))  xlim <- range(object$size)
    if (missing(ylim_rss)) {
        ylim_rss <- with(object$submodel, range(RSS, na.rm = TRUE))
    }

    plot.window(xlim, ylim_rss, ...)

    type_rss <- rep_len(type_rss, 2)
    lty_rss <- rep_len(lty_rss, 2)
    pch_rss <- rep_len(pch_rss, 2)
    col_rss <- rep_len(col_rss, 2)
    bg_rss <- rep_len(bg_rss, 2)

    for (j in object$size) {
        x <- rep_len(j, object$nbest)
        y <- with(object$submodel, rev(RSS[SIZE == j]))

        lines(x, y, type = type_rss[2], lty = lty_rss[2],
              pch = c(rep_len(pch_rss[2], object$nbest - 1), NA),
              col = col_rss[2], bg = bg_rss[2], ...)
    }

    x <- object$size
    y <- with(object$submodel, RSS[(SIZE %in% x) & (BEST == 1)])
    lines(x, y, type = "o", lty = 1, pch = 16, col = "black", bg = "white", ...)

    ## axes
    if (axes) {
        axis(1)
        axis(2)
        box()
    }

    ## annotations
    if (ann) {
        if (missing(main))  main <- "All subsets"
        if (missing(sub)) {
            sub <- if (object$nbest > 1)  paste("nbest = ", object$nbest)
                   else NULL
        }
        if (missing(xlab))  xlab <- "size"
        if (missing(ylab_rss))  ylab_rss <- "deviance"

        title(main = main, sub = sub, xlab = xlab, ylab = ylab_rss)

        if (missing(legend_rss))  legend_rss <- "RSS"
    }

    ## criterion
    if (!is.null(penalty)) {
        ic <- ic(penalty)

        submodel <- within(object$submodel, {
            IC <- eval_ic(ic, object, size = NULL, best = NULL, na.rm = FALSE)
        })

        if (missing(ylim_ic)) {
            ylim_ic <- with(submodel, range(IC, na.rm = TRUE))
        }
        if (missing(type_ic))  type_ic <- "o"

        type_ic <- rep_len(type_ic, 2)
        lty_ic <- rep_len(lty_ic, 2)
        pch_ic <- rep_len(pch_ic, 2)
        col_ic <- rep_len(col_ic, 2)
        bg_ic <- rep_len(bg_ic, 2)

        plot.window(xlim, ylim_ic, ...)

        for (j in object$size) {
            x <- rep_len(j, object$nbest)
            y <- with(submodel, rev(IC[SIZE == j]))

            lines(x, y, type = type_ic[2], lty = lty_ic[2],
                  pch = c(rep_len(pch_ic[2], object$nbest - 1), NA),
                  col = col_ic[2], bg = bg_ic[2], ...)
        }

        x <- object$size
        y <- with(submodel, IC[(SIZE %in% x) & (BEST == 1)])
        lines(x, y, type = type_ic[1], lty = lty_ic[1], pch = pch_ic[1],
              col = col_ic[1], bg = bg_ic[1], ...)

        ## axes
        if (axes) {
            axis(4)
        }

        ## annotations
        if (ann) {
            if (missing(ylab_ic))  ylab_ic <- "criterion"
            mtext(ylab_ic, side = 4, line = 3)

            if (missing(legend_ic))  legend_ic <- format_ic(ic)
        }
    }

    ## legend
    if (ann) {
        legend <- c(legend_rss, legend_ic)
        legend("topright", legend = legend, lty = c(lty_rss[1], lty_ic[1]),
               pch = c(pch_rss[1], pch_rss[1]), col = c(col_rss[1], col_ic[1]),
               pt.bg = c(bg_rss[1], bg_ic[1]), bty = "n")
    }

    ## done
    invisible(object)
}



#########################
##  EXTRACTOR METHODS  ##
#########################


## extract variable names
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|logical[,])
##
variable.names.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                                      drop = TRUE) {
    ## full model
    x_names <- dimnames(object$which)[[1]]

    ## 'size' processing
    if (missing(size)) {
        size <- object$size
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract subsets
    ans <-  with(object$submodel, {
        z <- (SIZE %in% size) & (BEST %in% best)

        ## NA action
        if (na.rm)  z <- z & !is.na(RSS)

        cbind(SIZE = SIZE[z], BEST = BEST[z], object$subset[z, ])
    })

    rownames(ans) <- NULL

    ## drop
    if (drop) {
        ans <- structure(as.matrix(ans[-c(1, 2)]),
                         SIZE = ans$SIZE,
                         BEST = ans$BEST)

        if ((nrow(ans) == 1) && missing(drop)) {
            ans <- colnames(ans)[ans]
        }
    }

    ## done
    ans
}


## extract formula
##
## Arguments:
##   x    - (lmSubsets)
##   size - (integer)
##   best - (integer)
##   ...  - ignored
##
## Result: (formula)
##
formula.lmSubsets <- function (x, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        f <- formula(x$terms)

        return (f)
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract subset
    ans <- variable.names(x, size = size, best = best, drop = TRUE)
    ans <- as.logical(ans)

    if (x$intercept) {
        ans <- ans[-1]
    }

    ans <- x$terms[ans]
    ans <- formula(ans)

    ## done
    ans
}


## extract model frame
##
## Arguments:
##   object - (lmSubsets)
##   ...    - ignored
##
## Result: (model.frame)
##
## See also:  model.frame.lm
##
model.frame.lmSubsets <- function(formula, ...) {
    ## further arguments
    args <- list(...)
    m <- c("data", "na.action", "subset")
    m <- match(m, names(args), 0L)
    args <- args[m]

    if ((length(args) == 0) && !is.null(mf <- formula[["model"]])) {
        return (mf)
    }

    ## forward call to 'model.frame'
    cl <- formula$call
    m <- c("formula", "data", "subset", "weights", "na.action", "offset")
    m <- match(m, names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(model.frame)
    cl$drop.unused.levels <- TRUE
    cl$xlev <- formula$xlevels
    cl$formula <- terms(formula)
    cl[names(args)] <- args

    env <- environment(formula$terms)
    if (is.null(env)) {
        env <- parent.frame()
    }

    eval(cl, env)
}


## extract model matrix
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (double[,])
##
model.matrix.lmSubsets <- function (object, size, best = 1, ...) {
    ## full model
    x <- object[["x"]]
    if (is.null(x)) {
        data <- model.frame(object, xlev = object$xlevels, ...)
        x <- NextMethod("model.matrix", data = data,
                        contrasts.arg = object$contrasts)
    }

    ## 'size' processing
    if (missing(size)) {
        return (x)
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract subset
    ans <- variable.names(object, size = size, best = best, drop = TRUE)
    ans <- as.logical(ans)
    ans <- structure(x[, ans, drop = FALSE],
                     assign = attr(x, "assign")[ans])

    ## done
    ans
}


## extract model response
##
## Arguments:
##   data - (lmSubsets)
##   ...  - ignored
##
## Result: (double[])
##
model_response.lmSubsets <- function (data, ...) {
    ## extract response
    y <- data[["y"]]
    if (is.null(y)) {
        mf <- model.frame(data)
        y <- model.response(mf)
    }

    ## done
    y
}


## refit
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - forwarded
##
## Result: (lm)
##
refit.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing (size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract formula
    f <- formula(object, size = size, best = best)

    ## forward call to 'lm'
    cl <- object$call
    m <- c("formula", "data", "subset",
           "weights", "na.action", "offset")
    m <- match(m, names(cl), 0L)
    cl <- cl[c(1L, m)]
    cl[[1L]] <- quote(lm)
    cl$formula <- f

    dots <- list(...)
    cl[names(dots)] <- dots

    env <- environment(object$terms)
    if (is.null(env)) {
        env <- parent.frame()
    }

    eval(cl, env)
}


## extract deviance
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|double[])
##
deviance.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                                drop = TRUE) {
    ## 'size' processing
    if (missing(size)) {
        size <- object$size
    } else if (is.null(size)) {
        size <- seq_len(object$nvar)
    }

    ## 'best' processing
    if (is.null(best)) {
        best <- seq_len(object$nbest)
    }

    ## extract RSS
    ans <- subset(object$submodel, {
        z <- (SIZE %in% size) & (BEST %in% best)

        ## NA action
        if (na.rm)  z <- z & !is.na(RSS)

        z
    })

    rownames(ans) <- NULL

    ## drop
    if (drop) {
        ans <- with(ans, {
            structure(RSS, SIZE = SIZE, BEST = BEST)
        })

        if ((length(ans) == 1) && missing(drop)) {
            ans <- as.vector(ans)
        }
    }

    ## done
    ans
}


## extract residual standard deviation
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|double[])
##
sigma.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                             drop = TRUE) {
    ## extract sigma
    ans <- deviance(object, size = size, best = best, drop = FALSE,
                    na.rm = na.rm)
    ans <- within(ans, {
        sigma <- stats_sigma2(RSS, object$nobs - SIZE)
        sigma <- sqrt(sigma)

        rm(RSS)
    })

    ## drop
    if (drop) {
        ans <- with(ans, {
            structure(sigma, SIZE = SIZE, BEST = BEST)
        })

        if ((length(ans) == 1) && missing(drop)) {
            ans <- as.vector(ans)
        }
    }

    ## done
    ans
}


## extract log-likelihood
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|logLik)
##
logLik.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                              drop = TRUE) {
    ## extract log-likelihood
    ans <- deviance(object, size = size, best = best, drop = FALSE,
                    na.rm = na.rm)
    ans <- within(ans, {
        df <- SIZE + 1
        logLik <- stats_log_lik(object$nobs, object$weights, RSS)

        rm(RSS)
    })

    ## drop
    if (drop) {
        ans <- with(ans, {
            structure(logLik, df = df, SIZE = SIZE, BEST = BEST)
        })

        class(ans) <- "logLik"
    }

    ## decorate
    attr(ans, "nobs") <- object$nobs

    ## done
    ans
}


## compute AIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   k      - (double) penalty
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|double[])
##
AIC.lmSubsets <- function (object, size, best = 1, ..., k = 2, na.rm = TRUE,
                           drop = TRUE) {
    ## extract AIC
    ans <- logLik(object, size = size, best = best, drop = FALSE,
                  na.rm = na.rm)
    ans <- within(ans, {
        AIC <- stats_aic(logLik, k, df)

        rm(logLik)
    })

    ## drop
    if (drop) {
        ans <- with(ans, {
            structure(AIC, df = df, SIZE = SIZE, BEST = BEST)
        })

        if ((length(ans) == 1) && missing(drop)) {
            ans <- as.vector(ans)
        }
    }

    ## done
    ans
}


## compute BIC
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|double[])
##
BIC.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                           drop = TRUE) {
    ## extract BIC
    ans <- logLik(object, size = size, best = best, drop = FALSE,
                  na.rm = na.rm)
    ans <- within(ans, {
        BIC <- stats_bic(logLik, object$nobs, df)

        rm(logLik)
    })

    ## drop
    if (drop) {
        ans <- with(ans, {
            structure(BIC, df = df, SIZE = SIZE, BEST = BEST)
        })

        if ((length(ans) == 1) && missing(drop)) {
            ans <- as.vector(ans)
        }
    }

    ## done
    ans
}


## extract ceofficients
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer[])
##   best   - (integer[])
##   ...    - ignored
##   na.rm  - (logical)
##   drop   - (logical)
##
## Result: (data.frame|double[,])
##
coef.lmSubsets <- function (object, size, best = 1, ..., na.rm = TRUE,
                            drop = TRUE) {
    y <- model_response(object)

    ## extract coefficients
    ans <- variable.names(object, size = size, best = best, drop = TRUE,
                          na.rm = na.rm)
    ans[] <- do.call(rbind, lapply(seq_len(nrow(ans)), function (i) {
        size <- attr(ans, "SIZE")[i]
        best <- attr(ans, "BEST")[i]

        if (with(object$submodel, {
            is.na(RSS[(SIZE == size) & (BEST == best)])
        })) {
            return (rep_len(NA, object$nvar))
        }

        x <- model.matrix(object, size = size, best = best)

        z <- rep_len(0, object$nvar)
        z[ans[i, ]] <- stats_coef(x, y, object$offset, object$weights)

        z
    }))

    ## drop
    if (drop) {
        ## return matrix

        if ((nrow(ans) == 1) && missing(drop)) {
            ans <- ans[, ans != 0]
        }
    } else {
        ans <- cbind(
            SIZE = attr(ans, "SIZE"),
            BEST = attr(ans, "BEST"),
            as.data.frame(ans)
        )
    }

    ## done
    ans
}


## extract variance-covariance matrix
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[,])
##
vcov.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_vcov(x, y, object$offset, object$weights)

    ## names
    x_names <- variable.names(object, size = size, best = best)
    dimnames(ans) <- list(x_names, x_names)

    ## done
    ans
}


## extract fitted values
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
fitted.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_fitted(x, y, object$offset, object$weights)

    ## NA action
    ans <- napredict(object$na.action, ans)

    ## done
    ans
}


## extract residuals
##
## Arguments:
##   object - (lmSubsets)
##   size   - (integer)
##   best   - (integer)
##   ...    - ignored
##
## Result: (double[])
##
residuals.lmSubsets <- function (object, size, best = 1, ...) {
    ## 'size' processing
    if (missing(size)) {
        stop ("missing argument: 'size'")
    }

    if (length(size) > 1) {
        warning ("'size' has length > 1: only the first element will be used")

        size <- size[1]
    }

    ## 'best' processing
    if (length(best) > 1) {
        warning ("'best' has length > 1: only the first element will be used")

        best <- best[1]
    }

    ## extract submodel
    x <- model.matrix(object, size = size, best = best)
    y <- model_response(object)

    ## solve
    ans <- stats_residuals(x, y, object$offset, object$weights)

    ## NA action
    ans <- naresid(object$na.action, ans)

    ## done
    ans
}
