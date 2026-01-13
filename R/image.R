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



image.lmSubsets <- function (x, size = NULL, best = 1, which = NULL, hilite,
                             hilite_penalty, main, sub, xlab = NULL, ylab,
                             ann = par("ann"), axes = TRUE,
                             col = c("gray40", "gray90"), lab = "lab",
                             col_hilite = cbind("red", "pink"),
                             lab_hilite = "lab", pad_size = 3, pad_best = 1,
                             pad_which = 3, axis_pos = -4, axis_tck = -4,
                             axis_lab = -10, ...) {
    object <- x;  x <- NULL

    ## DATA

    ## heatmap
    if (is.null(size))  size <- object$size
    if (is.null(which))  which <- seq_len(object$nvar)

    heatmap <- variable.names(object, size = size, best = best, na.rm = TRUE,
                              drop = TRUE)
    rownames(heatmap) <- attr(heatmap, "SIZE")
    colnames(heatmap)[object$include] <-
        paste0("+", colnames(heatmap)[object$include])
    colnames(heatmap)[object$exclude] <-
        paste0("-", colnames(heatmap)[object$exclude])
    heatmap <- heatmap[, which, drop = FALSE]

    ## highlight
    if (missing(hilite)) {
        hilite <- integer(0)
    } else if (is.null(hilite)) {
        hilite <- seq_len(nrow(heatmap))
    }

    if (!missing(hilite_penalty)) {
        val <- ic(hilite_penalty)
        val <- eval_ic(val, object, size = size, best = best,
                       na.rm = TRUE, drop = TRUE)

        hilite <- order(val)[hilite]
    }

    ## PLOT

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## plot window
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

    ## color map
    col <- matrix(col, ncol = 2)
    col <- col[rep_len(seq_len(nrow(col)), nrow(heatmap)), , drop = FALSE]

    if (!is.null(col_hilite)) {
        if (is.matrix(col_hilite)) {
            ix <- rep_len(nrow(col_hilite), length(hilite))
            col_hilite <- col_hilite[ix, ]

            col[hilite, ] <- col_hilite
        } else {
            col_hilite <- rep_len(col_hilite, length(hilite))

            col[hilite, 1] <- col_hilite
        }
    }

    col <- ifelse(heatmap, col[row(heatmap), 1], col[row(heatmap), 2])

    ## padding
    pad_which <- pix2usr(x = pad_which, carry = 0)
    pad_size <- pix2usr(y = pad_size, carry = 0)
    pad_best <- pix2usr(y = pad_best, carry = 0)

    group <- match(rownames(heatmap), unique(rownames(heatmap)))

    padcnt_which <- seq_len(ncol(heatmap)) - 1
    padcnt_which <- rep(padcnt_which, each = nrow(heatmap))

    padcnt_best <- cumsum(c(0, group[-nrow(heatmap)] == group[-1]))
    padcnt_best <- rep.int(padcnt_best, ncol(heatmap))

    padcnt_size <- cumsum(c(0, group[-nrow(heatmap)] != group[-1]))
    padcnt_size <- rep.int(padcnt_size, ncol(heatmap))

    ## coords
    w <- (1 - (ncol(heatmap) - 1) * pad_which) / ncol(heatmap)
    h <- (1 - (max(padcnt_best) * pad_best) - (max(padcnt_size) * pad_size)) /
        nrow(heatmap)

    x <- (col(heatmap) - 1) * w + padcnt_which * pad_which
    y <- 1 - (row(heatmap) - 1) * h - padcnt_best * pad_best -
        padcnt_size * pad_size

    ## plot
    rect(x, y - h, x + w, y, col = col, border = NA)

    ## axes
    if (axes) {
        axis_pos <- rep_len(axis_pos, 2)
        axis_tck <- rep_len(axis_tck, 2)
        axis_lab <- rep_len(axis_lab, 2)

        ## y axis
        right <- pix2usr(x = axis_pos[2], carry = +1)
        tick <- pix2usr(x = axis_tck[2], carry = +1)
        left <- right + tick

        top <- y[which(c(TRUE, group[-nrow(heatmap)] != group[-1]))]
        bottom <- c(top[-1] + pad_size, 0)

        segments(x0 = right, y0 = top, x1 = right, y1 = bottom,
                 lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = left, y1 = top,
                 lend = 1, xpd = TRUE)

        right <- pix2usr(x = axis_lab[2], carry = +1)
        text(x = right, y = top, labels = unique(rownames(heatmap)),
             adj = c(1, 1), cex = 0.9, xpd = TRUE)

        ## x axis
        left <- x[which(row(heatmap) == 1)]
        right <- c(left[-1] - pad_which, 1)

        top <- pix2usr(y = axis_pos[1], carry = +1)
        tick <- pix2usr(y = axis_tck[1], carry = +1)
        bottom <- top + tick

        segments(x0 = left, y0 = top, x1 = right, y1 = top,
                 lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = right, y1 = bottom,
                 lend = 1, xpd = TRUE)

        lab <- rep_len(lab, 2)
        lab0 <- parse(text = lab[1])[[1]]
        lab1 <- parse(text = lab[2])[[1]]

        lab_hilite <- parse(text = lab_hilite)[[1]]

        labels <- colnames(heatmap)
        labels <- sapply(seq_along(labels), function (j) {
            ans <- labels[j]

            if (any(heatmap[, j])) {
                ans <- do.call(substitute, list(lab1, list(lab = ans)))

                if (any(heatmap[hilite, j])) {
                    ans <- do.call(substitute, list(lab_hilite,
                                                    list(lab = ans)))
                }
            } else {
                ans <- do.call(substitute, list(lab0, list(lab = ans)))
            }

            as.expression(ans)
        })

        top <- pix2usr(y = axis_lab[1], carry = +1)
        text(x = right, y = top, labels = labels, adj = c(1, 1),
             srt = 45, cex = 0.9, xpd = TRUE)
    }

    ## annotations
    if (ann) {
        if (missing(main))  main <- "All subsets"
        if (missing(sub)) {
            if (length(best) > 1) {
                sub <- paste0(format_ordinal(best), collapse = ", ")
                sub <- paste0("best = ", sub)
            } else {
                sub <- NULL
            }
        }
        if (missing(ylab))  ylab <- if (length(best) > 1) quote(size %*% best)
                                    else quote(size)

        title(main = main, sub = sub, xlab = xlab, ylab = ylab)
    }

    ## done
    invisible(object)
}


image.lmSelect <- function (x, best = NULL, which = NULL, hilite,
                            hilite_penalty, main, sub = NULL, xlab = NULL,
                            ylab, ann = par("ann"), axes = TRUE,
                            col = c("gray40", "gray90"), lab = "lab",
                            col_hilite = cbind("red", "pink"),
                            lab_hilite = "lab", pad_best = 2, pad_which = 2,
                            axis_pos = -4, axis_tck = -4, axis_lab = -10, ...) {
    object <- x;  x <- NULL

    ## DATA

    ## heatmap
    if (is.null(best))  best <- seq_len(object$nbest)
    if (is.null(which))  which <- seq_len(object$nvar)

    heatmap <- variable.names(object, best = best, drop = TRUE,
                              na.rm = TRUE)
    rownames(heatmap) <- paste0(attr(heatmap, "BEST"))
    colnames(heatmap)[object$include] <-
        paste0("+", colnames(heatmap)[object$include])
    colnames(heatmap)[object$exclude] <-
        paste0("-", colnames(heatmap)[object$exclude])
    heatmap <- heatmap[seq_len(nrow(heatmap)), which, drop = FALSE]

    ## highlight
    if (missing(hilite)) {
        hilite <- integer(0)
    } else if (is.null(hilite)) {
        hilite <- seq_along(best)
    }

    if (!missing(hilite_penalty)) {
        val <- ic(hilite_penalty)
        val <- eval_ic(val, object, best = best, na.rm = TRUE, drop = TRUE)

        hilite <- order(val)[hilite]
    }

    ## PLOT

    ## new plot
    plot.new()

    ## pars
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))

    ## margins
    par(mar = c(5, 4, 4, 4) + 0.1)

    ## plot window
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

    ## color map
    col <- matrix(col, ncol = 2)
    col <- col[rep_len(seq_len(nrow(col)), nrow(heatmap)), , drop = FALSE]

    if (!is.null(col_hilite)) {
        if (is.matrix(col_hilite)) {
            ix <- rep_len(nrow(col_hilite), length(hilite))
            col_hilite <- col_hilite[ix, ]

            col[hilite, ] <- col_hilite
        } else {
            col_hilite <- rep_len(col_hilite, length(hilite))

            col[hilite, 1] <- col_hilite
        }
    }

    col <- ifelse(heatmap, col[row(heatmap), 1], col[row(heatmap), 2])

    ## padding
    pad_which <- pix2usr(x = pad_which, carry = 0)
    pad_best <- pix2usr(y = pad_best, carry = 0)

    padcnt_which <- seq_len(ncol(heatmap)) - 1
    padcnt_which <- rep(padcnt_which, each = nrow(heatmap))

    padcnt_best <- seq_len(nrow(heatmap)) - 1
    padcnt_best <- rep.int(padcnt_best, ncol(heatmap))

    ## coords
    w <- (1 - (ncol(heatmap) - 1) * pad_which) / ncol(heatmap)
    h <- (1 - (nrow(heatmap) - 1) * pad_best) / nrow(heatmap)

    x <- (col(heatmap) - 1) * w + padcnt_which * pad_which
    y <- 1 - (row(heatmap) - 1) * h - padcnt_best * pad_best

    ## plot
    col <- col[rev(seq_len(nrow(col))), , drop = FALSE]
    rect(x, y - h, x + w, y, col = col, border = NA)

    ## axes
    if (axes) {
        axis_pos <- rep_len(axis_pos, 2)
        axis_tck <- rep_len(axis_tck, 2)
        axis_lab <- rep_len(axis_lab, 2)

        ## y axis
        right <- pix2usr(x = axis_pos[2], carry = +1)
        tick <- pix2usr(x = axis_tck[2], carry = +1)
        left <- right + tick

        top <- y[seq_len(nrow(heatmap))]
        bottom <- c(top[-1] + pad_best, 0)

        segments(x0 = right, y0 = top, x1 = right, y1 = bottom,
                 lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = left, y1 = top,
                 lend = 1, xpd = TRUE)

        right <- pix2usr(x = axis_lab[2], carry = +1)
        text(x = right, y = top, labels = rev(rownames(heatmap)),
             adj = c(1, 1), cex = 0.9, xpd = TRUE)

        ## x axis
        left <- x[which(row(heatmap) == 1)]
        right <- c(left[-1] - pad_which, 1)

        top <- pix2usr(y = axis_pos[1], carry = +1)
        tick <- pix2usr(y = axis_tck[1], carry = +1)
        bottom <- top + tick

        segments(x0 = left, y0 = top, x1 = right, y1 = top,
                 lend = 1, xpd = TRUE)
        segments(x0 = right, y0 = top, x1 = right, y1 = bottom,
                 lend = 1, xpd = TRUE)

        lab <- rep_len(lab, 2)
        lab1 <- parse(text = lab[1])[[1]]
        lab0 <- parse(text = lab[2])[[1]]

        lab_hilite <- parse(text = lab_hilite)[[1]]

        labels <- colnames(heatmap)
        labels <- sapply(seq_along(labels), function (j) {
            ans <- labels[j]

            if (any(heatmap[, j])) {
                ans <- do.call(substitute, list(lab1, list(lab = ans)))

                if (any(heatmap[hilite, j])) {
                    ans <- do.call(substitute, list(lab_hilite,
                                                    list(lab = ans)))
                }
            } else {
                ans <- do.call(substitute, list(lab0, list(lab = ans)))
            }

            as.expression(ans)
        })

        top <- pix2usr(y = axis_lab[1], carry = +1)
        text(x = right, y = top, labels = labels, adj = c(1, 1),
             srt = 45, cex = 0.9, xpd = TRUE)
    }

    ## annotations
    if (ann) {
        if (missing(main))  main <- "Best subsets"
        if (missing(ylab))  ylab <- "best"

        title(main = main, sub = sub, xlab = xlab, ylab = ylab)
    }
}
