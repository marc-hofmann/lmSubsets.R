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



format_default <- function (x, na.encode = TRUE, ...) {
    ans <- format(x, na.encode = TRUE, ...)

    if (!na.encode) {
        ans <- ifelse(is.na(x), NA, ans)
    }

    ans
}


format_pval <- function (x, na.encode = TRUE, ...) {
    ans <- format.pval(x, ...)

    if (!na.encode) {
        ans <- ifelse(is.na(x), NA, ans)
    }

    ans
}


format_ordinal <- function (x) {
    ind <- c("th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "th", "th", "th", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th",
             "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th")
    paste0(x, ind[x %% 100 + 1])
}


format_which <- function (mask, seq_ind = "-", seq_sep = ",") {
    w <- which(mask)
    if (length(w) < 1) {
        return ("")
    }

    d <- c(diff(w), -1)
    z <- character(0)

    state <- 1
    pos <- 1
    while (state > 0) {
        switch (state, {
            ## state = 1
            if (d[pos] < 0) {
                z <- c(z, w[pos])
                state <- -1
            } else if (d[pos] == 1) {
                z <- c(z, w[pos], seq_ind)
                state <- 2
            } else {
                z <- c(z, w[pos], seq_sep)
                state <- 1
            }
        }, {
            ## state = 2
            if (d[pos] < 0) {
                z <- c(z, w[pos])
                state <- -1
            } else if (d[pos] == 1) {
                state <- 2
            } else {
                z <- c(z, w[pos], seq_sep)
                state <- 1
            }
        })

        pos <- pos + 1
    }

    paste(z, collapse = "")
}
