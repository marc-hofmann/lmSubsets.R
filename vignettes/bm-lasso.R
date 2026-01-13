library(lmSubsets)
library(glmnet)


## run lasso
lasso <- function (x, y) {
    ans <- glmnet(x = x, y = y, family = "gaussian", alpha = 1)
    ans$x <- x
    ans$y <- y

    class(ans) <- c("lasso", class(ans))
    ans
}


## compute sparse OLS estimator and return RSS
deviance.lasso <- function (object, index = NULL, ...) {
    if (is.null(index))  index <- seq_along(object$lambda)

    if (length(index) > 1) {
        ans <- sapply(index, function (i) deviance(object, index = i, ...))
        attr(ans, "ix") <- index

        return (ans)
    }


    x <- object$x
    y <- object$y

    beta <- object$beta[, index]
    which <- beta != 0
    size <- object$df[index]

    ans <- qr.R(qr(cbind(1, x[, which], y)))[size + 2, size + 2]^2
    attr(ans, "ix") <- index

    ans
}


## compute log-likelihood of sparse OLS estimator
logLik.lasso <- function (object, index = NULL, ...) {
    m <- object$nobs
    rss <- deviance(object, index = index)

    ans <- 0.5 * (-m * (log(2 * pi) + 1 - log(m) + log(rss)))
    attr(ans, "df") <- object$df[attr(rss, "ix")] + 2

    ans
}


## compute BIC of sparse OLS estimator
BIC.lasso <- function (object, index = NULL, ...) {
    ll <- logLik(object, index = index)
    m <- object$nobs
    k <- log(m)

    -2 * ll + k * attr(ll, "df")
}


## simulate dataset and run hook
simu <- function (seed, nobs, nvar, ntrue, sd, hook) {
    intercept <- TRUE

    set.seed(seed)

    x <- rnorm(nobs * nvar)
    dim(x) <- c(nobs, nvar)

    coefs <- rep_len(0, length.out = nvar)
    coefs[sample.int(nvar, size = ntrue)] <- 1

    error <- rnorm(nobs, sd = sd)

    y <- cbind(intercept, x) %*% c(1, coefs) + error

    hook(x, as.numeric(y))
}


## run lasso on simulated dataset
run_lasso <- function (seed, nobs, nvar, ntrue, sd, nbest = 1) {
    hook <- function (x, y) {
        time <- system.time(val <- lasso(x, y))

        bic <- sort(unique(BIC(val)))[seq.int(nbest)]
        bic <- bic[!is.na(bic)]

        list(tm = summary(time)["elapsed"], bic = bic)
    }

    simu(seed, nobs, nvar, ntrue, sd, hook)
}


## run 'lmSelect' on simulated dataset
run_lmSelect <- function (seed, nobs, nvar, ntrue, sd, nbest = 1) {
    hook <- function (x, y) {
        time <- system.time(value <- lmSelect(x, y, nbest = nbest))

        bic <- BIC(value, best = seq.int(nbest))

        list(tm = summary(time)["elapsed"], bic = bic)
    }

    simu(seed, nobs, nvar, ntrue, sd, hook)
}


## run simulation
run_benchmark <- function (seed = 1) {
    nbest <- 10
    tol <- 0.001

    NOBS <- 1000
    SD <- c(0.05, 0.1, 0.5, 1.0)
    NVAR <- c(60, 100, 140)
    NREP <- 5

    set.seed(seed)

    runs <- expand.grid(REP = 1:NREP, NVAR = NVAR, NOBS = NOBS, SD = SD)
    nruns <- nrow(runs)

    runs <- cbind(RUN_ID = 1:nruns,
                  SEED = sample.int(.Machine$integer.max, size = nruns),
                  runs, NTRUE = runs$NVAR %/% 2, TM_LMSELECT = NA,
                  BIC_LMSELECT = NA, TM_LASSO = NA, BIC_LASSO = NA,
                  HITS_LASSO = NA)

    message ("NRUNS: ", nruns)

    for (id in runs$RUN_ID) {
        ix <- with(runs, which(RUN_ID == id))

        seed <- runs[ix, "SEED"]
        nobs <- runs[ix, "NOBS"]
        nvar <- runs[ix, "NVAR"]
        ntrue <- runs[ix, "NTRUE"]
        sd <- runs[ix, "SD"]

        message ("RUN: ", id, " (SD: ", sd, ", NVAR: ", nvar, ")")

        r_lmSelect <- run_lmSelect(seed, nobs, nvar, ntrue, sd, nbest)

        message ("  lmSelect: ", format(r_lmSelect$tm, digits = 1, nsmall = 2),
                 "  ", format(r_lmSelect$bic[1], digits = 1, nsmall = 4))

        runs[ix, "BIC_LMSELECT"] <- r_lmSelect$bic[1]
        runs[ix, "TM_LMSELECT"] <- r_lmSelect$tm

        r_lasso <- run_lasso(seed, nobs, nvar, ntrue, sd, nbest)

        hits <- sum(r_lasso$bic - tol < r_lmSelect$bic[seq_along(r_lasso$bic)] &
                    r_lasso$bic + tol > r_lmSelect$bic[seq_along(r_lasso$bic)])

        message ("  lasso: ", format(r_lasso$tm, digits = 1, nsmall = 2),
                 "  ", format(r_lasso$bic[1], digits = 1, nsmall = 4),
                 "  ", hits)

        runs[ix, "BIC_LASSO"] <- r_lasso$bic[1]
        runs[ix, "TM_LASSO"] <- r_lasso$tm

        runs[ix, "HITS_LASSO"] <- hits
    }

    runs
}


report_benchmark <- function () {
    T_SUM <- bm
    T_SUM <- split(T_SUM, T_SUM[, c("SD", "NVAR")])
    T_SUM <- lapply(T_SUM, function (grp) with(grp, {
        cbind(SD = SD[1], NVAR = NVAR[1],
              LMSELECT = mean(TM_LMSELECT),
              LASSO = mean(TM_LASSO),
              HITS = mean(HITS_LASSO))
    }))
    T_SUM <- do.call(rbind, T_SUM)

    as.data.frame(T_SUM)
}


f <- "bm-lasso.RData"
if (file.exists(f)) {
    load(f)
} else {
    bm <- run_benchmark()
    save(bm, file = f)
}
