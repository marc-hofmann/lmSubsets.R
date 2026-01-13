TOL <- 5e-10

## load package
library(lmSubsets)

## load data
data(AirPollution)

nall <- nrow(AirPollution)



## fundamental algorithms
local({
    message ("Fundamental algorithms...")

    ## DCA
    dca <- lmSelect(mortality ~ ., data = AirPollution, .algo = "dca")
    stopifnot(dca$.nodes == 16384)

    ## BBA
    bba <- lmSelect(mortality ~ ., data = AirPollution, pradius = 0,
                    .algo = "bba")
    stopifnot(bba$.nodes == 604)

    ## BBA1
    bba1 <- lmSelect(mortality ~ ., data = AirPollution, pradius = 1,
                     .algo = "bba")
    stopifnot(bba1$.nodes == 313)

    ## ABBA
    abba <- lmSelect(mortality ~ ., data = AirPollution,
                     tolerance = 0.1, pradius = 0, .algo = "abba")
    stopifnot(abba$.nodes == 411)

    ## ABBA1
    abba1 <- lmSelect(mortality ~ ., data = AirPollution,
                      tolerance = 0.1, pradius = 1, .algo = "abba")
    stopifnot(abba1$.nodes == 169)
})



## canonical case
local({
    message ("Canonical case...")

    best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10)

    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    rf5 <- refit(best, best = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)           ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(best, best = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)           ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(best, best = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)           ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(best, best = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)           ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(best, best = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)           ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(best, best = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)           ) < TOL)
    stopifnot(abs(coef(lm5) - coef(best, best = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)           ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(best, best = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)           ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(best, best = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)           ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(best, best = 5)) < TOL)
})



## coercion
local({
    message ("Coercion...")

    all <- lmSubsets(mortality ~ ., data = AirPollution, nbest = 10)

    best0 <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10)
    best1 <- lmSelect(all)

    stopifnot(abs(with(best0$submodel, IC) - with(best1$submodel, IC)) < TOL)
})



## no intercept
local({
    message ("No intercept...")

    best <- lmSelect(mortality ~ . - 1, data = AirPollution, nbest = 10)

    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5)]
    lm5 <- lm(mortality ~ . - 1, data = lm5)

    rf5 <- refit(best, best = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)           ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(best, best = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)           ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(best, best = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)           ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(best, best = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)           ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(best, best = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)           ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(best, best = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)           ) < TOL)
    stopifnot(abs(coef(lm5) - coef(best, best = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)           ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(best, best = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)           ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(best, best = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)           ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(best, best = 5)) < TOL)
})



## weights and offset
local({
    message ("Weights and offset...")

    w <- abs(rnorm(nall))
    o <- rnorm(nall)

    best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10,
                     weights = w, offset = o)

    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5, weights = w, offset = o)

    rf5 <- refit(best, best = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)           ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(best, best = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)           ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(best, best = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)           ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(best, best = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)           ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(best, best = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)           ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(best, best = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)           ) < TOL)
    stopifnot(abs(coef(lm5) - coef(best, best = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)           ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(best, best = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)           ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(best, best = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)           ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(best, best = 5)) < TOL)
})



## zero weights
local({
    message ("Zero weights...")

    w <- rep_len(1, nall)
    w[sample(nall, 10)] <- 0

    best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10,
                     weights = w)

    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5, weights = w)

    rf5 <- refit(best, best = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)           ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(best, best = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)           ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(best, best = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)           ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(best, best = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)           ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(best, best = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)           ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(best, best = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)           ) < TOL)
    stopifnot(abs(coef(lm5) - coef(best, best = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)           ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(best, best = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)           ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(best, best = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)           ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(best, best = 5)) < TOL)
})



## matrix interface
local({
    message ("Matrix interface...")

    best <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10)

    mat <- lmSelect(as.matrix(AirPollution), y = "mortality", nbest = 10)

    lm5 <- variable.names(mat, best = 5)
    lm5 <- paste0("mortality ~ ", paste(lm5[-1], collapse = " + "))
    lm5 <- lm(as.formula(lm5), data = AirPollution)

    rf5 <- refit(mat, best = 5)

    ## deviance
    stopifnot(abs(deviance(lm5)           - deviance(best, best = 5)) < TOL)
    stopifnot(abs(deviance(rf5)           - deviance(best, best = 5)) < TOL)
    stopifnot(abs(deviance(mat, best = 5) - deviance(best, best = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5)           - sigma(best, best = 5)) < TOL)
    stopifnot(abs(sigma(rf5)           - sigma(best, best = 5)) < TOL)
    stopifnot(abs(sigma(mat, best = 5) - sigma(best, best = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5)           - logLik(best, best = 5)) < TOL)
    stopifnot(abs(logLik(rf5)           - logLik(best, best = 5)) < TOL)
    stopifnot(abs(logLik(mat, best = 5) - logLik(best, best = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5)           - AIC(best, best = 5)) < TOL)
    stopifnot(abs(AIC(rf5)           - AIC(best, best = 5)) < TOL)
    stopifnot(abs(AIC(mat, best = 5) - AIC(best, best = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5)           - BIC(best, best = 5)) < TOL)
    stopifnot(abs(BIC(rf5)           - BIC(best, best = 5)) < TOL)
    stopifnot(abs(BIC(mat, best = 5) - BIC(best, best = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5)           - coef(best, best = 5)) < TOL)
    stopifnot(abs(coef(rf5)           - coef(best, best = 5)) < TOL)
    stopifnot(abs(coef(mat, best = 5) - coef(best, best = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5)           - vcov(best, best = 5)) < TOL)
    stopifnot(abs(vcov(rf5)           - vcov(best, best = 5)) < TOL)
    stopifnot(abs(vcov(mat, best = 5) - vcov(best, best = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5)           - fitted(best, best = 5)) < TOL)
    stopifnot(abs(fitted(rf5)           - fitted(best, best = 5)) < TOL)
    stopifnot(abs(fitted(mat, best = 5) - fitted(best, best = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5)           - residuals(best, best = 5)) < TOL)
    stopifnot(abs(residuals(rf5)           - residuals(best, best = 5)) < TOL)
    stopifnot(abs(residuals(mat, best = 5) - residuals(best, best = 5)) < TOL)
})



## custom criterion
local({
    message ("Custom criterion...")

    k <- 2

    ll <- function (rss) {
        -nall/2 * (log(2 * pi) - log(nall) + log(rss) + 1)
    }

    aic <- function (size, rss) {
        -2 * ll(rss) + k * (size + 1)
    }

    best_aic <- lmSelect(mortality ~ ., data = AirPollution, penalty = "AIC",
                         .algo = "abba")
    best_fun <- lmSelect(mortality ~ ., data = AirPollution, penalty = aic,
                         .algo = "abba")

    stopifnot(best_aic$.nodes == best_fun$.nodes)  # 137 nodes
    stopifnot(abs(with(best_aic$submodel, IC) - with(best_fun$submodel, IC))
              < TOL)
})



## ABBA/HBBA
local({
    message ("ABBA vs. HBBA...")

    k <- 2

    ll <- function (rss) {
        -nall/2 * (log(2 * pi) - log(nall) + log(rss) + 1)
    }

    aic <- function (size, rss) {
        -2 * ll(rss) + k * (size + 1)
    }

    aic0 <- aic(16, 53680.02)

    fun <- function (size, rss) {
        aic(size, rss) - aic0
    }

    abba <- lmSelect(mortality ~ ., data = AirPollution, penalty = "AIC",
                     .algo = "abba")
    hbba <- lmSelect(mortality ~ ., data = AirPollution, penalty = fun,
                     .algo = "hbba")

    stopifnot(abba$.nodes == hbba$.nodes)  # 137 nodes
    stopifnot(abs(with(abba$submodel, IC) - with(hbba$submodel, IC + aic0))
              < TOL)
})



## NSE 1
local({
    message ("Non-standard evaluation I...")

    lm0 <- lm(mortality ~ ., data = AirPollution)

    best <- lmSelect(lm0, nbest = 10)
    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    nse <- lmSelect(mortality ~ ., data = AirPollution, nbest = 10,
                    model = FALSE)
    rf5 <- refit(nse, best = 5)

    ## model frame
    stopifnot(abs(model.frame(nse)           - model.frame(lm0)) < TOL)
    stopifnot(abs(model.frame(rf5)           - model.frame(lm5)) < TOL)
})



## NSE 2
local({
    message ("Non-standard evaluation II...")

    lm0 <- lm(mortality ~ ., data = AirPollution)

    best <- lmSelect(lm0, nbest = 10)
    lm5 <- variable.names(best, best = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    nse <- paste(names(AirPollution)[-16], collapse = " + ")
    nse <- paste0("mortality ~ ", nse)
    nse <- lmSelect(as.formula(nse), nbest = 10, model = FALSE)
    rf5 <- refit(nse, best = 5)

    ## model frame
    stopifnot(abs(model.frame(nse)           - model.frame(lm0)) < TOL)
    stopifnot(abs(model.frame(rf5)           - model.frame(lm5)) < TOL)
}, env = AirPollution)
