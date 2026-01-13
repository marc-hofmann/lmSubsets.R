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
    dca <- lmSubsets(mortality ~ ., data = AirPollution, .algo = "dca")
    stopifnot(dca$.nodes == 16384)

    ## BBA
    bba <- lmSubsets(mortality ~ ., data = AirPollution, pradius = 0,
                     .algo = "bba")
    stopifnot(bba$.nodes == 712)

    ## BBA1
    bba1 <- lmSubsets(mortality ~ ., data = AirPollution, pradius = 1,
                      .algo = "bba")
    stopifnot(bba1$.nodes == 381)

    ## BBA15
    bba15 <- lmSubsets(mortality ~ ., data = AirPollution,
                       pradius = 15, .algo = "bba")
    stopifnot(bba15$.nodes == 146)

    ## HBBA
    hbba <- lmSubsets(mortality ~ ., data = AirPollution,
                      tolerance = 0.1, pradius = 0, .algo = "hbba")
    stopifnot(hbba$.nodes == 353)

    ## HBBA1
    hbba1 <- lmSubsets(mortality ~ ., data = AirPollution,
                       tolerance = 0.1, pradius = 1, .algo = "hbba")
    stopifnot(hbba1$.nodes == 124)

    ## ABBA
    abba <- lmSubsets(mortality ~ ., data = AirPollution,
                      tolerance = 0.1, pradius = 0, .algo = "abba")
    stopifnot(abba$.nodes == 612)

    ## ABBA1
    abba1 <- lmSubsets(mortality ~ ., data = AirPollution,
                       tolerance = 0.1, pradius = 1, .algo = "abba")
    stopifnot(abba1$.nodes == 295)
})



## canonical case
local({
    message ("Canonical case...")

    all <- lmSubsets(mortality ~ ., data = AirPollution)

    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    rf5 <- refit(all, size = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)          ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(all, size = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)          ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(all, size = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)          ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(all, size = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)          ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(all, size = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)          ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(all, size = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)          ) < TOL)
    stopifnot(abs(coef(lm5) - coef(all, size = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)          ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(all, size = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)          ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(all, size = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)          ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(all, size = 5)) < TOL)
})



## no intercept
local({
    message ("No intercept...")

    all <- lmSubsets(mortality ~ . - 1, data = AirPollution)

    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5)]
    lm5 <- lm(mortality ~ . - 1, data = lm5)

    rf5 <- refit(all, size = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)          ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(all, size = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)          ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(all, size = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)          ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(all, size = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)          ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(all, size = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)          ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(all, size = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)          ) < TOL)
    stopifnot(abs(coef(lm5) - coef(all, size = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)          ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(all, size = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)          ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(all, size = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)          ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(all, size = 5)) < TOL)
})



## weights and offset
local({
    message ("Weights and offset...")

    w <- abs(rnorm(nall))
    o <- rnorm(nall)

    all <- lmSubsets(mortality ~ ., data = AirPollution,
                     weights = w, offset = o)

    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5, weights = w, offset = o)

    rf5 <- refit(all, size = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)          ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(all, size = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)          ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(all, size = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)          ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(all, size = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)          ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(all, size = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)          ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(all, size = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)          ) < TOL)
    stopifnot(abs(coef(lm5) - coef(all, size = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)          ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(all, size = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)          ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(all, size = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)          ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(all, size = 5)) < TOL)
})



## zero weights
local({
    message ("Zero weights...")

    w <- rep_len(1, nall)
    w[sample(nall, 10)] <- 0

    all <- lmSubsets(mortality ~ ., data = AirPollution, weights = w)

    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5, weights = w)

    rf5 <- refit(all, size = 5)

    ## deviance
    stopifnot(abs(deviance(lm5) - deviance(rf5)          ) < TOL)
    stopifnot(abs(deviance(lm5) - deviance(all, size = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5) - sigma(rf5)          ) < TOL)
    stopifnot(abs(sigma(lm5) - sigma(all, size = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5) - logLik(rf5)          ) < TOL)
    stopifnot(abs(logLik(lm5) - logLik(all, size = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5) - AIC(rf5)          ) < TOL)
    stopifnot(abs(AIC(lm5) - AIC(all, size = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5) - BIC(rf5)          ) < TOL)
    stopifnot(abs(BIC(lm5) - BIC(all, size = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5) - coef(rf5)          ) < TOL)
    stopifnot(abs(coef(lm5) - coef(all, size = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5) - vcov(rf5)          ) < TOL)
    stopifnot(abs(vcov(lm5) - vcov(all, size = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5) - fitted(rf5)          ) < TOL)
    stopifnot(abs(fitted(lm5) - fitted(all, size = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5) - residuals(rf5)          ) < TOL)
    stopifnot(abs(residuals(lm5) - residuals(all, size = 5)) < TOL)
})



## matrix interface
local({
    message ("Matrix interface...")

    all <- lmSubsets(mortality ~ ., data = AirPollution)

    mat <- lmSubsets(as.matrix(AirPollution), y = "mortality")

    lm5 <- variable.names(mat, size = 5)
    lm5 <- paste0("mortality ~ ", paste(lm5[-1], collapse = " + "))
    lm5 <- lm(as.formula(lm5), data = AirPollution)

    rf5 <- refit(mat, size = 5)

    ## deviance
    stopifnot(abs(deviance(lm5)           - deviance(all, size = 5)) < TOL)
    stopifnot(abs(deviance(rf5)           - deviance(all, size = 5)) < TOL)
    stopifnot(abs(deviance(mat, size = 5) - deviance(all, size = 5)) < TOL)

    ## sigma
    stopifnot(abs(sigma(lm5)           - sigma(all, size = 5)) < TOL)
    stopifnot(abs(sigma(rf5)           - sigma(all, size = 5)) < TOL)
    stopifnot(abs(sigma(mat, size = 5) - sigma(all, size = 5)) < TOL)

    ## logLik
    stopifnot(abs(logLik(lm5)           - logLik(all, size = 5)) < TOL)
    stopifnot(abs(logLik(rf5)           - logLik(all, size = 5)) < TOL)
    stopifnot(abs(logLik(mat, size = 5) - logLik(all, size = 5)) < TOL)

    ## AIC
    stopifnot(abs(AIC(lm5)           - AIC(all, size = 5)) < TOL)
    stopifnot(abs(AIC(rf5)           - AIC(all, size = 5)) < TOL)
    stopifnot(abs(AIC(mat, size = 5) - AIC(all, size = 5)) < TOL)

    ## BIC
    stopifnot(abs(BIC(lm5)           - BIC(all, size = 5)) < TOL)
    stopifnot(abs(BIC(rf5)           - BIC(all, size = 5)) < TOL)
    stopifnot(abs(BIC(mat, size = 5) - BIC(all, size = 5)) < TOL)

    ## coef
    stopifnot(abs(coef(lm5)           - coef(all, size = 5)) < TOL)
    stopifnot(abs(coef(rf5)           - coef(all, size = 5)) < TOL)
    stopifnot(abs(coef(mat, size = 5) - coef(all, size = 5)) < TOL)

    ## vcov
    stopifnot(abs(vcov(lm5)           - vcov(all, size = 5)) < TOL)
    stopifnot(abs(vcov(rf5)           - vcov(all, size = 5)) < TOL)
    stopifnot(abs(vcov(mat, size = 5) - vcov(all, size = 5)) < TOL)

    ## fitted
    stopifnot(abs(fitted(lm5)           - fitted(all, size = 5)) < TOL)
    stopifnot(abs(fitted(rf5)           - fitted(all, size = 5)) < TOL)
    stopifnot(abs(fitted(mat, size = 5) - fitted(all, size = 5)) < TOL)

    ## residuals
    stopifnot(abs(residuals(lm5)           - residuals(all, size = 5)) < TOL)
    stopifnot(abs(residuals(rf5)           - residuals(all, size = 5)) < TOL)
    stopifnot(abs(residuals(mat, size = 5) - residuals(all, size = 5)) < TOL)
})



## NSE 1
local({
    message ("Non-standard evaluation I...")

    lm0 <- lm(mortality ~ ., data = AirPollution)

    all <- lmSubsets(lm0)
    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    nse <- lmSubsets(mortality ~ ., data = AirPollution, model = FALSE)
    rf5 <- refit(nse, size = 5)

    ## model frame
    stopifnot(abs(model.frame(nse)           - model.frame(lm0)) < TOL)
    stopifnot(abs(model.frame(rf5)           - model.frame(lm5)) < TOL)
})



## NSE 2
local({
    message ("Non-standard evaluation II...")

    lm0 <- lm(mortality ~ ., data = AirPollution)

    all <- lmSubsets(lm0)
    lm5 <- variable.names(all, size = 5)
    lm5 <- AirPollution[c("mortality", lm5[-1])]
    lm5 <- lm(mortality ~ ., data = lm5)

    nse <- paste(names(AirPollution)[-16], collapse = " + ")
    nse <- paste0("mortality ~ ", nse)
    nse <- lmSubsets(as.formula(nse), model = FALSE)
    rf5 <- refit(nse, size = 5)

    ## model frame
    stopifnot(abs(model.frame(nse)           - model.frame(lm0)) < TOL)
    stopifnot(abs(model.frame(rf5)           - model.frame(lm5)) < TOL)
}, env = AirPollution)
