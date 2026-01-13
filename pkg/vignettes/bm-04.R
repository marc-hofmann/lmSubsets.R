source("benchmark.R")


NAME <- "bm-04"


## Load existing benchmark.
##
bm <- load_benchmark(NAME)


## Initialize benchmark.
##
if (is.null(bm)) {


    ## New benchmark.
    ##
    ## Repetitions per case:  5
    ##
    bm <- benchmark(name = NAME, nrep = 5, seed = 1, timeout = 600)


    ## Setup distributions.
    ##
    ## DIST_ID |   SD | Distribution
    ## =============================
    ##       1 | 1.00 |       normal
    ##
    bm <- dist_benchmark(bm, SD = 1.00)


    ## Setup datasets.
    ##
    ## DATA_ID | NOBS | NVAR | NTRUE | DIST_ID | INTERCEPT
    ## ===================================================
    ##       1 | 1000 |  100 |    50 |       1 |       yes
    ## ---------------------------------------------------
    ##       2 | 1000 |  200 |   100 |       1 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 100, NTRUE =  50, DIST_ID = 1)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 200, NTRUE = 100, DIST_ID = 1)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT      | DATA_ID | IC  | NBEST
    ## ======================================================
    ##     1-5 | LM_SELECT |       1 | AIC | 1, 5, 10, 15, 20
    ## -----------------------------------------------------
    ##    6-10 | LM_SELECT |       2 | AIC | 1, 5, 10, 15, 20
    ## ------------------------------------------------------
    ##   11-15 | LM_SELECT |       1 | BIC | 1, 5, 10, 15, 20
    ## ------------------------------------------------------
    ##   16-20 | LM_SELECT |       2 | BIC | 1, 5, 10, 15, 20
    ##
    bm <- lmSelect_benchmark(bm, DATA_ID = 1, IC = "AIC", NBEST = c(1, 5, 10, 15, 20))
    #bm <- lmSelect_benchmark(bm, DATA_ID = 2, IC = "AIC", NBEST = c(1, 5, 10, 15, 20))
    bm <- lmSelect_benchmark(bm, DATA_ID = 1, IC = "BIC", NBEST = c(1, 5, 10, 15, 20))
    bm <- lmSelect_benchmark(bm, DATA_ID = 2, IC = "BIC", NBEST = c(1, 5, 10, 15, 20))


}


## Run benchmark.
##
bm <- run_benchmark(bm)


## Summarize benchmark.
##
bm <- summary_benchmark(bm)


## Save benchmark.
##
save_benchmark(bm)



###################
##  consolidate  ##
###################


report_benchmark <- function () {
    T_REPORT <- with(bm$DATA, data.frame(NVAR, DATA_ID = ID))

    T_REPORT <- merge(T_REPORT, with(bm$CASE, {
        data.frame(DATA_ID, CASE_ID = ID)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$LM_SELECT, {
        data.frame(CASE_ID, IC, NBEST)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$SUMMARY, {
        data.frame(CASE_ID, LM_SELECT = TM_AVG)
    }))

    with(T_REPORT, data.frame(NVAR, IC, NBEST, LM_SELECT))
}
