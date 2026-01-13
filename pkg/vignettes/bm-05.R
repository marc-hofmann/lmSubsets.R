source("benchmark.R")


NAME <- "bm-05"


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
    ##       1 | 1000 |   80 |    40 |       1 |       yes
    ## ---------------------------------------------------
    ##       2 | 1000 |  100 |    50 |       1 |       yes
    ## ---------------------------------------------------
    ##       3 | 1000 |  120 |    60 |       1 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 1000, NVAR =  80, NTRUE =  40, DIST_ID = 1)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 100, NTRUE =  50, DIST_ID = 1)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 120, NTRUE =  60, DIST_ID = 1)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT      | DATA_ID | IC
    ## ===============================================================
    ##     1-5 | LM_SELECT |       1 |  1.0, 2.0, 4.0, 8.0, 16.0, 32.0
    ## ---------------------------------------------------------------
    ##    6-10 | LM_SELECT |       2 |  1.0, 2.0, 4.0, 8.0, 16.0, 32.0
    ## ---------------------------------------------------------------
    ##   11-15 | LM_SELECT |       3 |  1.0, 2.0, 4.0, 8.0, 16.0, 32.0
    ##
    bm <- lmSelect_benchmark(bm, DATA_ID = 1, IC = c(1, 2, 4, 8, 16, 32))
    bm <- lmSelect_benchmark(bm, DATA_ID = 2, IC = c(1, 2, 4, 8, 16, 32))
    bm <- lmSelect_benchmark(bm, DATA_ID = 3, IC = c(1, 2, 4, 8, 16, 32))


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
        data.frame(CASE_ID, IC)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$SUMMARY, {
        data.frame(CASE_ID, LM_SELECT = TM_AVG)
    }))

    with(T_REPORT, data.frame(NVAR, IC, LM_SELECT))
}
