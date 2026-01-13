source("benchmark.R")


NAME <- "bm-01"


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
    bm <- benchmark(name = NAME, nrep = 5, seed = 1)


    ## Setup distributions.
    ##
    ## DIST_ID |   SD | Distribution
    ## =============================
    ##       1 | 0.05 |       normal
    ##       2 | 0.10 |       normal
    ##       3 | 0.50 |       normal
    ##       4 | 1.00 |       normal
    ##       5 | 5.00 |       normal
    ##
    bm <- dist_benchmark(bm, SD = c(0.05, 0.10, 0.50, 1.00, 5.00))


    ## Setup datasets.
    ##
    ## DATA_ID | NOBS | NVAR | NTRUE | DIST_ID | INTERCEPT
    ## ===================================================
    ##     1-5 | 1000 |   20 |    10 |     1-5 |       yes
    ## ---------------------------------------------------
    ##    6-10 | 1000 |   25 |    12 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   11-15 | 1000 |   30 |    15 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   16-20 | 1000 |   35 |    17 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   21-25 | 1000 |   40 |    20 |     1-5 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 20, NTRUE = 10, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 25, NTRUE = 12, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 30, NTRUE = 15, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 35, NTRUE = 17, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 40, NTRUE = 20, DIST_ID = 1:5)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID
    ## ==============================
    ##    1-25 | LM_SUBSETS |    1-25
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:25)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT  | DATA_ID | PREORD
    ## ==================================
    ##   26-50 | LEAPS |    1-25 |     no
    ## ----------------------------------
    ##   51-75 | LEAPS |    1-25 |    yes
    ##
    bm <- leaps_benchmark(bm, DATA_ID = 1:25, PREORD = FALSE)
    bm <- leaps_benchmark(bm, DATA_ID = 1:25, PREORD = TRUE )


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


check_benchmark <- function (tol = 1e-5) {
    T_CHECK <- do.call(rbind, {
        x <- with(bm$LM_SUBSETS,
                  data.frame(CASE_ID, WHAT = "LM_SUBSETS"))

        y <- with(subset(bm$LEAPS, !PREORD),
                  data.frame(CASE_ID, WHAT = "LEAPS"))

        z <- with(subset(bm$LEAPS, PREORD),
                  data.frame(CASE_ID, WHAT = "LEAPS1"))

        list(x, y, z)
    })

    T_CHECK <- merge(T_CHECK, with(bm$CASE, {
        data.frame(CASE_ID = ID, DATA_ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$RUN, {
        data.frame(CASE_ID, REP, RUN_ID = ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$VALUE, {
        data.frame(RUN_ID, RANK, RSS)
    }))

    T_CHECK <- {
        x <- with(subset(T_CHECK, WHAT == "LM_SUBSETS"),
                  data.frame(DATA_ID, REP, RANK,
                             LM_SUBSETS_CASE_ID = CASE_ID,
                             LM_SUBSETS_RSS = RSS))

        y <- with(subset(T_CHECK, WHAT == "LEAPS"),
                  data.frame(DATA_ID, REP, RANK,
                             LEAPS_CASE_ID = CASE_ID,
                             LEAPS_RSS = RSS))

        z <- with(subset(T_CHECK, WHAT == "LEAPS1"),
                  data.frame(DATA_ID, REP, RANK,
                             LEAPS1_CASE_ID = CASE_ID,
                             LEAPS1_RSS = RSS))

        merge(x, merge(y, z))
    }

    T_CHECK <- within(T_CHECK, {
        RSS_OK <- abs(LM_SUBSETS_RSS - LEAPS_RSS) <= tol
        RSS1_OK <- abs(LM_SUBSETS_RSS - LEAPS1_RSS) <= tol
    })

    with(T_CHECK, {
        stopifnot(all(RSS_OK))
        stopifnot(all(RSS1_OK))
    })

    invisible(T_CHECK)
}


report_benchmark <- function () {
    check_benchmark()

    T_REPORT <- with(bm$DIST, data.frame(SD, DIST_ID = ID))

    T_REPORT <- merge(T_REPORT, with(bm$DATA, {
        data.frame(DIST_ID, NVAR, DATA_ID = ID)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$CASE, {
        data.frame(DATA_ID, CASE_ID = ID)
    }))

    T_REPORT <- merge(T_REPORT, do.call(rbind, {
        x <- with(bm$LM_SUBSETS,
                  data.frame(CASE_ID, WHAT = "LM_SUBSETS"))

        y <- with(subset(bm$LEAPS, !PREORD),
                  data.frame(CASE_ID = CASE_ID, WHAT = "LEAPS"))

        z <- with(subset(bm$LEAPS, PREORD),
                  data.frame(CASE_ID = CASE_ID, WHAT = "LEAPS1"))

        list(x, y, z)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$SUMMARY, {
        data.frame(CASE_ID, TM_AVG)
    }))

    T_REPORT <- {
        x <- with(subset(T_REPORT, WHAT == "LM_SUBSETS"),
                  data.frame(SD, NVAR, LM_SUBSETS = TM_AVG))

        y <- with(subset(T_REPORT, WHAT == "LEAPS"),
                  data.frame(SD, NVAR, LEAPS = TM_AVG))

        z <- with(subset(T_REPORT, WHAT == "LEAPS1"),
                  data.frame(SD, NVAR, LEAPS1 = TM_AVG))

        merge(x, merge(y, z))
    }

    T_REPORT
}
