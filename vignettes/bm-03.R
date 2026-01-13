source("benchmark.R")


NAME <- "bm-03"


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
    ##    6-10 | 1000 |   40 |    20 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   11-15 | 1000 |   60 |    30 |     1-5 |       yes
    ## ---------------------------------------------------
    ##   16-20 | 1000 |   80 |    40 |     1-5 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 20, NTRUE = 10, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 40, NTRUE = 20, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 60, NTRUE = 30, DIST_ID = 1:5)
    bm <- data_benchmark(bm, NOBS = 1000, NVAR = 80, NTRUE = 40, DIST_ID = 1:5)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE
    ## ==========================================
    ##    1-15 | LM_SUBSETS |    1-15 |       0.0
    ## ------------------------------------------
    ##   16-30 | LM_SUBSETS |    1-15 |       0.1
    ## ------------------------------------------
    ##   31-45 | LM_SELECT  |    1-15 |       0.0
    ## ------------------------------------------
    ##   46-60 | LM_SELECT  |    1-15 |       0.1
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.0)
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.1)
    bm <-  lmSelect_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.0)
    bm <-  lmSelect_benchmark(bm, DATA_ID = 1:15, TOLERANCE = 0.1)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX
    ## ========================================================
    ##   61-65 | LM_SUBSETS |   16-20 |       0.0 |   30 |   50
    ## --------------------------------------------------------
    ##   66-70 | LM_SUBSETS |   16-20 |       0.1 |   30 |   50
    ## --------------------------------------------------------
    ##   71-75 | LM_SELECT  |   16-20 |       0.0 |    - |    -
    ## --------------------------------------------------------
    ##   76-80 | LM_SELECT  |   16-20 |       0.1 |    - |    -
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.0, NMIN = 30, NMAX = 50)
    bm <- lmSubsets_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.1, NMIN = 30, NMAX = 50)
    bm <-  lmSelect_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.0)
    bm <-  lmSelect_benchmark(bm, DATA_ID = 16:20, TOLERANCE = 0.1)


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


check0_benchmark <- function (tol = 1e-5) {
    T_CHECK <- do.call(rbind, {
        x <- with(subset(bm$LM_SUBSETS, TOLERANCE == 0.0),
                  data.frame(CASE_ID, WHAT = "LM_SUBSETS"))

        y <- with(subset(bm$LM_SELECT, TOLERANCE == 0.0),
                  data.frame(CASE_ID, WHAT = "LM_SELECT"))

        list(x, y)
    })

    T_CHECK <- merge(T_CHECK, with(bm$CASE, {
        data.frame(CASE_ID = ID, DATA_ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$RUN, {
        data.frame(CASE_ID, REP, RUN_ID = ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$VALUE, {
        data.frame(RUN_ID, RANK, BIC)
    }))

    T_CHECK <- {
        x <- with(subset(T_CHECK, WHAT == "LM_SUBSETS"),
                  data.frame(DATA_ID, REP,
                             LM_SUBSETS_CASE_ID = CASE_ID,
                             LM_SUBSETS_RANK = RANK,
                             LM_SUBSETS_BIC = BIC))

        y <- with(subset(T_CHECK, WHAT == "LM_SELECT"),
                  data.frame(DATA_ID, REP,
                             LM_SELECT_CASE_ID = CASE_ID,
                             LM_SELECT_RANK = RANK,
                             LM_SELECT_BIC = BIC))

        merge(x, y)
    }

    T_CHECK <- within(T_CHECK, {
        OK <- (LM_SELECT_BIC - LM_SUBSETS_BIC) <= tol
    })

    with(T_CHECK, {
        stopifnot(all(OK))
    })

    invisible(T_CHECK)
}


check1_benchmark <- function (tol = 1e-5) {
    T_CHECK <- with(bm$LM_SUBSETS, data.frame(CASE_ID, TOLERANCE))

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
        x <- with(subset(T_CHECK, TOLERANCE > 0.0),
                  data.frame(DATA_ID, REP, RANK,
                             CASE_ID, RSS, TOLERANCE))

        y <- with(subset(T_CHECK, TOLERANCE == 0.0),
                  data.frame(DATA_ID, REP, RANK,
                             CASE0_ID = CASE_ID, RSS0 = RSS))

        merge(x, y)
    }

    T_CHECK <- merge(T_CHECK, with(bm$SIMU, {
        data.frame(DATA_ID, REP, RSS_FUL)
    }))

    T_CHECK <- within(T_CHECK, {
        OK <- ((RSS - RSS_FUL) - tol) <= ((1 + TOLERANCE) * (RSS0 - RSS_FUL))
    })

    stopifnot(with(T_CHECK, all(OK)))

    invisible(T_CHECK)
}


check2_benchmark <- function (tol = 1e-5) {
    T_CHECK <- with(bm$LM_SELECT, data.frame(CASE_ID, TOLERANCE))

    T_CHECK <- merge(T_CHECK, with(bm$CASE, {
        data.frame(CASE_ID = ID, DATA_ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$RUN, {
        data.frame(CASE_ID, REP, RUN_ID = ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$VALUE, {
        data.frame(RUN_ID, BIC)
    }))

    T_CHECK <- {
        x <- with(subset(T_CHECK, TOLERANCE > 0.0),
                  data.frame(DATA_ID, REP,
                             CASE_ID, BIC, TOLERANCE))

        y <- with(subset(T_CHECK, TOLERANCE == 0.0),
                  data.frame(DATA_ID, REP,
                             CASE0_ID = CASE_ID, BIC0 = BIC))

        merge(x, y)
    }

    T_CHECK <- merge(T_CHECK, with(bm$SIMU, {
        data.frame(DATA_ID, REP, BIC_INF)
    }))

    T_CHECK <- within(T_CHECK, {
        OK <- ((BIC - BIC_INF) - tol) <= ((1 + TOLERANCE) * (BIC0 - BIC_INF))
    })

    with(T_CHECK, {
        stopifnot(all(OK))
    })

    invisible(T_CHECK)
}


check_benchmark <- function (tol = 1e-5) {
    check0_benchmark(tol);
    check1_benchmark(tol);
    check2_benchmark(tol);

    invisible()
}


report_benchmark <- function () {
    check_benchmark()

    T_REPORT <- with(bm$DIST, data.frame(SD, DIST_ID = ID))

    T_REPORT <- merge(T_REPORT, with(bm$DATA, {
        data.frame(DIST_ID, DATA_ID = ID, NVAR)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$CASE, {
        data.frame(DATA_ID, CASE_ID = ID)
    }))

    T_REPORT <- merge(T_REPORT, do.call(rbind, {
        x <- with(bm$LM_SUBSETS,
                  data.frame(CASE_ID, TOLERANCE, NMIN, NMAX,
                             WHAT = "LM_SUBSETS"))

        y <- with(bm$LM_SELECT,
                  data.frame(CASE_ID, TOLERANCE, NMIN = NA, NMAX = NA,
                             WHAT = "LM_SELECT"))

        list(x, y)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$SUMMARY, {
        data.frame(CASE_ID, TM_AVG)
    }))

    T_REPORT <- {
        x <- with(subset(T_REPORT, WHAT == "LM_SUBSETS"),
                  data.frame(SD, NVAR, TOLERANCE, NMIN, NMAX,
                             LM_SUBSETS = TM_AVG))

        y <- with(subset(T_REPORT, WHAT == "LM_SELECT"),
                  data.frame(SD, NVAR, TOLERANCE,
                             LM_SELECT = TM_AVG))

        merge(x, y)
    }

    T_REPORT
}
