source("benchmark.R")


NAME <- "bm-02"


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
    ## CASE_ID | WHAT      | DATA_ID | IC
    ## ===================================
    ##    1-25 | LM_SELECT |    1-25 | AIC
    ## -----------------------------------
    ##   26-50 | LM_SELECT |    1-25 | BIC
    ##
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:25, IC = "AIC")
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:25, IC = "BIC")


    ## Setup cases.
    ##
    ## CASE_ID | WHAT    | DATA_ID | IC  | PREORD
    ## ==========================================
    ##   51-75 | BESTGLM |    1-25 | AIC |     no
    ## ------------------------------------------
    ##  76-100 | BESTGLM |    1-25 | AIC |    yes
    ## ------------------------------------------
    ## 101-125 | BESTGLM |    1-25 | BIC |     no
    ## ------------------------------------------
    ## 126-150 | BESTGLM |    1-25 | BIC |    yes
    ##
    bm <- bestglm_benchmark(bm, DATA_ID = 1:25, IC = "AIC", PREORD = FALSE)
    bm <- bestglm_benchmark(bm, DATA_ID = 1:25, IC = "AIC", PREORD = TRUE )
    bm <- bestglm_benchmark(bm, DATA_ID = 1:25, IC = "BIC", PREORD = FALSE)
    bm <- bestglm_benchmark(bm, DATA_ID = 1:25, IC = "BIC", PREORD = TRUE )


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
        x <- with(bm$LM_SELECT,
                  data.frame(CASE_ID, IC, WHAT = "LM_SELECT"))

        y <- with(subset(bm$BESTGLM, !PREORD),
                  data.frame(CASE_ID, IC, WHAT = "BESTGLM"))

        z <- with(subset(bm$BESTGLM, PREORD),
                  data.frame(CASE_ID, IC, WHAT = "BESTGLM1"))

        list(x, y, z)
    })

    T_CHECK <- merge(T_CHECK, with(bm$CASE, {
        data.frame(CASE_ID = ID, DATA_ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$RUN, {
        data.frame(CASE_ID, REP, RUN_ID = ID)
    }))

    T_CHECK <- merge(T_CHECK, with(bm$VALUE, {
        data.frame(RUN_ID, AIC, BIC)
    }))

    T_CHECK <- {
        x <- with(subset(T_CHECK, WHAT == "LM_SELECT"),
                  data.frame(DATA_ID, REP, IC,
                             LM_SELECT_CASE_ID = CASE_ID,
                             LM_SELECT_AIC = AIC,
                             LM_SELECT_BIC = BIC))

        y <- with(subset(T_CHECK, WHAT == "BESTGLM"),
                  data.frame(DATA_ID, REP, IC, BESTGLM_CASE_ID = CASE_ID,
                             BESTGLM_AIC = AIC, BESTGLM_BIC = BIC))

        z <- with(subset(T_CHECK, WHAT == "BESTGLM1"),
                  data.frame(DATA_ID, REP, IC, BESTGLM1_CASE_ID = CASE_ID,
                             BESTGLM1_AIC = AIC, BESTGLM1_BIC = BIC))

        merge(x, merge(y, z))
    }

    T_CHECK <- within(T_CHECK, {
        AIC_OK <- abs(LM_SELECT_AIC - BESTGLM_AIC) <= tol
        AIC1_OK <- abs(LM_SELECT_AIC - BESTGLM1_AIC) <= tol

        BIC_OK <- abs(LM_SELECT_BIC - BESTGLM_BIC) <= tol
        BIC1_OK <- abs(LM_SELECT_BIC - BESTGLM1_BIC) <= tol
    })

    with(T_CHECK, {
        stopifnot(all(AIC_OK))
        stopifnot(all(AIC1_OK))

        stopifnot(all(BIC_OK))
        stopifnot(all(BIC1_OK))
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
        x <- with(bm$LM_SELECT,
                  data.frame(CASE_ID, IC, WHAT = "LM_SELECT"))

        y <- with(subset(bm$BESTGLM, !PREORD),
                  data.frame(CASE_ID = CASE_ID, IC, WHAT = "BESTGLM"))

        z <- with(subset(bm$BESTGLM, PREORD),
                  data.frame(CASE_ID = CASE_ID, IC, WHAT = "BESTGLM1"))

        list(x, y, z)
    }))

    T_REPORT <- merge(T_REPORT, with(bm$SUMMARY, {
        data.frame(CASE_ID, TM_AVG)
    }))

    T_REPORT <- {
        x <- with(subset(T_REPORT, WHAT == "LM_SELECT"),
                  data.frame(IC, SD, NVAR, LM_SELECT = TM_AVG))

        y <- with(subset(T_REPORT, WHAT == "BESTGLM"),
                  data.frame(IC, SD, NVAR, BESTGLM = TM_AVG))

        z <- with(subset(T_REPORT, WHAT == "BESTGLM1"),
                  data.frame(IC, SD, NVAR, BESTGLM1 = TM_AVG))

        merge(x, merge(y, z))
    }

    T_REPORT
}
