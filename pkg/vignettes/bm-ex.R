source("benchmark.R")


## benchmark name
NAME <- "bm-ex"


## cases to execute
## e.g.  EXEC_CASES = 1:3
EXEC_CASES <- 1:16


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
    ## DIST_ID | MEAN |  SD | RATE | MEANLOG | SDLOG | Distribution
    ## ============================================================
    ##       1 |  0.0 | 1.0 |    - |       - |     - |       normal
    ## ------------------------------------------------------------
    ##       2 |    - |   - |  1.0 |       - |     - |  exponential
    ##
    bm <- dist_benchmark(bm, SD = 1.0)
    bm <- dist_benchmark(bm, RATE = 1.0)


    ## Setup datasets.
    ##
    ## DATA_ID | NOBS | NVAR | NTRUE | DIST_ID | INTERCEPT
    ## ===================================================
    ##       1 |   25 |   20 |    10 |       1 |       yes
    ##       2 |   25 |   20 |    10 |       2 |       yes
    ## ---------------------------------------------------
    ##       3 |   50 |   40 |    20 |       1 |       yes
    ##       4 |   50 |   40 |    20 |       2 |       yes
    ##
    bm <- data_benchmark(bm, NOBS = 25, NVAR = 20, NTRUE = 10, DIST_ID = 1:2)
    bm <- data_benchmark(bm, NOBS = 50, NVAR = 40, NTRUE = 20, DIST_ID = 1:2)


    ## Setup cases.
    ##
    ## CASE_ID | WHAT       | DATA_ID | TOLERANCE | NMIN | NMAX | NBEST
    ## ================================================================
    ##       1 | LM_SUBSETS |       1 |       0.0 |    - |    - |     -
    ##       2 | LM_SUBSETS |       2 |       0.0 |    - |    - |     -
    ##       3 | LM_SUBSETS |       3 |       0.0 |    - |    - |     -
    ##       4 | LM_SUBSETS |       4 |       0.0 |    - |    - |     -
    ## ----------------------------------------------------------------
    ##       5 | LM_SUBSETS |       1 |       0.1 |    - |    - |     -
    ##       6 | LM_SUBSETS |       2 |       0.1 |    - |    - |     -
    ##       7 | LM_SUBSETS |       3 |       0.1 |    - |    - |     -
    ##       8 | LM_SUBSETS |       4 |       0.1 |    - |    - |     -
    ## ------------------------------------------ ---------------------
    ##       9 | LM_SELECT  |       1 |       0.0 |    - |    - |     -
    ##      10 | LM_SELECT  |       2 |       0.0 |    - |    - |     -
    ##      11 | LM_SELECT  |       3 |       0.0 |    - |    - |     -
    ##      12 | LM_SELECT  |       4 |       0.0 |    - |    - |     -
    ## ------------------------------------------ ---------------------
    ##      13 | LM_SELECT  |       1 |       0.1 |    - |    - |     -
    ##      14 | LM_SELECT  |       2 |       0.1 |    - |    - |     -
    ##      15 | LM_SELECT  |       3 |       0.1 |    - |    - |     -
    ##      16 | LM_SELECT  |       4 |       0.1 |    - |    - |     -
    ##
    ## Note:  Number of repetitions per case: nrep = 5
    ##
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:4, TOLERANCE = 0.0)
    bm <- lmSubsets_benchmark(bm, DATA_ID = 1:4, TOLERANCE = 0.1)
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:4, TOLERANCE = 0.0)
    bm <- lmSelect_benchmark(bm, DATA_ID = 1:4, TOLERANCE = 0.1)


}


## Run benchmark.
##
bm <- run_benchmark(bm, cases = EXEC_CASES)


## Summarize benchmark.
##
bm <- summary_benchmark(bm)


## Save benchmark.
##
save_benchmark(bm)
