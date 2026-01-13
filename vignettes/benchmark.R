library("parallel")
library("tools")

library("lmSubsets")
library("leaps")
library("bestglm")



data_frame <- function (...) {
    data.frame(..., stringsAsFactors = FALSE)
}


expand_frame <- function (...) {
    args <- rev(list(...))
    args$stringsAsFactors <- FALSE

    ans <- do.call(expand.grid, args)
    rev(ans)
}


frame_id <- function (df_new, df) {
    ans <- seq_len(NROW(df_new))

    if (is.null(df))  return (ans)

    ans + df[NROW(df), "ID"]
}


preorder_data <- function (x, y, intercept) {
    n <- ncol(x)

    xy <- cbind(x, y)
    if (intercept)  xy <- cbind(1, xy)

    rss <- numeric(n)
    for (j in seq_len(n)) {
        rz <- qr.R(qr(xy[, -(j + intercept)]))
        rss[j] <- rz[n + intercept, n + intercept]^2
    }

    x[, order(rss, decreasing = TRUE)]
}



## Initialize benchmark.
##
## Arguments:
##   name    - (character)
##   nrep    - (integer) number of repetitions
##   timeout - (integer) timeout in seconds
##   seed    - (integer)
##
## Result: (list)
##   benchmark
##
benchmark <- function (name = "bm", nrep = 1, timeout = NULL,
                       seed = NULL) {
    bm <- list()

    bm$name <- name
    bm$nrep <- nrep
    bm$timeout <- timeout
    bm$seed <- if (!is.null(seed)) seed else as.numeric(Sys.time())

    bm
}


## Save benchmark.
##
## Arguments:
##   bm   - (list) benchmark
##   name - (character) benchmark name
##
## Result: (list)
##   benchmark
##
save_benchmark <- function (bm, name = NULL) {
    if (is.null(name))  name <- bm$name

    file <- paste0(bm$name, ".RData")
    save(bm, file = file)

    invisible(bm)
}


## Load benchmark.
##
## Arguments:
##   name - (character) benchmark name
##
## Result: (list)
##   benchmark
##
load_benchmark <- function (name) {
    file <- paste0(name, ".RData")

    if (!file.exists(file)) {
        return (NULL)
    }

    load(file = file)

    bm
}


## Remove benchmark.
##
## Arguments:
##   bm   - (list) benchmark
##   name - (character) benchmark name
##
## Result: (logical)
##
remove_benchmark <- function (bm, name = NULL) {
    if (is.null(name))  name <- bm$name

    file <- paste0(name, ".RData")

    if (!file.exists(file)) {
        return (FALSE)
    }

    file.remove(file = file)
}


## Distribution.
##
## Arguments:
##   MEAN    - (double[])
##   SD      - (double[])
##   RATE    - (double[])
##   MEANLOG - (double[])
##   SDLOG   - (double[])
##
## Result: (list)
##   benchmark
##
dist_benchmark <- function (bm, MEAN = NA, SD = NA, RATE = NA, MEANLOG = NA,
                            SDLOG = NA) {
    ## table DIST:  distribution parameters
    T_DIST <- expand_frame(ID = NA, MEAN = MEAN, SD = SD, RATE = RATE,
                           MEANLOG = MEANLOG, SDLOG = SDLOG)
    T_DIST[, "ID"] <- frame_id(T_DIST, bm$DIST)
    bm$DIST <- rbind(bm$DIST, T_DIST)

    ## done
    bm
}



## Define datasets.
##
## Arguments:
##   bm        - (list) benchmark
##   NOBS      - (integer[]) number of observations
##   NVAR      - (integer[]) number of variables (regressors)
##   NTRUE     - (integer[]) number of regressors in "true" model
##   DIST_ID   - (integer[]) distribution ID
##   INTERCEPT - (logical[])
##
## Rval: (list)
##   benchmark
##
data_benchmark <- function (bm, NOBS, NVAR, NTRUE, DIST_ID,
                            INTERCEPT = TRUE) {
    ## table DATA: artificial data
    T_DATA <- expand_frame(ID = NA, NOBS = NOBS, NVAR = NVAR, NTRUE = NTRUE,
                           DIST_ID = DIST_ID, INTERCEPT = INTERCEPT)
    T_DATA[, "ID"] <- frame_id(T_DATA, bm$DATA)
    bm$DATA <- rbind(bm$DATA, T_DATA)

    ## table SIMU:  simulations
    T_SIMU <- expand_frame(ID = NA, DATA_ID = T_DATA[, "ID"],
                           REP = seq_len(bm$nrep), SEED = NA, RSS_FUL = NA,
                           AIC_FUL = NA, BIC_FUL = NA, AIC_INF = NA,
                           BIC_INF = NA, RSS_TRU = NA, AIC_TRU = NA,
                           BIC_TRU = NA, WHICH = NA)
    T_SIMU[, "ID"] <- frame_id(T_SIMU, bm$SIMU)
    bm$SIMU <- rbind(bm$SIMU, T_SIMU)

    ## seeding
    set.seed(bm$seed)
    bm$SIMU[, "SEED"] <- sample.int(.Machine$integer.max, size = NROW(bm$SIMU))

    ## simulations
    for (simu_id in T_SIMU[, "ID"]) {
        simu <- simu_benchmark(bm, simu_id)

        simu_ix <- with(bm$SIMU, which(ID == simu_id))
        bm$SIMU[simu_ix, "RSS_FUL"] <- simu$rss_full
        bm$SIMU[simu_ix, "AIC_FUL"] <- simu$aic_full
        bm$SIMU[simu_ix, "BIC_FUL"] <- simu$bic_full
        bm$SIMU[simu_ix, "AIC_INF"] <- simu$aic_inf
        bm$SIMU[simu_ix, "BIC_INF"] <- simu$bic_inf
        bm$SIMU[simu_ix, "RSS_TRU"] <- simu$rss_true
        bm$SIMU[simu_ix, "AIC_TRU"] <- simu$aic_true
        bm$SIMU[simu_ix, "BIC_TRU"] <- simu$bic_true
        bm$SIMU[simu_ix, "WHICH"] <- paste0(as.integer(simu$true), collapse = "")
    }

    ## done
    bm
}


## Define 'lmSubsets' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   DATA_ID   - (integer[]) relevant datasets
##   NMIN      - (integer) 'lmSubsets' arg
##   NMAX      - (integer) 'lmSubsets' arg
##   TOLERANCE - (numeric) 'lmSubsets' arg
##   NBEST     - (integer) 'lmSubsets' arg
##
## Result: (list)
##   benchmark
##
lmSubsets_benchmark <- function (bm, DATA_ID, NMIN = NA, NMAX = NA,
                                 TOLERANCE = NA, NBEST = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, NMIN = NMIN, NMAX = NMAX,
                        TOLERANCE = TOLERANCE, NBEST = NBEST)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "LM_SUBSETS",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table LM_SUBSETS
    T_LM_SUBSETS <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                               NMIN = TMP[, "NMIN"], NMAX = TMP[, "NMAX"],
                               TOLERANCE = TMP[, "TOLERANCE"],
                               NBEST = TMP[, "NBEST"])
    T_LM_SUBSETS[, "ID"] <- frame_id(T_LM_SUBSETS, bm$LM_SUBSETS)
    bm$LM_SUBSETS <- rbind(bm$LM_SUBSETS, T_LM_SUBSETS)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    for (run_id in T_RUN[, "ID"]) {
        case_id <- with(T_RUN, CASE_ID[ID == run_id])
        data_id <- with(T_CASE, DATA_ID[ID == case_id])
        nvar <- with(bm$DATA, NVAR[ID == data_id])

        T_VALUE <- expand_frame(ID = NA, RUN_ID = run_id, RANK = seq_len(nvar),
                                RSS = NA, AIC = NA, BIC = NA, WHICH = NA,
                                HIT = NA)
        T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
        bm$VALUE <- rbind(bm$VALUE, T_VALUE)
    }

    ## done
    bm
}


## Define 'lmSelect' cases and runs.
##
## Arguments:
##   bm        - (list) benchmark
##   DATA_ID   - (integer[]) relevant datasets
##   IC        - (character) 'lmSelect' arg ('penalty')
##   TOLERANCE - (numeric) 'lmSelect' arg
##   NBEST     - (integer) 'lmSelect' arg
##
## Result: (list)
##   benchmark
##
lmSelect_benchmark <- function (bm, DATA_ID, IC = NA, TOLERANCE = NA,
                                NBEST = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, IC = IC, TOLERANCE = TOLERANCE,
                        NBEST = NBEST)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "LM_SELECT",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table LM_SELECT
    T_LM_SELECT <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                              IC = TMP[, "IC"], TOLERANCE = TMP[, "TOLERANCE"],
                              NBEST = TMP[, "NBEST"])
    T_LM_SELECT[, "ID"] <- frame_id(T_LM_SELECT, bm$LM_SELECT)
    bm$LM_SELECT <- rbind(bm$LM_SELECT, T_LM_SELECT)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    T_VALUE <- expand_frame(ID = NA, RUN_ID = T_RUN[, "ID"], RANK = NA,
                            RSS = NA, AIC = NA, BIC = NA, WHICH = NA, HIT = NA)
    T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
    bm$VALUE <- rbind(bm$VALUE, T_VALUE)

    ## done
    bm
}


## Define 'leaps' cases and runs.
##
## Arguments:
##   bm      - (list) benchmark
##   DATA_ID - (integer[]) relevant datasets
##   NMAX    - (integer) 'leaps' arg ('nvmax')
##   NBEST   - (integer) 'leaps' arg
##   PREORD  - (logical)  preordering
##
## Result: (list)
##   benchmark
##
leaps_benchmark <- function (bm, DATA_ID, NMAX = NA, NBEST = NA,
                             PREORD = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, NMAX = NMAX,
                        NBEST = NBEST, PREORD = PREORD)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "LEAPS",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table LEAPS
    T_LEAPS <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          NMAX = TMP[, "NMAX"], NBEST = TMP[, "NBEST"],
                          PREORD = TMP[, "PREORD"])
    T_LEAPS[, "ID"] <- frame_id(T_LEAPS, bm$LEAPS)
    bm$LEAPS <- rbind(bm$LEAPS, T_LEAPS)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    for (run_id in T_RUN[, "ID"]) {
        case_id <- with(T_RUN, CASE_ID[ID == run_id])
        data_id <- with(T_CASE, DATA_ID[ID == case_id])
        nvar <- with(bm$DATA, NVAR[ID == data_id])

        T_VALUE <- expand_frame(ID = NA, RUN_ID = run_id, RANK = seq_len(nvar),
                                RSS = NA, AIC = NA, BIC = NA, WHICH = NA,
                                HIT = NA)
        T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
        bm$VALUE <- rbind(bm$VALUE, T_VALUE)
    }

    ## done
    bm
}


## Define 'bestglm' cases and runs.
##
## Arguments:
##   bm      - (list) benchmark
##   DATA_ID - (integer[]) relevant datasets
##   IC      - (character) 'bestglm' arg ('IC')
##   NMAX    - (integer) 'bestglm' arg ('nvmax')
##   NBEST   - (integer) 'bestglm' arg ('TopModels')
##   PREORD  - (logical)  preordering
##
## Result: (list)
##   benchmark
##
bestglm_benchmark <- function (bm, DATA_ID, IC = NA, NMAX = NA, NBEST = NA,
                               PREORD = NA) {
    if (missing(DATA_ID)) {
        stop ("missing argument 'DATA_ID'")
    }

    TMP <- expand_frame(DATA_ID = DATA_ID, IC = IC, NMAX = NMAX,
                        NBEST = NBEST, PREORD = PREORD)

    ## table CASE
    T_CASE <- data_frame(ID = NA, WHAT = "BESTGLM",
                         DATA_ID = TMP[, "DATA_ID"])
    T_CASE[, "ID"] <- frame_id(T_CASE, bm$CASE)
    bm$CASE <- rbind(bm$CASE, T_CASE)

    ## table BESTGLM
    T_BESTGLM <- data_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                            IC = TMP[, "IC"], NMAX = TMP[, "NMAX"],
                            NBEST = TMP[, "NBEST"], PREORD = TMP[, "PREORD"])
    T_BESTGLM[, "ID"] <- frame_id(T_BESTGLM, bm$BESTGLM)
    bm$BESTGLM <- rbind(bm$BESTGLM, T_BESTGLM)

    ## table RUN
    T_RUN <- expand_frame(ID = NA, CASE_ID = T_CASE[, "ID"],
                          REP = seq_len(bm$nrep), EXECUTED = FALSE,
                          INTERRUPTED = NA, USER = NA, SYSTEM = NA,
                          ELAPSED = NA)
    T_RUN[, "ID"] <- frame_id(T_RUN, bm$RUN)
    bm$RUN <- rbind(bm$RUN, T_RUN)

    ## table VALUE
    T_VALUE <- expand_frame(ID = NA, RUN_ID = T_RUN[, "ID"], RANK = NA,
                            RSS = NA, AIC = NA, BIC = NA, WHICH = NA, HIT = NA)
    T_VALUE[, "ID"] <- frame_id(T_VALUE, bm$VALUE)
    bm$VALUE <- rbind(bm$VALUE, T_VALUE)

    ## done
    bm
}


## Run benchmark.
##
## Arguments:
##   bm    - (list) benchmark
##   what  - (character[]) restrict execution
##   cases - (integer[]) restrict execution
##   reps  - (integer[]) restrict execution
##   runs  - (integer[]) restrict execution
##   force - (logical) force execution
##   save  - (logical) save benchmark
##
## Result: (list)
##   benchmark
##
run_benchmark <- function (bm, what, cases, reps, runs, force = FALSE,
                           save = TRUE) {
    message ("Running benchmark '", bm$name, "'...",
             format(Sys.time(), "  <%H:%M>"))

    if (missing(what)) {
        what <- c("LM_SUBSETS", "LM_SELECT", "LEAPS", "BESTGLM")
    }

    if (missing(cases)) {
        cases <- with(bm$CASE, ID[WHAT %in% what])
    }

    if (missing(reps)) {
        reps <- seq_len(bm$nrep)
    }

    if (missing(runs)) {
        runs <- with(bm$RUN, ID[(CASE_ID %in% cases) & (REP %in% reps)])
    }

    run_ixs <- with(bm$RUN, which(ID %in% runs))
    nrun <- NROW(run_ixs)

    for (i in seq_along(run_ixs)) {
        run_ix <- run_ixs[i]

        run_id <- bm$RUN[run_ix, "ID"]
        case_id <- bm$RUN[run_ix, "CASE_ID"]
        rep <- bm$RUN[run_ix, "REP"]

        message ("(", i, "/", nrun, ")  Run ", run_id,
                 ": Case ", case_id, " [", rep, "]...",
                 format(Sys.time(), "  <%H:%M>"))

        if (!force && bm$RUN[run_ix, "EXECUTED"]) {
            message ("  ... skipping.")

            next
        }

        what <- with(bm$CASE, WHAT[ID == case_id])
        if (what == "LM_SUBSETS") {
            bm <- exec_lmSubsets_benchmark(bm, case_id, rep)
        } else if (what == "LM_SELECT") {
            bm <- exec_lmSelect_benchmark(bm, case_id, rep)
        } else if (what == "LEAPS") {
            bm <- exec_leaps_benchmark(bm, case_id, rep)
        } else if (what == "BESTGLM") {
            bm <- exec_bestglm_benchmark(bm, case_id, rep)
        } else {
            warning ("unknown strategy: ", what)
        }

        if (save)  save_benchmark(bm)

        message ("  ... done.  (",
                 format(bm$RUN[run_ix, "ELAPSED"], digits = 1, nsmall = 2),
                 "s)")
    }

    message ("... done.  (benchmark '", bm$name, "')",
             format(Sys.time(), "  <%H:%M>"))

    bm
}


## Execute 'lmSubsets' case.
##
## Arguments:
##   bm      - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Result: (list)
##   benchmark
##
exec_lmSubsets_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    lmSubsets_ix <- with(bm$LM_SUBSETS, which(CASE_ID == case_id))
    lmSubsets_id <- bm$LM_SUBSETS[lmSubsets_ix, "ID"]
    nmin <- bm$LM_SUBSETS[lmSubsets_ix, "NMIN"]
    nmax <- bm$LM_SUBSETS[lmSubsets_ix, "NMAX"]
    tolerance <- bm$LM_SUBSETS[lmSubsets_ix, "TOLERANCE"]
    nbest <- bm$LM_SUBSETS[lmSubsets_ix, "NBEST"]

    hook <- function (x, y, intercept) {
        cl <- call("lmSubsets")
        if (!is.na(nmin))  cl$nmin <- nmin + intercept
        if (!is.na(nmax))  cl$nmax <- nmax + intercept
        if (!is.na(tolerance))  cl$tolerance <- tolerance
        if (!is.na(nbest))  cl$nbest <- nbest

        cl$formula <- quote(x)
        cl$y <- quote(y)
        cl$intercept <- intercept

        time <- system.time(value <- eval(cl))

        list(time = summary(time), value = value)
    }

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, hook)

    run_ix <- with(bm$RUN, which((CASE_ID == case_id) & (REP == rep)))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ixs <- with(bm$VALUE, which(RUN_ID == run_id))

    for (val_ix in val_ixs) {
        rank <- bm$VALUE[val_ix, "RANK"]
        size <- rank + simu$intercept

        if (with(simu$value$submodel, is.na(RSS[SIZE == size])))  next

        which <- variable.names(simu$value, size = size, drop = TRUE)
        if (simu$intercept)  which <- which[-1]

        bm$VALUE[val_ix, "RSS"] <- deviance(simu$value, size = size)
        bm$VALUE[val_ix, "AIC"] <- AIC(simu$value, size = size)
        bm$VALUE[val_ix, "BIC"] <- BIC(simu$value, size = size)
        bm$VALUE[val_ix, "WHICH"] <- paste(as.integer(which), collapse = "")
        bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)
    }

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$time["user"]
    bm$RUN[run_ix, "SYSTEM"] <- simu$time["system"]
    bm$RUN[run_ix, "ELAPSED"] <- simu$time["elapsed"]

    bm
}


## Execute 'lmSelect' case.
##
## Arguments:
##   bm    - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Return: (list)
##   benchmark
##
exec_lmSelect_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    lmSelect_ix <- with(bm$LM_SELECT, which(CASE_ID == case_id))
    lmSelect_id <- bm$LM_SELECT[lmSelect_ix, "ID"]
    ic <- bm$LM_SELECT[lmSelect_ix, "IC"]
    tolerance <- bm$LM_SELECT[lmSelect_ix, "TOLERANCE"]
    nbest <- bm$LM_SELECT[lmSelect_ix, "NBEST"]

    hook <- function (x, y, intercept) {
        cl <- call("lmSelect")
        if (!is.na(ic))  cl$penalty <- ic
        if (!is.na(tolerance))  cl$tolerance <- tolerance
        if (!is.na(nbest))  cl$nbest <- nbest

        cl$formula <- quote(x)
        cl$y <- quote(y)
        cl$intercept <- intercept

        time <- system.time(value <- eval(cl))

        list(time = summary(time), value = value)
    }

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, hook)

    which <- variable.names(simu$value, drop = TRUE)
    if (simu$intercept)  which <- which[-1]

    run_ix <- with(bm$RUN, (CASE_ID == case_id) & (REP == rep))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ix <- with(bm$VALUE, which(RUN_ID == run_id))

    bm$VALUE[val_ix, "RANK"] <- sum(which)
    bm$VALUE[val_ix, "RSS"] <- deviance(simu$value)
    bm$VALUE[val_ix, "AIC"] <- AIC(simu$value)
    bm$VALUE[val_ix, "BIC"] <- BIC(simu$value)
    bm$VALUE[val_ix, "WHICH"] <- paste0(as.integer(which), collapse = "")
    bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$time["user"]
    bm$RUN[run_ix, "SYSTEM"] <- simu$time["system"]
    bm$RUN[run_ix, "ELAPSED"] <- simu$time["elapsed"]

    bm
}


## Execute 'leaps' case.
##
## Arguments:
##   bm      - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Result: (list)
##   benchmark
##
exec_leaps_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    leaps_ix <- with(bm$LEAPS, which(CASE_ID == case_id))
    leaps_id <- bm$LEAPS[leaps_ix, "ID"]
    nmax <- bm$LEAPS[leaps_ix, "NMAX"]
    nbest <- bm$LEAPS[leaps_ix, "NBEST"]
    preord <- bm$LEAPS[leaps_ix, "PREORD"]

    hook <- function (x, y, intercept) {
        if (preord) {
            x <- preorder_data(x, y, intercept)
        }

        cl <- call("regsubsets", nbest = 1, nvmax = NULL)
        if (!is.na(nbest))  cl$nbest <- nbest
        if (!is.na(nmax))  cl$nvmax <- nmax

        cl$x <- quote(x)
        cl$y <- quote(y)
        cl$intercept <- intercept

        time <- system.time(value <- eval(cl))

        list(time = summary(time), value = value)
    }

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, hook, fork = TRUE)

    sum <- summary(simu$value)

    run_ix <- with(bm$RUN, which((CASE_ID == case_id) & (REP == rep)))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ixs <- with(bm$VALUE, which(RUN_ID == run_id))

    for (val_ix in val_ixs) {
        rank <- bm$VALUE[val_ix, "RANK"]
        size <- rank

        if (is.na(sum$rss[size]))  next

        which <- sum$which[size, ]
        if (simu$intercept)  which <- which[-1]

        bm$VALUE[val_ix, "RSS"] <- sum$rss[size]
        bm$VALUE[val_ix, "AIC"] <- NA
        bm$VALUE[val_ix, "BIC"] <- sum$bic[size]
        bm$VALUE[val_ix, "WHICH"] <- paste(as.integer(which), collapse = "")
        bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)
    }

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$time["user"]
    bm$RUN[run_ix, "SYSTEM"] <- simu$time["system"]
    bm$RUN[run_ix, "ELAPSED"] <- simu$time["elapsed"]

    bm
}


## Execute 'bestglm' case.
##
## Arguments:
##   bm      - (list) benchmark
##   case_id - (integer) case
##   rep     - (integer) repetition
##
## Result: (list)
##   benchmark
##
exec_bestglm_benchmark <- function (bm, case_id, rep = 1) {
    if (missing(case_id)) {
        stop ("missing argument 'case_id'")
    }

    case_ix <- with(bm$CASE, which(ID == case_id))
    data_id <- bm$CASE[case_ix, "DATA_ID"]

    bestglm_ix <- with(bm$BESTGLM, which(CASE_ID == case_id))
    bestglm_id <- bm$BESTGLM[bestglm_ix, "ID"]
    ic <- bm$BESTGLM[bestglm_ix, "IC"]
    nmax <- bm$BESTGLM[bestglm_ix, "NMAX"]
    nbest <- bm$BESTGLM[bestglm_ix, "NBEST"]
    preord <- bm$BESTGLM[bestglm_ix, "PREORD"]

    hook <- function (x, y, intercept) {
        if (preord) {
            x <- preorder_data(x, y, intercept)
        }

        cl <- call("bestglm", TopModels = 1, nvmax = NULL)
        if (!is.na(ic))  cl$IC <- ic
        if (!is.na(nbest))  cl$TopModels <- nbest
        if (!is.na(nmax))  cl$nvmax <- nmax

        cl$Xy <- quote(as.data.frame(cbind(x, y)))
        cl$intercept <- intercept

        time <- system.time(value <- eval(cl))

        list(time = summary(time), value = value)
    }

    simu_id <- with(bm$SIMU, ID[(DATA_ID == data_id) & (REP == rep)])
    simu <- simu_benchmark(bm, simu_id, hook, fork = TRUE)

    which <- simu$value$Subsets[simu$value$ModelReport$Bestk + 1, ]
    which <- as.logical(which)
    which <- head(which, -2)
    if (simu$intercept)  which <- which[-1]

    run_ix <- with(bm$RUN, which((CASE_ID == case_id) & (REP == rep)))
    run_id <- bm$RUN[run_ix, "ID"]
    val_ix <- with(bm$VALUE, which(RUN_ID == run_id))

    bm$VALUE[val_ix, "RANK"] <- sum(which)
    bm$VALUE[val_ix, "RSS"] <- deviance(simu$value$BestModel)
    bm$VALUE[val_ix, "AIC"] <- AIC(simu$value$BestModel)
    bm$VALUE[val_ix, "BIC"] <- BIC(simu$value$BestModel)
    bm$VALUE[val_ix, "WHICH"] <- paste0(as.integer(which), collapse = "")
    bm$VALUE[val_ix, "HIT"] <- all(which == simu$true)

    bm$RUN[run_ix, "EXECUTED"] <- TRUE
    bm$RUN[run_ix, "INTERRUPTED"] <- simu$interrupted
    bm$RUN[run_ix, "USER"] <- simu$time["user"]
    bm$RUN[run_ix, "SYSTEM"] <- simu$time["system"]
    bm$RUN[run_ix, "ELAPSED"] <- simu$time["elapsed"]

    bm
}


## Simulate and run.
##
## Arguments:
##   bm      - (list) benchmark
##   simu_id - (integer) simulation ID
##   hook    - (function)
##   x, y    - (logical)
##   fork    - (logical)
##
## Result: (list)
##
simu_benchmark <- function (bm, simu_id, hook, x = FALSE, y = FALSE,
                            fork = FALSE) {
    ret_x <- x;  x <- NULL
    ret_y <- y;  y <- NULL

    rdist <- function (n, MEAN = 0, SD = 1, RATE = 1, MEANLOG = 0,
                       SDLOG = 1) {
        if (!missing(MEAN) || !missing(SD)) {
            return (rnorm(n, mean = MEAN, sd = SD))
        }

        if (!missing(RATE)) {
            return (rexp(n, rate = RATE))
        }

        if (!missing(SDLOG) || !missing(MEANLOG)) {
            return (rlnorm(n, meanlog = MEANLOG, sdlog = SDLOG))
        }

        stop ("undefined distribution")
    }

    simu_ix <- with(bm$SIMU, which(ID == simu_id))

    ## seeding
    set.seed(bm$SIMU[simu_ix, "SEED"])

    ## args
    data_id <- bm$SIMU[simu_ix, "DATA_ID"]
    data_ix <- with(bm$DATA, which(ID == data_id))

    nobs <- bm$DATA[data_ix, "NOBS"]
    nvar <- bm$DATA[data_ix, "NVAR"]
    ntrue <- bm$DATA[data_ix, "NTRUE"]
    intercept <- bm$DATA[data_ix, "INTERCEPT"]
    dist_id <- bm$DATA[data_ix, "DIST_ID"]


    ## coefficients
    coefs <- rep_len(0, length.out = nvar)
    coefs[sample.int(nvar, size = ntrue)] <- 1

    ## error
    dist_ix <- with(bm$DIST, which(ID == dist_id))
    dist <- as.list(bm$DIST[dist_ix, ])
    dist <- dist[(names(dist) != "ID") & !is.na(dist)]
    e <- do.call(rdist, c(n = nobs, dist))

    ## independent variables
    x <- rnorm(nobs * nvar)
    dim(x) <- c(nobs, nvar)

    ## dependent variable
    y <- cbind(intercept, x) %*% c(1, coefs) + e

    ## answer
    ans <- list()

    ans$nobs <- nobs
    ans$nvar <- nvar
    ans$intercept <- intercept
    ans$true <- as.logical(coefs)

    ## call
    if (!missing(hook)) {
        ret <- NULL

        if (fork) {
            job <- mcparallel({
                hook(x, y, intercept)
            }, silent = FALSE)

            if (is.null(bm$timeout)) {
                ret <- mccollect(job, wait = TRUE)[[1]]
            } else {
                ret <- mccollect(job, wait = FALSE, timeout = bm$timeout)[[1]]
            }

            pskill(pid = job$pid, signal = SIGKILL)
            pskill(pid = -1 * job$pid, signal = SIGKILL)

            mccollect(job, wait = FALSE)
        } else {
            tryCatch({
                setTimeLimit(elapsed = bm$timeout)

                ret <- hook(x, y, intercept)
            }, finally = {
                setTimeLimit(elapsed = NULL)
            })
        }

        if (is.null(ret)) {
            ans$interrupted <- TRUE
            ans$value <- NA
            ans$time <- NA
        } else {
            ans$interrupted <- FALSE
            ans$value <- ret$value
            ans$time <- ret$time
        }
    } else {
        ## stats
        x_full <- cbind(intercept, x)
        r_full <- qr.resid(qr(x_full), y)
        rss_full <- sum(r_full^2)
        ll_full <- lmSubsets:::stats_log_lik(nobs, NULL, rss_full)

        ans$rss_full <- rss_full
        ans$aic_full <- lmSubsets:::stats_aic(ll_full, 2, nvar + intercept + 1)
        ans$bic_full <- lmSubsets:::stats_bic(ll_full, nobs, nvar + intercept + 1)

        ans$aic_inf <- lmSubsets:::stats_aic(ll_full, 2, intercept + 1)
        ans$bic_inf <- lmSubsets:::stats_bic(ll_full, nobs, intercept + 1)

        x_true <- cbind(intercept, x[, ans$true])
        r_true <- qr.resid(qr(x_true), y)
        rss_true <- sum(r_true^2)
        ll_true <- lmSubsets:::stats_log_lik(nobs, NULL, rss_true)

        ans$rss_true <- rss_true
        ans$aic_true <- lmSubsets:::stats_aic(ll_true, 2, ntrue + intercept + 1)
        ans$bic_true <- lmSubsets:::stats_bic(ll_true, nobs, ntrue + intercept + 1)
    }

    if (ret_x)  ans$x <- x
    if (ret_y)  ans$y <- y

    ## done
    ans
}


## Summarize benchmark.
##
## Arguments:
##   bm - (list) benchmark
##
## Result: (list)
##   benchmark
##
summary_benchmark <- function (bm) {
    T_SUM <- with(subset(bm$RUN, EXECUTED & !INTERRUPTED),
                  data.frame(CASE_ID, TM = ELAPSED))
    T_SUM <- split(T_SUM, T_SUM[, "CASE_ID"])
    T_SUM <- lapply(T_SUM, function (grp) with(grp, {
        cbind(CASE_ID = CASE_ID[1], TM_MIN = min(TM),
              TM_MAX = max(TM), TM_AVG = mean(TM))
    }))
    T_SUM <- do.call(rbind, T_SUM)

    bm$SUMMARY <- as.data.frame(T_SUM)

    bm
}



## extract models (quick and dirty)



get_simu_benchmark <- function (bm, data_id, rep, which) {
    simu_ix <- with(bm$SIMU, base::which((DATA_ID == data_id) & (REP == rep)))
    simu_id <- bm$SIMU[simu_ix, "ID"]

    simu <- simu_benchmark(bm, simu_id, x = TRUE, y = TRUE)

    if (missing(which)) {
        ## return full model
        which <- TRUE  # all
    } else if (is.logical(which) && which) {
        ## return "true" model
        which <- simu$true
    } else if (is.character(which)) {
        ## return submodel
        which <- unlist(strsplit(which, split = ""))
        which <- as.logical(as.integer(which))
    }

    intercept <- with(bm$DATA, INTERCEPT[ID == data_id])
    if (intercept) {
        ans <- lm(simu$y ~ simu$x[, which] + 1)
    } else {
        ans <- lm(simu$y ~ simu$x[, which] + 0)
    }

    ans <- list(simu_id = simu_id, lm = ans)
    ans$rss <- deviance(ans$lm)
    ans$bic <- BIC(ans$lm)

    ans
}


get_full_benchmark <- function (bm, data_id, rep) {
    get_simu_benchmark(bm, data_id, rep)
}


get_true_benchmark <- function (bm, data_id, rep) {
    get_simu_benchmark(bm, data_id, rep, TRUE)
}


get_which_benchmark <- function (bm, case_id, rep, rank) {
    run_id <- with(bm$RUN, ID[(CASE_ID == case_id) & (REP == rep)])

    if (missing(rank)) {
        rank <- with(bm$VALUE, RANK[RUN_ID == run_id])
    }

    if (length(rank) > 1) {
        stop ("undetermined rank")
    }

    data_id <- with(bm$CASE, DATA_ID[ID == case_id])
    val_ix <- with(bm$VALUE, which((RUN_ID == run_id) & (RANK == rank)))

    ans <- get_simu_benchmark(bm, data_id, rep, bm$VALUE[val_ix, "WHICH"])

    ans$rep <- rep
    ans$rank <- rank
    ans$which <- bm$VALUE[val_ix, "WHICH"]
    ans$what <- with(bm$CASE, WHAT[ID == case_id])
    ans$run_id <- run_id
    ans$data_id <- data_id
    ans$val_id <- bm$VALUE[val_ix, "ID"]

    ans
}
