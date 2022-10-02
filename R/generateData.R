##' Generate simulated data
##'
##' Generate multistate data simulating time-varying treatment assignment, baseline and time-varying confounding, and death.
##' @param n Number of individuals
##' @import data.table
##' @export
##' @examples
##' n <- 100
##' generateData(n)
##' dt <- readRDS(paste0("data/sim_data_n", n, ".rds")
##'
##' # Transitions
##' survival::survcheck(survival::Surv(from, to, to.state.fac)~1, dt, id = id)
##' # State occupancy probabilities
##' sfit <- survival::survfit(survival::Surv(from, to, to.state.fac)~X, dt, id = id)
##' # Solid lines are for individuals with X = 0, dashed lines are for those with X = 1
##' plot(sfit, col = rep(c(1, 3, 2), each = 2), lty = 1:2)
##' text(rep(5, 3), c(0.1, 0.6, 0.3), c("conf.t", "treat", "death"), col = c(1, 3, 2))
##'
generateData <- function(n) {

  if (!requireNamespace("miscFunctions", quietly = TRUE)) {
    stop(
      "Package \"miscFunctions\" must be installed to use this function.
      \nInstall it by devtools::install_github('palryalen/miscFunctions')",
      call. = FALSE
    )
  }

  set.seed(1)

  # Generate multistate data with 5 states ----------------------------------

  # Transition rates
  a13 <- 0.22; a14 <- 0.7; a15 <- 0.08
  a23 <- 0.45; a24 <- 0.45; a25 <- 0.1
  a34 <- 0.17; a35 <- 0.33
  a45 <- 0.02

  # Time grid
  end_time <- 10
  xx <- seq(0, end_time, length.out = 20)

  # It is possible to specify time-varying transition rates. Here, they are
  # constant in time
  spl13 <- stats::smooth.spline(xx, seq(a13, a13, length.out = length(xx)))
  spl14 <- stats::smooth.spline(xx, seq(a14, a14, length.out = length(xx)))
  spl15 <- stats::smooth.spline(xx, seq(a15, a15, length.out = length(xx)))
  spl23 <- stats::smooth.spline(xx, seq(a23, a23, length.out = length(xx)))
  spl24 <- stats::smooth.spline(xx, seq(a24, a24, length.out = length(xx)))
  spl25 <- stats::smooth.spline(xx, seq(a25, a25, length.out = length(xx)))
  spl34 <- stats::smooth.spline(xx, seq(a34, a34, length.out = length(xx)))
  spl35 <- stats::smooth.spline(xx, seq(a35, a35, length.out = length(xx)))
  spl45 <- stats::smooth.spline(xx, seq(a45, a45, length.out = length(xx)))

  # Transition matrix
  transition_matrix <- vector("list", length = 5)
  transition_matrix[[1]]$smoothspline[[1]] <- list(NULL)
  transition_matrix[[1]]$smoothspline[[2]] <- list(NULL)
  transition_matrix[[1]]$smoothspline[[3]] <- list(spl13)
  transition_matrix[[1]]$smoothspline[[4]] <- list(spl14)
  transition_matrix[[1]]$smoothspline[[5]] <- list(spl15)
  transition_matrix[[1]]$neighbours <- c(3, 4, 5)
  transition_matrix[[2]]$smoothspline[[1]] <- list(NULL)
  transition_matrix[[2]]$smoothspline[[2]] <- list(NULL)
  transition_matrix[[2]]$smoothspline[[3]] <- list(spl23)
  transition_matrix[[2]]$smoothspline[[4]] <- list(spl24)
  transition_matrix[[2]]$smoothspline[[5]] <- list(spl25)
  transition_matrix[[2]]$neighbours <- c(3, 4, 5)
  transition_matrix[[3]]$smoothspline[[1]] <- list(NULL)
  transition_matrix[[3]]$smoothspline[[2]] <- list(NULL)
  transition_matrix[[3]]$smoothspline[[3]] <- list(NULL)
  transition_matrix[[3]]$smoothspline[[4]] <- list(spl34)
  transition_matrix[[3]]$smoothspline[[5]] <- list(spl35)
  transition_matrix[[3]]$neighbours <- c(4, 5)
  transition_matrix[[4]]$smoothspline[[1]] <- list(NULL)
  transition_matrix[[4]]$smoothspline[[2]] <- list(NULL)
  transition_matrix[[4]]$smoothspline[[3]] <- list(NULL)
  transition_matrix[[4]]$smoothspline[[4]] <- list(NULL)
  transition_matrix[[4]]$smoothspline[[5]] <- list(spl45)
  transition_matrix[[4]]$neighbours <- c(5)
  transition_matrix[[5]]$smoothspline[[1]] <- list(NULL)
  transition_matrix[[5]]$smoothspline[[2]] <- list(NULL)
  transition_matrix[[5]]$smoothspline[[3]] <- list(NULL)
  transition_matrix[[5]]$smoothspline[[4]] <- list(NULL)
  transition_matrix[[5]]$smoothspline[[5]] <- list(NULL)
  transition_matrix[[5]]$neighbours <- NULL

  # Starting states: 1 and 2
  dfr <- miscFunctions::generateFrame(n, end_time, transition_matrix, c(1, 2))

  # Terminating state: 5
  dfr <- dfr[!(dfr$from.state == 5),]

  # Create labels for each state --------------------------------------------

  # state --> label
  # 1 --> diag
  # 2 --> diag
  # 3 --> conf.t
  # 4 --> treat
  # 5 --> dead

  dt <- data.table::as.data.table(dfr)

  # X: baseline confounder
  dt[, X := ifelse(from.state[[1]] == 1, 0, 1), by = id]

  # L: time-varying confounder
  dt[, L := ifelse(to.state == 3, 1, 0)]
  dt[, L := cumsum(L), by = id]
  # L(t) = I(T_L < t) is left-continous, while Lcurrent(t) = I(T_L <= t) is
  # right-continuous
  dt[, Lcurrent := L]
  dt[to.state == 3, L := 0]

  # A: time-varying treatment
  dt[, A := ifelse(to.state == 4, 1, 0)]
  dt[, A := cumsum(A), by = id]
  dt[to.state == 4, A := A - 1]

  # D: death
  dt[, D := 0]
  dt[to.state == 5, D := 1]

  dt[from.state %in% c(1, 2), from.state.str := "diag"]
  dt[from.state == 3, from.state.str := "conf.t"]
  dt[from.state == 4, from.state.str := "treat"]

  dt[to.state == 3, to.state.str := "conf.t"]
  dt[to.state == 4, to.state.str := "treat"]
  dt[to.state == 5, to.state.str := "dead"]

  # Censored for death
  dt[to == 10, to.state := 0]
  dt[to == 10, to.state.str := "cens"]

  # Create a factor variable of states which are transitioned to. This is only
  # relevant for checking state occupancy probabilities
  dt[, to.state.fac := factor(to.state, levels = c(0, 3, 4, 5),
                              labels = c("cens", "conf.t", "treat", "dead"))]

  # Delete and rename variable
  dt[, c("from.state", "to.state") := NULL]
  data.table::setnames(dt, "from.state.str", "from.state")
  data.table::setnames(dt, "to.state.str", "to.state")
  data.table::setcolorder(dt, c("id", "from", "to", "from.state", "to.state", "to.state.fac"))

  ## Create an at-risk for treatment indicator variable
  dt[, atRiskTreat := 1]
  dt[, tmp := cumsum(from.state == "treat"), by = id]
  dt[, tmp := ifelse(tmp > 1, 1, tmp)]
  dt[, atRiskTreat := atRiskTreat - tmp]
  dt[, tmp := NULL]

  # Save simulated data
  filename <- "data/sim_data_n"
  saveRDS(dt, file = paste0(filename, n, ".rds"))

}
