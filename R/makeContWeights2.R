#' Make continuous-time weights for a time-change intervention
#'
#' This function is based on `makeContWeights` from the `ahw`-package to implement a time-change intervention.
#'
#' @param faFit An object of type 'aalen' returned by `timereg::aalen`, describing the hazard rate of treatment in the observed scenario
#' @param dataFr A data.table on long format containing multiple observations per individual of state transitions times
#' @param eventState A string of the name of the event state of interest (treatment)
#' @param plotWeights Plot the estimated weights? (TRUE/FALSE)
#'
#' @return A data.table of estimated weights for each individual patient. As long as an individual is at risk for `eventState` (treatment), his weight process jumps whenever `eventState` (treatment) occurs in the patient population. The jumps occur at times `tstop`.
#' @export
#'
#' @examples
#' dataFr <- readRDS("data/sim_data_n100.rds")
#'
#' faFit <- timereg::aalen(survival::Surv(from, to, to.state == "treat") ~ 1 + X + L, # confounders
#'   data = dt[atRiskTreat == 1],
#'   id = dt[atRiskTreat == 1, id])
#'
#' frame <- makeContWeights2(faFit, dataFr, plotWeights = T)
#'
makeContWeights2 <- function(faFit, dataFr, eventState = "treat", plotWeights = FALSE) {

  stopifnot("atRiskTreat" %in% names(dataFr))
  stopifnot("thetaFunc" %in% names(dataFr))

  # Weight calculation ------------------------------------------------------

  # Data relevant for weight calculations: as long as an individual is at risk
  # for treatment.
  wtFrame <- dataFr[atRiskTreat == 1]
  # predTimes: The treatment times in the population
  predTimes <- wtFrame[to.state == eventState, to]
  predTimeIds <- wtFrame[to.state == eventState, id]

  # Predictions from the additive model for treatment, used to calculate the
  # integrator (dK) the integral equation for the weights R_t = R_0 + \int_0^t
  # R_{t-} dK_t
  fPred <- timereg::predict.aalen(faFit, newdata = wtFrame, n.sim = 0, se = F, resample.iid = 0,
                                  times = c(0, sort(predTimes)))

  # Weight trajectories for each individual, evaluated at each treatment time as
  # long as the individual himself is at risk for treatment
  predTable <- weightPredict2(fPred, wtFrame, predTimes, predTimeIds, eventState, plotWeights)

  # Individuals at risk -----------------------------------------------------

  # Expanded data, where each individual's follow-up time is expanded to include
  # the treatment times in the population, as long as the individual himself is
  # at risk for treatment
  Table <- refineTable2(dataFr, eventState = eventState)

  # Merge in the calculated weights -----------------------------------------

  # Merge so that Table includes the weight estimates and theta. Weights are
  # merged so that they evaluate to their value just before the event time
  # specified by 'to', while theta attains its value at the event time 'to'
  setnames(predTable, "t", "to")
  Table <- merge(Table,predTable[, .(id, to, W, theta)],by.x=c("id","from"), by.y = c("id", "to"),all.x=T)
  # Table <- merge(Table,predTable[, .(id, to, theta)],by=c("id","to"),all.x=T)
  setcolorder(Table, c("id", "from", "to"))

  # Fill in missing values of theta from the values copied from 'baseTable'.
  # These values are for event times after an individual has been treated
  Table[is.na(theta), theta := thetaFunc]
  Table[, thetaFunc := NULL]

  # Set weights to the last available value. This should apply to weights after
  # treatment time (this check is not implemented here!) This also applies to
  # weight estimates at times after max.time, which are the times after which the
  # design matrix in the treatment models (obs/hyp) could be non-invertible
  # (again, this check is also not implemented here!)
  Table[,W:=zoo::na.locf0(W),by=id]

  return(Table)

}

weightPredict2 <- function(fPred, wtFrame, predTimes, predTimeIds, eventState, plotWeights = FALSE) {

  # Predicted marginal cum. hazard based on (observed) additive model for
  # eventState, used in the calculation of the weights
  dA_f <- as.vector(apply(fPred$S0,1,function(rw)-diff(c(0,log(rw)))))

  # Prediction times (totTimes): Treatment times, in addition to time zero
  totTimes <- unique(c(0, sort(predTimes)))
  nTimes <- length(totTimes)
  mtf <- match(c(0,predTimes),totTimes) # the order of prediction times as sorted by id

  # Create a new data table which contains repeated instances of the vector of
  # prediction times; it is repeated for each row of wtFrame. Hence, if an
  # individual is represented by three rows in wtFrame, that is, if he
  # experiences two events prior to treatment/censoring, the vector of
  # prediction times is stacked three times for this individual
  predTable <- data.table(rowNum = rep(1:nrow(wtFrame), each = nTimes), # rows of wtFrame that are repeated
                          id = rep(wtFrame$id, each = nTimes), # id, from, to and theta are from wtFrame
                          from = rep(wtFrame$from, each = nTimes),
                          to = rep(wtFrame$to, each = nTimes),
                          theta = rep(wtFrame$thetaFunc, each = nTimes),
                          t = rep(totTimes, nrow(wtFrame)),  # prediction times (not from wtFrame!)
                          dA_f = dA_f)

  # In case of several risk states per individual ---------------------------

  # The number of rows in wtFrame per individual
  numRepId <- as.numeric(table(wtFrame$id))

  # Row number in wtFrame
  wtFrame$rowNum <- 1:nrow(wtFrame)

  # The first and last rowNum per individual
  predTable[, c("rowNumFirst", "rowNumLast") := .(rowNum[[1]], rowNum[[.N]]), by = id]

  # Keep rows where the prediction times fall in the interval of exposure times
  predTable[, keep := 1 * (from <= t & t < to)]
  # In addition, also keep rows of prediction times at the very beginning and at
  # the very end, which happen to fall outside the range of the exposure time
  # for an indivdual
  predTable[keep == 0 & rowNum == rowNumFirst, keep := 1 * (t <= min(to)), by = id]
  predTable[keep == 0 & rowNum == rowNumLast, keep := 1 * (max(to) <= t), by = id]

  predTable <- predTable[keep==1]
  predTable <- predTable[order(id, t)]
  predTable[, keep := NULL]

  # update last row
  predTable[, rowNumLast := rowNum[[.N]], by = id]

  # Define the at risk for treatment indicator variable ---------------------

  # Indicator variable which denotes whether an individual received treatment
  predTable[, event := 1 * (id %in% predTimeIds), by = id]
  # The treatment time for an individual (defined over all his rows)
  predTable[event == 1, eventTime := predTimes[predTimeIds == id], by = id]
  # An indicator variable denoting the row where the treatment event occurs
  predTable[, eventRow := 1 * (eventTime == t)]

  # In case the individual did not receive treatment, we define a variable with
  # the censoring time
  predTable[event == 0, censoredTime := to[rowNum == rowNumLast][[1]], by = id]
  predTable[, c("rowNum", "rowNumFirst", "rowNumLast", "from", "to") := NULL]

  # The at risk indicator variable
  predTable[event==1, atRisk := as.numeric(t <= eventTime)]
  predTable[event==0, atRisk := as.numeric(t <= censoredTime)]
  predTable[, c("event", "eventTime", "censoredTime") := NULL]

  # Calculation of weights --------------------------------------------------

  # Calculate the increments of the integrator (dK), using the integral equation
  # for the weights R_t = R_0 + \int_0^t R_{t-} dK_t
  predTable[is.na(eventRow), eventRow := 0]
  predTable[, dK := eventRow - atRisk * dA_f]

  predTable[, W := cumprod(1 + (theta - 1) * dK), by = id]

  predTable[, c("dA_f", "eventRow", "dK") := NULL]

  # Plot weights
  if(plotWeights) {
    # y-limit
    Wmax <- max(predTable$W)
    # transparent color
    blackTrans <- rgb(red = 0, green = 0, blue = 0, alpha = 0.2)

    # plot weights of the first patient
    firstId <- predTable$id[[1]]
    plot(predTable[id == firstId, .(t, W)],
         type = "s", col = blackTrans, ylim = c(0, Wmax + 0.5))

    # keep summing the weights for each treatment time, over all individuals.
    # Note: here we assume that the weights of every individual is represented
    # at each treatment time, and that predTable is sorted by id and t
    weightSum <- predTable[id == firstId, W]

    # plot the weights for all other patients, except the first
    for(i in sort(unique(predTable$id)[-1])){
      lines(predTable[id == i, .(t, W)], col = blackTrans, type = "s")

      # keep adding weights
      weightSum <- weightSum + predTable[id == i, W]
    }

    # calculate and plot the average weight at each treatment time
    weightMean <- weightSum / length(unique(predTable$id))
    lines(sort(unique(predTable$t)), weightMean, col = "blue", type = "s")
  }

  return(predTable)
}

refineTable2 <- function(dataFr, eventState){

  stopifnot("atRiskTreat" %in% names(dataFr))

  baseTable <- copy(dataFr) # deep copy

  # Refine at treatment times in the population when at risk
  # (atRiskTreat == 1 is not really needed, since every individual is at
  # risk for treatment before being treated)
  eventTimes <- baseTable[atRiskTreat == 1 & to.state == eventState, to]

  # Defining rownumber
  baseTable[, rowNumber := 1:.N]

  # numRep is the number of times a given row should be replicated
  baseTable[, numRep := 1]

  # When the subjects are at risk of being treated, they should have one row for
  # each treatment time in the population. Note that by = rowNumber is needed
  # because we are comparing each 'from' and 'to' to a vector of 'eventTimes'
  baseTable[atRiskTreat == 1,
            numRep := 1 + as.numeric(sum(from < eventTimes & eventTimes < to)),
            by = rowNumber]

  eventTimesSorted <- sort(eventTimes)

  # Expanding the data.table. Table is a copy of baseTable, except that it has
  # repeated rows
  Table <- baseTable[rep(1:nrow(baseTable), times = baseTable$numRep),
                     .(id, from, to, from.state, to.state, atRiskTreat, rowNumber, thetaFunc)]

  # Preparing the event times in the population that falls within the interval
  # where the subject is at risk. Realize that rowNumber refers to the rows of
  # 'baseTable', that is, it is not unique in 'Table', which is why we need to
  # specify from[[1]] and to[[1]], otherwise the latter are vectors of varying
  # sizes not comparable to the size of eventTimesSorted
  Table[atRiskTreat == 1,
        putEventTimes := c(eventTimesSorted[from[[1]] < eventTimesSorted & eventTimesSorted < to[[1]]], to[[1]]),
        by = rowNumber]

  # Setting the new event times from the population
  Table[!is.na(putEventTimes), to := putEventTimes]
  Table[, from := c(from[[1]], to[-.N]), by = id]

  # Setting to.state = eventState.pop as long as subjects are at risk of treatment
  Table[atRiskTreat == 1,
        to.state := c(rep(paste0(eventState,".pop"), .N - 1), baseTable[rowNumber, to.state]),
        by = rowNumber]
  # Matching from.state to the newly set to.state
  Table[atRiskTreat == 1,
        from.state := c(from.state[[1]], to.state[-.N]),
        by = rowNumber]

  # Removing redundant columns
  Table <- subset(Table,select = !(names(Table) %in% c("numRep","rowNumber","putEventTimes")))

  return(Table)
}
