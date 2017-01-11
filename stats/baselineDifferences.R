library(plyr) # for function 'count'

# import the data
tDCSdata = read.csv(file = "T:/_CogNeurPhD/03-Experiment_tDCS/tDCSraw.csv", header = T)

# ---- CALCULATE NECCESSARY VARIABLES
# Mean sensations
tDCSdata$maxMean = rowMeans(subset(tDCSdata, select = c(max1, max2, max3, max4)))
tDCSdata$affectedMean = rowMeans(subset(tDCSdata, select = c(
  affected1, affected2, affected3, affected4
)))
tDCSdata$stoppedMean = rowMeans(subset(tDCSdata, select = c(stop1, stop2, stop3, stop4)))

# Baseline perfrmance
tDCSdata$ITPre = (tDCSdata$ITPreTR + tDCSdata$ITPreUN) / 2
tDCSdata$ETPre = (tDCSdata$ETPreTR + tDCSdata$ETPreUN) / 2
tDCSdata$ErrPre = (tDCSdata$ErrPreTR + tDCSdata$ErrPreUN) / 2

# ---- BASELINE DIFFERENCES
scores = c("age",
           "ITPre",
           "ETPre",
           "ErrPre",
           "maxMean",
           "affectedMean",
           "stoppedMean")
groupCount = count(tDCSdata, 'stimulation')
shamN = subset(groupCount, stimulation == 0)$freq
activeN = subset(groupCount, stimulation == 1)$freq

# Gender proportion
genderCount = count(tDCSdata, c('stimulation', 'gender'))

meansSham = lapply(scores, function(x)
  mean(subset(tDCSdata[, x], tDCSdata$stimulation == 0)))
meansActive = lapply(scores, function(x)
  mean(subset(tDCSdata[, x], tDCSdata$stimulation == 1)))
sdSham = lapply(scores, function(x)
  sd(subset(tDCSdata[, x], tDCSdata$stimulation == 0)))
sdActive = lapply(scores, function(x)
  sd(subset(tDCSdata[, x], tDCSdata$stimulation == 1)))

genderDiffpval = wilcox.test(
  gender ~ stimulation,
  data = tDCSdata,
  alternative = "two.sided",
  correct = TRUE,
  exact = 1
)$p.value
genderDiffpvalTxt = ifelse(genderDiffpval < 0.001, "< 0.001", round(genderDiffpval, 3))

baselineDiffpval = lapply(scores, function(x)
  wilcox.test(
    tDCSdata[, x] ~ stimulation,
    data = tDCSdata,
    alternative = "two.sided",
    correct = TRUE,
    exact = 1
  )$p.value)
baselineDiffpvalTxt = lapply(baselineDiffpval, function(x)
  ifelse(x < 0.001, "< 0.001", round(x, 3)))

resultsGender = matrix(
  c(
    paste(
      subset(genderCount, stimulation == 0 &
               gender == 2)$freq,
      ":",
      subset(genderCount, stimulation == 0 & gender == 1)$freq,
      sep = ""
    ),
    # sham
    paste(
      subset(genderCount, stimulation == 1 &
               gender == 2)$freq,
      ":",
      subset(genderCount, stimulation == 1 & gender == 1)$freq,
      sep = ""
    ),
    # active
    genderDiffpvalTxt
  ),
  nrow = 1,
  ncol = 3,
  dimnames = list("gender M:F", c(
    paste("Sham (N =", shamN, ")"),
    paste("Active (N =", activeN, ")"),
    "Difference p-val"
  ))
)

# other differences
resultsOther = matrix(
  c(
    paste(
      lapply(meansSham, function(x)
        round(x, 2)),
      "+-",
      lapply(sdSham, function(x)
        round(x, 2))
    ),
    paste(
      lapply(meansActive, function(x)
        round(x, 2)),
      "+-",
      lapply(sdActive, function(x)
        round(x, 2))
    ),
    baselineDiffpvalTxt
  ),
  nrow = 7,
  ncol = 3,
  dimnames = list(scores, c(
    paste("Sham (N =", shamN, ")"),
    paste("Active (N =", activeN, ")"),
    "Difference p-val"
  ))
)
resultsBaseline = rbind(resultsGender, resultsOther)

print(resultsBaseline)
write.table(resultsBaseline, file = "T:/_CogNeurPhD/03-Experiment_tDCS/resultsBaseline.txt")
