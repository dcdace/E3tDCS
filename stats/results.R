library(plyr)
library(tidyr)
library(ez)
library(BayesFactor)

# import the data
tDCSdata = read.csv(file = "T:/_CogNeurPhD/03-Experiment_tDCS/tDCSraw.csv", header = T)

# BF. support for the effect (H1), based on Jeffreys
support = c("anecdotal", "substantial", "strong", "very strong", "decisive")
cuts = c(1, 3, 10, 30, 100, Inf)

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

# UN TR difference
tDCSdata$ITPreDiff = tDCSdata$ITPreUN / tDCSdata$ITPreTR - 1
tDCSdata$ITPostDiff = tDCSdata$ITPostUN / tDCSdata$ITPostTR - 1
tDCSdata$ITRetDiff = tDCSdata$ITRetUN / tDCSdata$ITRetTR - 1

tDCSdata$ETPreDiff = tDCSdata$ETPreUN / tDCSdata$ETPreTR - 1
tDCSdata$ETPostDiff = tDCSdata$ETPostUN / tDCSdata$ETPostTR - 1
tDCSdata$ETRetDiff = tDCSdata$ETRetUN / tDCSdata$ETRetTR - 1

tDCSdata$ErrPreDiff = tDCSdata$ErrPreUN - tDCSdata$ErrPreTR
tDCSdata$ErrPostDiff = tDCSdata$ErrPostUN - tDCSdata$ErrPostTR
tDCSdata$ErrRetDiff = tDCSdata$ErrRetUN - tDCSdata$ErrRetTR

# ---- BASELINE DIFFERENCES
scores = c("age","ITPre","ETPre","ErrPre","maxMean","affectedMean","stoppedMean")
groupCount = count(tDCSdata, 'stimulation')
shamN = subset(groupCount, stimulation==0)$freq
activeN = subset(groupCount, stimulation==1)$freq

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
genderDiffpvalTxt = ifelse(genderDiffpval < 0.001, "< 0.001", round(genderDiffpval,3))

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
    paste(subset(genderCount, stimulation==0 & gender==2)$freq, ":", subset(genderCount, stimulation==0 & gender==1)$freq, sep = ""),
    # sham
    paste(subset(genderCount, stimulation==1 & gender==2)$freq, ":", subset(genderCount, stimulation==1 & gender==1)$freq, sep = ""),
    # active
    genderDiffpvalTxt
  ),
  nrow = 1,
  ncol = 3,
  dimnames = list("gender M:F", c(paste("Sham (N =", shamN, ")"), paste("Active (N =", activeN, ")"), "Difference p-val"))
)

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
  dimnames = list(scores, c(paste("Sham (N =", shamN, ")"), paste("Active (N =", activeN, ")"), "Difference p-val"))
)
resultsBaseline = rbind(resultsGender, resultsOther)

#---- OBSERvATION AccURACY ANOVA (DAY,STIMULATION)
# Convert it to long format
data_longAcc = gather(tDCSdata, trDay, Acc, D1Acc:D4Acc)
data_longAcc$trDay = factor(data_longAcc$trDay)
data_longAcc$sID = factor(data_longAcc$sID)
data_longAcc$stimulation = factor(data_longAcc$stimulation)

# get ANOVAs and Bayes Factors
accANOVA = ezANOVA(
  data = data_longAcc,
  dv = Acc,
  wid = sID,
  within = .(trDay),
  between = stimulation,
  type = 2
)

accBFANOVA = anovaBF(
  Acc ~ stimulation * trDay + sID,
  data = data_longAcc,
  whichRandom = "sID"
)

txtSigDay = ifelse(accANOVA$ANOVA$p[2] < 0.05, "Significant", "No significant")
txtpvalDay = ifelse(accANOVA$ANOVA$p[2] < 0.001, "p < 0.001", paste("p = ", round(accANOVA$ANOVA$p[2], 3), sep = ""))

txtSigGroup = ifelse(accANOVA$ANOVA$p[1] < 0.05, "Significant", "No significant")
txtpvalGroup = ifelse(accANOVA$ANOVA$p[1] < 0.001, "p < 0.001", paste("p = ", round(accANOVA$ANOVA$p[1], 3), sep = ""))

bfAnova10 = exp(1) ^ accBFANOVA@bayesFactor$bf
bfAnova01 = 1 / exp(1) ^ accBFANOVA@bayesFactor$bf

txtBFDay = ifelse(
  bfAnova10[2] > 1
  , paste(as.character(cut(bfAnova10[2], breaks = cuts, labels = support)), "evidence for")
  , paste(as.character(cut(bfAnova01[2], breaks = cuts, labels = support)), "evidence against")
)
txtBFGroup = ifelse(
  bfAnova10[1] > 1
  , paste(as.character(cut(bfAnova10[1], breaks = cuts, labels = support)), "evidence for")
  , paste(as.character(cut(bfAnova01[1], breaks = cuts, labels = support)), "evidence against")
)

resTxt = c(
  "Accuracy: %s main effect of Day, F(%d,%d) = %.2f, %s, ng = %.2f; %s the effect (BF10/BF01 = %.2f/%.2f). \n"
  ,
  "%s main effect of Group, F(%d,%d) = %.2f, %s, ng = %.2f; %s the effect (BF10/BF01 = %.2f/%.2f). \n"
)
resultsAcc = sprintf(
  paste(resTxt, collapse = "")
  ,
  txtSigDay,
  accANOVA$ANOVA$DFn[2],
  accANOVA$ANOVA$DFd[2],
  round(accANOVA$ANOVA$F[2], 2),
  txtpvalDay,
  round(accANOVA$ANOVA$ges[2], 2),
  txtBFDay,
  round(bfAnova10[2], 2),
  round(bfAnova01[2], 2)
  ,
  txtSigGroup,
  accANOVA$ANOVA$DFn[1],
  accANOVA$ANOVA$DFd[1],
  round(accANOVA$ANOVA$F[1], 2),
  txtpvalGroup,
  round(accANOVA$ANOVA$ges[1], 2),
  txtBFGroup,
  round(bfAnova10[1], 2),
  round(bfAnova01[1], 2)
)

#---- EXECUTiON ANOVA (DAY,STIMULATION) ------
scores = c("IT", "ET", "Err")
# Convert it to long format
data_longIT = gather(tDCSdata, test, Diff, ITPreDiff:ITRetDiff)
data_longIT$test = factor(data_longIT$test)
data_longIT$sID = factor(data_longIT$sID)
data_longIT$stimulation = factor(data_longIT$stimulation)

data_longET = gather(tDCSdata, test, Diff, ETPreDiff:ETRetDiff)
data_longET$test = factor(data_longET$test)
data_longET$sID = factor(data_longET$sID)
data_longET$stimulation = factor(data_longET$stimulation)

data_longErr = gather(tDCSdata, test, Diff, ErrPreDiff:ErrRetDiff)
data_longErr$test = factor(data_longErr$test)
data_longErr$sID = factor(data_longErr$sID)
data_longErr$stimulation = factor(data_longErr$stimulation)

# put all converted data frames into list
dataList = list(data_longIT, data_longET, data_longErr)

# get ANOVAs and Bayes Factors
resANOVA = lapply(dataList,  function(x)
  ezANOVA(
    data = x,
    dv = Diff,
    wid = sID,
    within = .(test),
    between = stimulation,
    type = 2
  ))
bfANOVA = lapply(dataList, function(x)
  anovaBF(
    Diff ~ stimulation * test + sID,
    data = x,
    whichRandom = "sID"
  ))

# format results
# for each data frame (score)

txtSigDay = lapply(resANOVA, function(x)
  ifelse(x$ANOVA$p[2] < 0.05, "Significant", "No significant"))
txtpvalDay = lapply(resANOVA, function(x)
  ifelse(x$ANOVA$p[2] < 0.001, "p < 0.001", paste("p = ", round(x$ANOVA$p[2], 3), sep = "")))

txtSigGroup = lapply(resANOVA, function(x)
  ifelse(x$ANOVA$p[1] < 0.05, "Significant", "No significant"))
txtpvalGroup = lapply(resANOVA, function(x)
  ifelse(x$ANOVA$p[1] < 0.001, "p < 0.001", paste("p = ", round(x$ANOVA$p[1], 3), sep = "")))
bfAnova10 = lapply(bfANOVA, function(x)
  exp(1) ^ x@bayesFactor$bf)
bfAnova01 = lapply(bfANOVA, function(x)
  1 / exp(1) ^ x@bayesFactor$bf)

txtBFDay = mapply(function(x, y)
{
  ifelse(
    x[2] > 1
    , paste(as.character(cut(x[2], breaks = cuts, labels = support)), "evidence for")
    , paste(as.character(cut(y[2], breaks = cuts, labels = support)), "evidence against")
  )
}
, bfAnova10, bfAnova01)

txtBFGroup = mapply(function(x, y)
{
  ifelse(
    x[1] > 1
    ,
    paste(as.character(cut(
      x[1], breaks = cuts, labels = support
    )), "evidence for")
    ,
    paste(as.character(cut(
      y[1], breaks = cuts, labels = support
    )), "evidence against")
  )
}
, bfAnova10, bfAnova01)

resTxt = c(
  "%s: %s main effect of Day, F(%d,%d) = %.2f, %s, ng = %.2f; %s the effect (BF10/BF01 = %.2f/%.2f). \n",
  "%s main effect of Group, F(%d,%d) = %.2f, %s, ng = %.2f; %s the effect (BF10/BF01 = %.2f/%.2f). \n"
)
resultsANOVA = mapply(
  function(scores,txtSigDay,resANOVA,txtpvalDay,txtBFDay,bfAnova10,bfAnova01,
           txtSigGroup,txtpvalGroup,txtBFGroup)
  {
    sprintf(
      paste(resTxt, collapse = "")
      ,
      scores,
      txtSigDay,
      resANOVA$ANOVA$DFn[2],
      resANOVA$ANOVA$DFd[2],
      round(resANOVA$ANOVA$F[2], 2),
      txtpvalDay,
      round(resANOVA$ANOVA$ges[2], 2),
      txtBFDay,
      round(bfAnova10[2], 2),
      round(bfAnova01[2], 2)
      ,
      txtSigGroup,
      resANOVA$ANOVA$DFn[1],
      resANOVA$ANOVA$DFd[1],
      round(resANOVA$ANOVA$F[1], 2),
      txtpvalGroup,
      round(resANOVA$ANOVA$ges[1], 2),
      txtBFGroup,
      round(bfAnova10[1], 2),
      round(bfAnova01[1], 2)
    )
  }
  ,
  scores,txtSigDay,resANOVA,txtpvalDay,txtBFDay,bfAnova10,bfAnova01,txtSigGroup,txtpvalGroup,txtBFGroup
)

#---- T-TEST: independent samples, two-tailed ----------------

# get the effect size that the design could detect with 80% power and p < 0.05
d = power.t.test(n = 25, NULL, 1, .05, .8, "two.sample", "two.sided")$d

# which scores (columns) to process
scores <-
  c("ITPostDiff","ITRetDiff","ETPostDiff","ETRetDiff","ErrPostDiff","ErrRetDiff")

ttest = lapply(scores, function(x)
  t.test(tDCSdata[, x] ~ stimulation, data = tDCSdata, var.eq = TRUE))
txtSig = lapply(ttest, function(x)
  ifelse(round(x$p.value, 3) < 0.05, "Significant", "No significant"))

# BAYES FACTOR
x = lapply(scores, function(x)
  subset(tDCSdata[, x], tDCSdata$stimulation == 0))
y = lapply(scores, function(x)
  subset(tDCSdata[, x], tDCSdata$stimulation == 1))

bf = mapply(function(x, y) {
  ttestBF(x, y, rscale = d)
}, x, y)

bf10 = lapply(bf, function(x)
  round(exp(1) ^ x@bayesFactor$bf, 2))
bf01 = lapply(bf, function(x)
  round(1 / exp(1) ^ x@bayesFactor$bf, 2))

txt = mapply(function(x, y)
{
  ifelse(
    x > 1
    ,
    paste(as.character(cut(
      x, breaks = cuts, labels = support
    )), "evidence for")
    ,
    paste(as.character(cut(
      y, breaks = cuts, labels = support
    )), "evidence against")
  )
}
, bf10, bf01)

# EFFECT SIZE: cohen's d'
cd = mapply(function(tt, x, y)
  tt$statistic * sqrt(1 / length(x) + 1 / length(y)), ttest, x, y)
# 95% CI
ciSham = lapply(x, function(x)
  sd(x) / sqrt(length(x)) * 1.96)
ciActive = lapply(y, function(y)
  sd(y) / sqrt(length(y)) * 1.96)

# RESULtS
resTxt = c(
  "%s: %s difference between sham (M = %.2f, 95%% CI [%.2f, %.2f])",
  "and active (M = %.2f, 95%% CI [%.2f, %.2f]) groups;",
  "t(%d) = %.3f, p = %.3f, d = %.2f;",
  "%s the effect (BF10/BF01 = %.2f/%.2f) \n"
)
resultsTTests = mapply(function(scores,txtSig,x,ciSham,y,ciActive,ttest,cd,txt,bf10,bf01)
{
  sprintf(
    paste(resTxt, collapse = " ")
    ,
    scores,
    txtSig,
    round(mean(x), 2),
    round(mean(x) - ciSham, 2),
    round(mean(x) + ciSham, 2),
    round(mean(y), 2),
    round(mean(y) - ciActive, 2),
    round(mean(y) + ciActive, 2),
    dim(tDCSdata)[1] - 2,
    round(ttest$statistic, 3),
    round(ttest$p.value, 3),
    round(abs(cd), 2),
    txt,
    bf10,
    bf01
  )
}
, scores,txtSig,x,ciSham,y,ciActive,ttest,cd,txt,bf10,bf01)

resultsEffects = c(resultsAcc, resultsANOVA, resultsTTests)

print(resultsBaseline)
print(resultsEffects)

write.table(resultsBaseline, file = "T:/_CogNeurPhD/03-Experiment_tDCS/resultsBaseline.txt")
write.table(resultsEffects, file = "T:/_CogNeurPhD/03-Experiment_tDCS/resultsEffects.txt")
