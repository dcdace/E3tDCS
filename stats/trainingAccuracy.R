library(tidyr)
library(ez)
library(BayesFactor)

# import the data
tDCSdata = read.csv(file = "T:/_CogNeurPhD/03-Experiment_tDCS/tDCSraw.csv", header = T)

# BF. support for the effect (H1), based on Jeffreys
support = c("anecdotal", "substantial", "strong", "very strong", "decisive")
cuts = c(1, 3, 10, 30, 100, Inf)

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

print(resultsAcc)
write.table(resultsAcc, file = "T:/_CogNeurPhD/03-Experiment_tDCS/resultsTrainingAccuracy.txt")