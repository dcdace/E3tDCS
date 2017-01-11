library(plyr) # for function 'count'
library(BayesFactor)

# import the data
tDCSdata = read.csv(file = "T:/_CogNeurPhD/03-Experiment_tDCS/tDCSraw.csv", header = T)

# BF. support for the effect (H1), based on Jeffreys
support = c("anecdotal", "substantial", "strong", "very strong", "decisive")
cuts = c(1, 3, 10, 30, 100, Inf)

# ---- CALCULATE NECCESSARY VARIABLES
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

# which scores (columns) to process
scoresPre =
  c("ITPreDiff",
    "ITPreDiff",
    "ETPreDiff",
    "ETPreDiff",
    "ErrPreDiff",
    "ErrPreDiff")
scoresPost =
  c("ITPostDiff",
    "ITRetDiff",
    "ETPostDiff",
    "ETRetDiff",
    "ErrPostDiff",
    "ErrRetDiff")

scN = length(scoresPost) # how many scores

# group sample sizes
groupCount = count(tDCSdata, 'stimulation')
shamN = subset(groupCount, stimulation == 0)$freq
activeN = subset(groupCount, stimulation == 1)$freq

## OVERALL TRAINING EFFECT

trEffect = mapply(function(post, pre, stim)
{
  f = paste(post, "~", pre)
  summary(lm(as.formula(f), data = subset(tDCSdata, stimulation ==
                                            stim)))
},
scoresPost, scoresPre, c(rep(0, scN), rep(1, scN)))

trEffect_tvals = sapply(1:dim(trEffect)[2], function(x)
  round(coef(trEffect[, x])[, "t value"][1], 2))

trEffect_pvals = sapply(1:dim(trEffect)[2], function(x)
  round(coef(trEffect[, x])[, "Pr(>|t|)"][1], 3))

trEffect_pvalsTXT = lapply(trEffect_pvals, function(x)
  ifelse(x < 0.001, "< 0.001", paste("=", round(x, 3))))

# Bayes Factor
# calculate Bayes Factor from t-test results

# get the effect size that the design could detect with 80% power and p < 0.05
dSham = power.t.test(n = shamN, NULL, 1, .05, .8, "one.sample", "two.sided")$d
dActive = power.t.test(n = activeN, NULL, 1, .05, .8, "one.sample", "two.sided")$d

bf10 = mapply(
  function(tt, n, d)
    ttest.tstat(tt, n, rscale = d, simple = TRUE),
  trEffect_tvals,
  c(rep(shamN, scN), rep(activeN, scN)),
  c(rep(dSham, scN), rep(dActive, scN))
)
bf01 = 1 / bf10

supportTXT = mapply(function(x, y)
{
  ifelse(
    x > 1
    ,
    paste(as.character(cut(
      x, breaks = cuts, labels = support
    )), "evidence for the effect")
    ,
    paste(as.character(cut(
      y, breaks = cuts, labels = support
    )), "evidence against the effect")
  )
}
, bf10, bf01)

# format results as txt strings
resTrainingTXT = mapply(
  function(df, t, p, supptxt, bf10, bf01)
    sprintf(
      "t(%d) = %.2f, p %s, %s (BF10/BF01 = %.2f/%.2f).",
      df,
      t,
      p,
      supptxt,
      bf10,
      bf01
    ),
  c(rep(shamN, scN), rep(activeN, scN)),
  trEffect_tvals,
  trEffect_pvalsTXT,
  supportTXT,
  bf10,
  bf01
)
# put sham and active results in separate columns
resTrainingTXT2 = matrix(
  resTrainingTXT,
  nrow = scN,
  ncol = 2,
  dimnames = list(scoresPost, c("sham", "active"))
)

## ---- STIMULATION EFFeCT
# Comparing Sham and Active Post~Pre Intercepts
# See: http://rcompanion.org/rcompanion/e_04.htm

# degrees of freedom = total(N) - N(groups) - N(slopes)
df = dim(tDCSdata)[1] - 2 - 2

# t and p
lmPost =  mapply(function(post, pre)
{
  f = paste(post, "~", pre, "* stimulation")
  summary(lm(as.formula(f), data = tDCSdata))
},
scoresPost, scoresPre)

pvals = sapply(1:dim(lmPost)[2], function(x)
  coef(lmPost[, x])[, "Pr(>|t|)"][3])
tvals = sapply(1:dim(lmPost)[2], function(x)
  coef(lmPost[, x])[, "t value"][3])

pvalsTXT = lapply(pvals, function(x)
  ifelse(x < 0.001, "< 0.001", paste("=", round(x, 3))))

# Bayes Factors
BFlmPost = as.data.frame(mapply(function(post, pre)
{
  f = paste(post, "~", pre, "+ stimulation")
  regressionBF(as.formula(f), data = tDCSdata)
},
scoresPost, scoresPre))

BF10vals = BFlmPost[2, c(TRUE, rep(FALSE, 3))]
BF01vals = 1 / BF10vals

supportTXT = mapply(function(x, y)
{
  ifelse(
    x > 1
    ,
    paste(as.character(cut(
      x, breaks = cuts, labels = support
    )), "evidence for the effect")
    ,
    paste(as.character(cut(
      y, breaks = cuts, labels = support
    )), "evidence against the effect")
  )
}
, BF10vals, BF01vals)

# format results as txt strings
resTXT = mapply(
  function(df, t, p, supptxt, bf10, bf01)
    sprintf(
      "t(%d) = %.2f, p %s, %s (BF10/BF01 = %.2f/%.2f).",
      df,
      t,
      p,
      supptxt,
      bf10,
      bf01
    ),
  df,
  as.numeric(tvals),
  pvalsTXT,
  supportTXT,
  as.numeric(BF10vals),
  as.numeric(BF01vals)
)
# put in matrix
resTXT2 = matrix(
  resTXT,
  nrow = scN,
  ncol = 1,
  dimnames = list(scoresPost, "stimulation effect (group difference)")
)

# put training results and stimulation results in the same matrix
results = cbind(resTrainingTXT2, resTXT2)

print(results)
write.table(results, file = "T:/_CogNeurPhD/03-Experiment_tDCS/resultsIntercepts.txt")
