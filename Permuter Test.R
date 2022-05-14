library(permuter)
head(pneumovac)

# randomization-based inference to determine there was a difference in bp_episodes between trt vs ctrl
# test statistic is the estimated log IRR using a poisson GLM
# will calculate the p-value based on 1000 permutations
test <- permtest_glm(bpepisodes~spnvac, trtname="spnvac", runit="randunit", family=poisson, data=pneumovac, nperm=1000, ncores=6)
print(c(test$coef, test$pval))
print(c(exp(test$coef), test$pval)) 
plot(test)

# alternative way to run this permutation test is to specify your poisson model (m1) and pass m1 into permtest function

# get randomization-based 95% CI 
ci <- permci_glm(bpepisodes~spnvac, trtname="spnvac", runit="randunit", family=poisson, data=pneumovac, nperm=1000, ncores=2, seed=445,
                 level=0.95, initmethod="perm")

ci$ci


# Time to event data + matched pairs
library(survival)
test <- permtest_survreg(Surv(left, right, type = 'interval2') ~ treat, 
                         trtname = "treat", runit = "group", strat = "pair.id", 
                         data = bcpp, nperm = 10, ncores = 3, seed = 446)
ci <- permci(m1, trtname = "treat", runit  = "group", strat = "pair.id", data = bcpp, nperm = 10, ncores = 2, seed = 447)
print(c(test$coef, test$pval))
ci$ci

ds <- gendata_crt(family = binomial, nclus = c(10, 10), size = c(30, 50),
                  theta = log(1.5), mu = qlogis(0.25), sigma = 0.2)

