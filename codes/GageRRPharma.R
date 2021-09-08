############################################################################################################################################
### This code was written to generate the content of this publication: 																                                	 ###
### https://towardsdatascience.com/are-my-bio-pharmaceutical-assay-performances-reliable-only-probability-of-success-counts-9f85f27cb208 ###
### Author: Thomas de Marchin                                                    											                            			 ###
############################################################################################################################################

rm(list=ls())

library(openxlsx)
library(MCMCpack)
library(tidyverse)
library(VCA)
library(broom.mixed)
library(brms)
library(rstan)
library(lme4)
options(mc.cores = parallel::detectCores())

# Specs
specs <- c(90,110)

#####################
### GAGE R&R DATA ###
#####################

### Simulate data ###
batches <- data.frame(Batch=c("A", "B", "C"), ShiftBatch=c(85, 100, 110)) 
laboratories <- data.frame(Laboratory=c("Lab A", "Lab B", "Lab C", "Lab D", "Lab E"), 
                           ShiftLab=c(2.2129925, 1.7, 1.2,  2.6645,  4.15))
day <- data.frame(Day=c("1", "2"))
analyst <- data.frame(Analyst=c("1", "2"))
device <- data.frame(Device=c("1", "2"))

dataGage <- laboratories %>% 
  tidyr::crossing(day) %>% 
  tidyr::crossing(analyst) %>% 
  tidyr::crossing(device) %>% 
  tidyr::crossing(batches)

nRep <- 3
value <- c()
for(i in 1:nrow(dataGage)){
  valueTemp <- dataGage$ShiftBatch[i] +
    dataGage$ShiftLab[i] +
    rnorm(1, mean=0, sd=0.6) + #Day
    rnorm(1, mean=0, sd=1.5) + #Analyst
    rnorm(1, mean=0, sd=0.3) + #Device
    rnorm(n=nRep,mean=0,sd=0.5) #Repeatability
  value <- c(value, valueTemp)
}

dataGage <- dataGage[rep(1:nrow(dataGage), each=nRep),] %>%
  mutate(Value=value)

write.xlsx(dataGage, "data/dataGage.xlsx", overwrite = T)
# dataGage <- read.xlsx("data/dataGage.xlsx")

dataGage <- dataGage %>% mutate(Laboratory=factor(Laboratory),
                                Day=factor(Day),
                                Analyst=factor(Analyst),
                                Device=factor(Device),
                                Batch=factor(Batch))
glimpse(dataGage)

# Plot "Variability charts as a function of laboratory, day, analyst, device and batch"
png("results/varPlot.png", height=480, width=800, res=100)
varPlot(Value~(Laboratory+Day)/Analyst/Device, 
        Data=as.data.frame(dataGage),
        Mean=NULL,
        Points=list(pch=list(var="Batch", pch=c(21, 21, 21)), 
                    bg =list(var="Batch", bg=c("blue", "red", "green"))),
        Join=NULL)
dev.off()

### Model ###
dataGage <- dataGage %>% mutate(Day=paste(Day, Laboratory, sep="_"),
                                Analyst=paste(Analyst, Laboratory, sep="_"),
                                Device=paste(Device, Laboratory, sep="_"))

# First try a frequentist approach with lmer to adjust the model and then do it in bayesian (slower).
modelGageLmer <- lmer(formula=Value ~ Batch + (1|Laboratory) + (1|Day) + (1|Analyst) + (1|Device), data=dataGage)
tidy(modelGageLmer) %>%
  filter(effect=="ran_pars") %>%
  mutate(`Pct of total`=estimate^2*100/(sum(estimate^2)))

# Bayesian model with brms
modelGage <- brm(formula=Value ~ Batch + (1|Laboratory) + (1|Day) + (1|Analyst) + (1|Device), 
                 data=dataGage, chains=3, iter=50000)
saveRDS(modelGage, "results/modelGage.rds")
# modelGage <- readRDS("results/modelGage.rds")

summaryGage <- tidy(modelGage, conf.method = "HPDinterval", robust=TRUE)

summaryGage <- summaryGage %>% 
  filter(effect=="ran_pars") %>% 
  mutate(`Pct of total`=estimate^2*100/(sum(estimate^2))) %>%
  select(-effect, -term, -std.error, -component ) %>%
  rename(`Random effect`=`group`, `SD estimate`=`estimate`, `Lower 95% CI`=`conf.low`, `Upper 95% CI`=`conf.high`)
write.xlsx(summaryGage, "results/summaryGage.xlsx", overwrite = T)

chainsGage <- as.data.frame(posterior_samples(modelGage))

### laboratory precision ###
#Multiple laboratory
sigmaForTol <- sqrt(chainsGage$sd_Analyst__Intercept^2 +
                      chainsGage$sd_Day__Intercept^2 +
                      chainsGage$sd_Device__Intercept^2 +
                      chainsGage$sd_Laboratory__Intercept^2 +
                      chainsGage$sigma^2)

simulations <- rnorm(n=length(sigmaForTol), mean=0, sd=sigmaForTol) 
multipleLabPrecision <- quantile(simulations, probs=c(0.025, 0.975))

#Single laboratory
sigmaForTol <- sqrt(chainsGage$sd_Analyst__Intercept^2 +
                      chainsGage$sd_Day__Intercept^2 +
                      chainsGage$sd_Device__Intercept^2 +
                      chainsGage$sigma^2)

simulations <- rnorm(n=length(sigmaForTol), mean=0, sd=sigmaForTol) 
singleLabPrecision <- quantile(simulations, probs=c(0.025, 0.975))

labPrecisions <- tibble(name=c("Multiple lab precision (95% coverage):", "Single lab precision (95% coverage):"), 
                        Precision=paste0("X±", round(c(multipleLabPrecision[2], singleLabPrecision[2]),2)))
write.xlsx(labPrecisions, "results/labPrecisions.xlsx", overwrite = T)

### Simulate the PoS to be within spec for several batches true value ###
predictedMS <- list()
means <- seq(specs[1], specs[2], by=0.05)

predictedMS <- lapply(means,FUN=function(mean){
  
  pred <- mean + 
    rnorm(n=nrow(chainsGage), 0, sd = chainsGage$sd_Laboratory__Intercept) + 
    rnorm(n=nrow(chainsGage), 0, sd = chainsGage$sd_Day__Intercept) +
    rnorm(n=nrow(chainsGage), 0, sd = chainsGage$sd_Analyst__Intercept) +
    rnorm(n=nrow(chainsGage), 0, sd = chainsGage$sd_Device__Intercept) +
    rnorm(n=nrow(chainsGage), 0, sd = chainsGage$sigma)
  
  PoS <- mean(specs[1] <= pred & pred <= specs[2])
  
  return(data.frame(Mean=mean,PoS=PoS))
})

predictedMS <- bind_rows(predictedMS)

# Plot "Probability of success to be within specification as a function of the batch true value"
pMS <- ggplot(data=predictedMS, aes(x=Mean, y=PoS)) + theme_bw() +
  labs(x="'True' Batch Value", y="Probability of Success to be within specifications") +
  geom_vline(xintercept = specs[1],colour="red",
             linetype="dashed") +
  geom_vline(xintercept = specs[2],colour="red",
             linetype="dashed") +
  geom_path() + ylim(c(0.5,1)) 
ggsave("results/pMS.png", pMS, height=10, width=14, unit="cm")

### Generate a prior for CPV data ###
randomEffectsWithoutLaboratory <- chainsGage$sd_Analyst__Intercept^2 + chainsGage$sd_Day__Intercept^2 + chainsGage$sd_Device__Intercept^2 + chainsGage$sigma^2
shape=mean(randomEffectsWithoutLaboratory)**2/var(randomEffectsWithoutLaboratory)+2
scale=mean(randomEffectsWithoutLaboratory)**3/var(randomEffectsWithoutLaboratory)+mean(randomEffectsWithoutLaboratory)

# check if the inverse gamma fits well the posterior
png("results/prior.png", height=480, width=800)
#plot the original chain
plot(density(randomEffectsWithoutLaboratory))
#plot the computed inverse gamma
lines(density(rinvgamma(10000, shape=shape, scale=scale)), col="red")
dev.off()

################
### CPV DATA ###
################

### Simulate data ###
dataCPV <- data.frame(Batch=LETTERS, Value=rnorm(26, mean=104, sd=2.5)) 
write.xlsx(dataCPV, "data/dataCPV.xlsx", overwrite = T)
# dataCPV <- read.xlsx("data/dataCPV.xlsx")

controlLimitsCPV <- c(mean(dataCPV$Value)+3*sd(dataCPV$Value), mean(dataCPV$Value)-3*sd(dataCPV$Value))

# plot "Measured value of batches produced"
pCPV <- ggplot(aes(x=Batch, y=Value), data=dataCPV) + 
  geom_point() + theme_bw() +
  geom_hline(yintercept=c(specs[1], specs[2]), col="red") +
  ggtitle("CPV data") +
  geom_hline(yintercept=c(controlLimitsCPV[1], controlLimitsCPV[2]), col="green", linetype="dashed") 
ggsave("results/pCPV.png", pCPV, height=10, width=14, unit="cm")

### Model ###

# Bayesian model with stan. We use stan and not brms as we can fine tune the model more easily with stan.
modelCPV <- '
data {
  int<lower=0> N;
  //vector[N] x;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> sigma2e;
  real<lower=0> sigma2b; //batch
}
model {
sigma2e ~ inv_gamma(10.64, 35.71);
y ~ normal(alpha, sqrt(sigma2b+sigma2e));
}
'
B=list()
formula = Value ~ 1
B$x = model.matrix(formula,dataCPV)
B$y = dataCPV$Value
B$N = nrow(dataCPV)

iter <- 40000
chains <- 3
thin <- 1
warmup <- 20000

modelCPV <- stan(model_code = modelCPV, data = B, thin = thin, warmup = warmup, iter = warmup+(iter/chains))
saveRDS(modelCPV, "results/modelCPV.rds")
# modelCPV <- readRDS("results/modelCPV.rds")

summaryCPV <- tidy(modelCPV, conf.int = T, conf.level = 0.95) %>% select(-std.error) %>%
  mutate(estimate = ifelse(grepl("sigma2", term), sqrt(estimate), estimate)) %>%
  mutate(conf.low = ifelse(grepl("sigma2", term), sqrt(conf.low), conf.low)) %>%
  mutate(conf.high = ifelse(grepl("sigma2", term), sqrt(conf.high), conf.high)) %>%
  mutate(term=recode(term, alpha="Intercept", sigma2b="Random Batch (SD)", sigma2e="Residual (SD)")) %>%
  rename(`Effect`=`term`, `SD estimate`=`estimate`, `Lower 95% CI`=`conf.low`, `Upper 95% CI`=`conf.high`)
write.xlsx(summaryCPV, "results/summaryCPV.xlsx", overwrite = T)

traceplot(modelCPV)
chainsCPV <- extract(modelCPV)

### predictions ### 

# only the batch to batch included
predictionsCPVOnlyBatch <- rnorm(n=nrow(chainsCPV$alpha), mean=chainsCPV$alpha, sd=sqrt(chainsCPV$sigma2b))
qtCPVOnlyBatch <- as.data.frame(t(quantile(predictionsCPVOnlyBatch, probs = c(0.025,0.975))))

# batch to batch + measurement variability from Gage R&R
# samples as chains from CPV and Gage R&R do not have the same length
samples <- sample(1:length(randomEffectsWithoutLaboratory), length(chainsCPV$sigma2b), replace=TRUE)
predictionsCPVBatchAndResidual <- rnorm(n=nrow(chainsCPV$alpha), mean=chainsCPV$alpha, sd=sqrt(randomEffectsWithoutLaboratory[samples]+chainsCPV$sigma2b))
qtCPVBatchAndResidual <- as.data.frame(t(quantile(predictionsCPVBatchAndResidual, probs = c(0.025,0.975))))

# Plot "CPV prediction"
pCPVPred <- ggplot(aes(x=Batch, y=Value), data=dataCPV) + 
  geom_point() + theme_bw() +
  geom_hline(yintercept=c(specs[1], specs[2]), col="red") +
  ggtitle("CPV data") + 
  geom_ribbon(aes(x=as.numeric(factor(LETTERS)), ymin=qtCPVBatchAndResidual$`2.5%`, ymax=qtCPVBatchAndResidual$`97.5%`), alpha=0.3, fill="blue") +
  geom_ribbon(aes(x=as.numeric(factor(LETTERS)), ymin=qtCPVOnlyBatch$`2.5%`, ymax=qtCPVOnlyBatch$`97.5%`), alpha=0.3, fill="red")
ggsave("results/pCPVPred.png", pCPVPred, height=10, width=14, unit="cm")

# make predictions for each lab from Gage R&R
chainsGageLabs <- chainsGage[,grepl("r_Laboratory", names(chainsGage))]

simulationsPerLab <- lapply(1:ncol(chainsGageLabs), function(i){
  simulations <- chainsGageLabs[samples,i] + #intercept from each lab from Gage R&R
    chainsCPV$alpha + #intercept from CPV
    rnorm(n=nrow(chainsCPV$alpha), mean=0, sd=sqrt(randomEffectsWithoutLaboratory[samples])) + #run-to-run from Gage R&R
    rnorm(n=nrow(chainsCPV$alpha), mean=0, sd=sqrt(chainsCPV$sigma2b)) #batch-to-batch from CPV
  
  return(data.frame(simulations))
})

names(simulationsPerLab) <- names(chainsGageLabs)

simulationsPerLab <- bind_rows(simulationsPerLab, .id="Lab") %>% 
  mutate(Lab=str_remove_all(Lab, "r_Laboratory\\[")) %>% 
  mutate(Lab=str_remove_all(Lab, ",Intercept\\]"))

# Plot "Posterior predictive distribution of future measurement as a function of the measuring laboratory"
pMeasurementsInLab <- ggplot(simulationsPerLab, aes(simulations, colour = Lab)) +
  geom_density() +
  geom_vline(xintercept = specs[1],colour="red",
             linetype="dashed") +
  geom_vline(xintercept = specs[2],colour="red",
             linetype="dashed") + 
  theme_bw() + 
  labs(x="Range")
ggsave("results/pMeasurementsInLab.png", pMeasurementsInLab, height=10, width=14, unit="cm")

# compute PoS to be within spec per lab
computePoS <- function(x, specs) {
  mean(specs[1] <= x & x <= specs[2])
}

PoSLabs <- simulationsPerLab %>% 
  group_by(Lab) %>% summarise(computePoS(simulations, specs))
write.xlsx(PoSLabs, "results/PoSLabs.xlsx", overwrite = T)