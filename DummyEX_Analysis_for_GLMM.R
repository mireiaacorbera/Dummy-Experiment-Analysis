#Dummy Exp - Script for Analysis####
#import dataset####
setwd("C:/Users/User/Desktop/Galapagos/data/Field season 2024/dummy ex/final documents/analysis")

data <- read.csv("~/Desktop/dummyEX_28.06.2024.csv", stringsAsFactors = TRUE, na.strings = "", sep = ";", dec=".")

head(data)
#View(data)

#packages used####

#
source("~/Desktop/diagnostic_fcns.R")
#source("C:/Users/User/Desktop/Galapagos/glmm_stability.R")
#source("C:/Users/User/Desktop/Galapagos/glmmTMB_stability.R")

library(glmmTMB)
library(sjPlot)
library(DHARMa)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library(interactions)
#install.packages("kyotil")
library(kyotil)
library(ggeffects)
#install.packages("ggnewscale")
library(ggnewscale)
#install.packages("effects")
library(effects)
#install.packages("data.table")
library(data.table)
library(MuMIn)
library(performance)
#install.packages("censReg")
library(censReg)
#install.packages("lmerTest")
library(lmerTest)
library(dplyr)
library(MASS)
library(tidyr)
library(rcompanion)
library("writexl")

#data preparation####

#exclude rows with only NA 147-162
data <- data %>% slice(-(147:162))
#View(data)



# Convert the 'date' column from dd.mm.yyyy format to Date object
data$date <- as.Date(data$date, format = "%d.%m.%Y")
#reference_date <- as.POSIXct("01.01.2024", format = "%d.%m.%Y")
#data$duration.since.ref <- as.numeric(difftime(data$date, reference_date, units = "days"))

# Convert the 'time' column from hh:mm:ss format to POSIXct object
data$time <- as.POSIXct(data$time, format = "%H:%M:%S", tz = "UTC")
data$hour <- as.numeric(format(data$time, "%H"))

data$nest.ID<-as.factor(data$nest.ID)
data$dummy<-as.factor(data$dummy)
data$enemy.type<-as.factor(data$enemy.type)
data$pres.order<-as.factor(data$pres.order)
data$stage<-as.factor(data$stage)
data$Philornis<-as.numeric(data$Philornis)
data$latency.exit<-as.numeric(data$latency.exit)
data$PC<-as.factor(data$PC)
data$sgf.nb<-as.factor(data$sgf.nb)
data$ind.pres<-as.numeric(data$ind.pres)
data$parents.call.initial<-as.factor(data$parents.call.initial)
data$parents.alarm.1<-as.factor(data$parents.alarm.1)
data$parents.pres<-as.factor(data$parents.pres)

str(data)

#feed <- subset(data, stage=="Feeding")
parents <- subset(data, parents.alarm=="1")

#alarm calls (community)####

#some R function
#using diagnostics_fcns.R from Roger Mundry. We determine theoretically possible random slopes and remove incomplete cases and dummy code categorical predictors
names(data)
xx.fe.re<-fe.re.tab(fe.model="alarm.calls~dummy*stage+pres.order+visibility+number.nb+tot.ind.alarm", 
                    re="(1|nest.ID) + (1|date)",
                    data=data)

nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

#centering of dummy coded predictors (categorical)
model.data$dummy.Control=model.data$dummy.Control-mean(model.data$dummy.Control)
model.data$dummy.SGF=model.data$dummy.SGF-mean(model.data$dummy.SGF)
model.data$dummy.Owl=model.data$dummy.Owl-mean(model.data$dummy.Owl)

model.data$pres.order.2=model.data$pres.order.2-mean(model.data$pres.order.2)
model.data$pres.order.3=model.data$pres.order.3-mean(model.data$pres.order.3)
model.data$pres.order.4=model.data$pres.order.4-mean(model.data$pres.order.4)


#z-transforming of covariate predictors
model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))
model.data$z.tot.ind.alarm<-as.vector(scale(model.data$tot.ind.alarm))

#set dummy "Control"as the reference level
model.data$dummy <- relevel(model.data$dummy, ref = "Control")

max.mod.1<-glmmTMB(alarm.calls~dummy*stage+pres.order+z.visibility+z.number.nb+tot.ind.alarm+
                     (1+tot.ind.alarm|nest.ID)+
                     (1+dummy.Control+dummy.SGF+dummy.Owl+pres.order.2+pres.order.3+pres.order.4+z.visibility+z.number.nb+tot.ind.alarm|date),
                   data=model.data, ziformula=~1, family=poisson(link="log"))

#model has problems to converge due to complexity (random slopes)
#reduce random slopes
#double pipe to remove correlations in random slopes
full.mod.1<-glmmTMB(alarm.calls~dummy*stage+pres.order+z.visibility+z.number.nb+tot.ind.alarm+
                      (1|nest.ID),
                    data=model.data, family=poisson(link="log"))

#calculate lm model without interactions to check vifs (collinearity)
mod.1<-lm(alarm.calls~dummy+stage+pres.order+z.visibility+z.number.nb+tot.ind.alarm,
          data=model.data)

vif_values <- vif(mod.1)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.319*1.319
max.vif


#check model assumptions: dispersion, zero-inflation
hist(model.data$alarm.calls)

testDispersion(full.mod.1)
simulationOutput <- simulateResiduals(fittedModel = full.mod.1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

#zero inflation!!! include ziformula and exclude z.visibility to make model converge

full.mod.1<-glmmTMB(alarm.calls~dummy*stage+pres.order+z.number.nb+z.tot.ind.alarm+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=poisson(link="log"))

testDispersion(full.mod.1)
#test dispersion: test not significant
simulationOutput <- simulateResiduals(fittedModel = full.mod.1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
#qq-plot residuals okay, residual vs. predicted??
testZeroInflation(simulationOutput)
#zi test not significant

#Blup - best linear unbiased predictor: should be +- normally distributed and sholdn?t be skewed (+-)
#range of x-axis shouldn?t be >4 (for random slopes); no random slopes in my model
ranef.diagn.plot(full.mod.1)

null.mod.1<-glmmTMB(alarm.calls~pres.order+z.tot.ind.alarm+z.number.nb+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=poisson(link="log"))

as.data.frame(anova(full.mod.1, null.mod.1, test="Chisq"))
#full model explains variation better than null model (without interaction), lower AIC
#if null model doesn?t converge it has to be adapted (reducing covariates) and full model must also be adapted and diagnostics must be checked again (until both modls converge)

#model stability - problems with formula (not necessarily important to be checked but I will try to fix this)
#m.stab=glmm.model.stab(model.res=full.mod.1)-> only for lme4 models
#calculates many models reducing categories of predictor variables, to see how stable a predictor is when levels vary. The wider confidence interval of estimate the more unstable the model is
#full.stab=glmmTMB.stab(model.res=full.mod.1, data=model.data)
#round(full.stab$summary[, -1], 3)
#m.stab.plot(full.stab$summary[, -1])


#after diagnostics check out and model comparison is significant we look at the results

summary(full.mod.1)

r_sq_p = performance::r2(model = full.mod.1) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.1)
r_sq_p
r_sq_m

#manually change reference order to compare groups
model.data$pres.order = relevel(model.data$pres.order, ref = "1")
#significant effect for presentation order, wind, hour

emm <- emmeans(full.mod.1, ~ dummy*stage)
pairwise_comparisons <- contrast(emm, method = "pairwise")
pairwise_comparisons
#plot(emm, CIs = F, comparisons = TRUE)
#export emmeans summary and pairwise comparisons summaries
emm_df <- as.data.frame(emm)
pairwise_comparison_df <- as.data.frame(pairwise_comparisons)
write_xlsx(pairwise_comparison_df, path = "AC_community_pairwise_comparisons.xlsx")
write_xlsx(emm_df, path = "AC_community_emmeans_summary.xlsx")

plot(allEffects(full.mod.1))

#PLOT
ggplot(model.data, aes(x=stage, y=alarm.calls))+
  geom_boxplot(aes(fill=dummy))

# plot raw data with model estimates 
# get model estimates
plot(ggpredict(full.mod.1, terms = c("dummy", "stage")))
xxp = ggpredict(full.mod.1, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

p <- ggplot(df1, aes(x = dummy, y = alarm.calls)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of alarm calls", fill="dummy", color="dummy", title = "Number of alarm calls emitted by the community")

ggsave("AC_Community_plot.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

means_dt <- df1 %>%
  group_by(dummy, stage) %>%
  summarise(mean_alarm_calls = mean(alarm.calls, na.rm = TRUE))

p.1 <- ggplot(df1, aes(x = dummy, y = alarm.calls)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3, show.legend = FALSE)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  geom_point(data = means_df, aes(x = dummy, y = mean_alarm_calls), 
             color = "black", size = 3, shape = 18, 
             position = position_dodge(0.9)) +  # ensure points are dodged if necessary
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of alarm calls", fill="dummy", color="dummy", title = "Number of alarm calls emitted by the community")
print(p.1)
ggsave("AC_Community_simplified.png", plot = p.1, width = 12, height = 8, units = "in", dpi = 300)

#alarm calls (only parents)####
#we have to discuss if it is meaningful to exclude all cases other individuals were warning too.
#maybe warning of other individuals is part of the response and when excluding this we see only half of the truth??


#only cases where ind.alarm=0 (no other individuals were alarming)
data.ind<-subset(data, ind.alarm=="0")
str(data.ind)

xx.fe.re<-fe.re.tab(fe.model="alarm.calls~dummy*stage+pres.order+visibility+number.nb+total.time.in+parents.alarm", 
                    re="(1|nest.ID)",
                    data=data.ind)

nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)


full.mod.1<-glmmTMB(alarm.calls~dummy*stage+pres.order+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=poisson(link = "log"))

mod.1<-lm(alarm.calls~dummy+stage+pres.order,
          data=model.data)

vif_values <- vif(mod.1)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.018*1.018
max.vif

overdisp.test(full.mod.1)

testDispersion(full.mod.1)
simulationOutput <- simulateResiduals(fittedModel = full.mod.1, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.1)

null.mod.1<-glmmTMB(alarm.calls~pres.order+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=poisson(link = "log"))

as.data.frame(anova(full.mod.1, null.mod.1, test="Chisq"))
#no difference between full and null when family=nbinom2!

summary(full.mod.1)

r_sq_p = performance::r2(model = full.mod.1) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.1)
r_sq_p
r_sq_m

#manually change reference order to compare groups
model.data$pres.order = relevel(model.data$pres.order, ref = "1")
#significant effect for presentation order, wind, hour

emm <- emmeans(full.mod.1, ~ dummy*stage)
pairwise_comparisons <- contrast(emm, method = "pairwise")
pairwise_comparisons
#plot(emm, CIs = F, comparisons = TRUE)
#Export summary of emmeans and pairwise_comparisons
emm_df <- as.data.frame(emm)
pairwise_comparison_df <- as.data.frame(pairwise_comparisons)
write_xlsx(pairwise_comparison_df, path = "AC_parents_pairwise_comparisons.xlsx")
write_xlsx(emm_df, path = "AC_parents_emmeans_summary.xlsx")

plot(allEffects(full.mod.1))

#PLOT
ggplot(model.data, aes(x=stage, y=alarm.calls))+
  geom_boxplot(aes(fill=dummy))

# plot raw data with model estimates
# get model estimates
plot(ggpredict(full.mod.1, terms = c("dummy", "stage")))
xxp = ggpredict(full.mod.1, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

p <- ggplot(df1, aes(x = dummy, y = alarm.calls)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of alarm calls", fill="dummy", color="dummy", title = "Number of alarm calls emitted by parents only")
ggsave("AC_ParentsOnly.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

#plot only with mean alarm calls, boxes, median & data points

p.1 <- ggplot(df1, aes(x = dummy, y = alarm.calls)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3, show.legend = FALSE)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
 #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  geom_point(aes(x=dummy, y=predicted, fill=dummy), color = "black", size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of alarm calls", fill="dummy", color="dummy", title = "Number of alarm calls emitted by parents only")
print(p.1)
ggsave("AC_ParentsOnly_simplified.png", plot = p.1, width = 12, height = 8, units = "in", dpi = 300)

#latency alarm#####

#cases when no calls were emitted NA
#0/1 calling and compare that
#include all NAs with 300

#data$latency.alarm.adj=data$latency.alarm
#data$latency.alarm.adj[is.na(data$latency.alarm.adj)]=300

xx.fe.re<-fe.re.tab(fe.model="latency.alarm~dummy*stage+pres.order+visibility+number.nb", 
                    re="(1|nest.ID)+(1|date)",
                    data=data)

nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)

range(model.data$latency.alarm)
model.data$latency.alarm <- log(model.data$latency.alarm + 1)


model.data$latency.alarm <- model.data$latency.alarm + 1e-6 
#(add small value to latency to avoid 0 and make gaussian model run)

hist(model.data$latency.alarm)

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))

full.mod.2<-glmmTMB(latency.alarm~dummy*stage+pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    family= "gaussian"(link = "identity"), ziformula=~1, 
                    data=model.data)

mod.2<-lm(latency.alarm~dummy+stage+pres.order+z.visibility+z.number.nb,
          data=model.data)


vif_values <- vif(mod.2)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.048*1.048
max.vif

hist(model.data$latency.alarm)

testDispersion(full.mod.2)
simulationOutput <- simulateResiduals(fittedModel = full.mod.2, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.2)

null.mod.2<-glmmTMB(latency.alarm~pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    family= "gaussian"(link = "identity"),ziformula=~1,
                    data=model.data)

as.data.frame(anova(full.mod.2, null.mod.2, test="Chisq"))

summary(full.mod.2)


sub.mod.2<-glmmTMB(latency.alarm~dummy+stage+pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    family= "gaussian"(link = "identity"),ziformula = ~1,
                    data=model.data)

testDispersion(sub.mod.2)
simulationOutput <- simulateResiduals(fittedModel = sub.mod.2, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(sub.mod.2)

#plot(model.data$latency.alarm.adj, model.data$dummy)
summary(sub.mod.2)

#closest approach####

xx.fe.re<-fe.re.tab(fe.model="approach~dummy*stage+pres.order+visibility+number.nb+parents.pres", 
                    re="(1|nest.ID)+(1|date)",
                    data=data)

nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

model.data$approach

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))

full.mod.3<-glmmTMB(approach~dummy*stage+dummy+pres.order+z.visibility+number.nb+parents.pres+
                      (1|nest.ID),
                    data=model.data, family=gaussian(link="identity"))

mod.3<-lm(approach~dummy+stage+pres.order+z.visibility+number.nb+parents.pres,
          data=model.data)

vif_values <- vif(mod.3)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.051*1.051
max.vif

hist(model.data$approach)

testDispersion(full.mod.3)
simulationOutput <- simulateResiduals(fittedModel = full.mod.3, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.3)

null.mod.3<-glmmTMB(approach~pres.order+z.visibility+number.nb+parents.pres+
                      (1|nest.ID),
                    data=model.data, family=gaussian(link="identity"))

as.data.frame(anova(full.mod.3, null.mod.3, test="Chisq"))

summary(full.mod.3)

emm <- emmeans(full.mod.3, ~ dummy|stage)
pairwise_comparisons <- contrast(emm, method = "pairwise")
pairwise_comparisons

r_sq_p = performance::r2(model = full.mod.3) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.3)
r_sq_p
r_sq_m

#manually change reference order to compare groups
#model.data$pres.order = relevel(model.data$pres.order, ref = "1")


emmeans_results <- emmeans(full.mod.3, ~ dummy |stage)

# Pairwise comparisons
pairwise_comparisons <- pairs(emmeans_results)

# Print the results
print(pairwise_comparisons)

#export emmeans summary and pairwise comparisons summaries
emm_df <- as.data.frame(emmeans_results)
pairwise_comparison_df <- as.data.frame(pairwise_comparisons)
write_xlsx(pairwise_comparison_df, path = "Pairwise_comparison_ClosestApproach.xlsx")
write_xlsx(emm_df, path = "Emmeans_ClosestApproach.xlsx")


plot(allEffects(full.mod.3))

ggplot(model.data, aes(x=dummy, y=approach))+
  geom_boxplot(aes(fill=dummy))

plot(ggpredict(full.mod.3, terms = c("dummy", "stage")))
xxp = ggpredict(full.mod.3, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

p <- ggplot(df1, aes(x = dummy, y = approach)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Closest approach", fill="dummy", color="dummy", title = "Closest approach to dummy")
ggsave("Closest_approach.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)
print (p)
#closer approach to owl than to control

#Simplified plot without error bars and average indicated

p.1 <- ggplot(df1, aes(x = dummy, y = approach)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3, show.legend = FALSE)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  geom_point(aes(x=dummy, y=predicted, fill=dummy), color = "black", size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Closest approach", fill="dummy", color="dummy", title = "Closest approach to dummy")
ggsave("Closest_approach_simplified.png", plot = p.1, width = 12, height = 8, units = "in", dpi = 300)
print (p.1)

#peak frequency####

boxplot(parents$peak.freq~parents$dummy)

xx.fe.re<-fe.re.tab(fe.model="peak.freq~dummy*stage+pres.order+visibility+parents.call.initial", 
                    re="(1|nest.ID)+(1|date)",
                    data=data)

xx.fe.re<-na.omit(xx.fe.re)
nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$peak.freq <- log(model.data$peak.freq + 1)
#log transformed response! Better model fit and lower estimates and SE



full.mod.4<-glmmTMB(peak.freq~dummy*stage+z.visibility+parents.call.initial+
                      (1|nest.ID),
                    data=model.data, family=gaussian(link="identity"))

hist(model.data$peak.freq)

mod.4<-lm(peak.freq~dummy+stage+pres.order+z.visibility+parents.call.initial, data=model.data)

vif_values <- vif(mod.4)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.071*1.071
max.vif

testDispersion(full.mod.4)
simulationOutput <- simulateResiduals(fittedModel = full.mod.4, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.4)

null.mod.4<-glmmTMB(peak.freq~pres.order+z.visibility+parents.call.initial+
                      (1|nest.ID),
                    data=model.data, family=gaussian(link="identity"))


as.data.frame(anova(full.mod.4, null.mod.4, test="Chisq"))

#no difference between full and null

summary(full.mod.4)


#total time in (of at least 1 parent)####

xx.fe.re<-fe.re.tab(fe.model="total.time.in~dummy*stage+pres.order+visibility+parents.pres+number.nb", 
                    re="(1|nest.ID)+(1|date)",
                    data=data)

xx.fe.re<-na.omit(xx.fe.re)
nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)


model.data$total.time.in.prop<-(model.data$total.time.in/300)

model.data$total.time.in.prop = (model.data$total.time.in.prop*(nrow(data)-1) + 0.5)/nrow(data)

plot(model.data$total.time.in, model.data$total.time.in.prop)

range(model.data$total.time.in)
#new range of proportions 0.0362-0.9969 - data was squeezed a bit to the centre to avoid 0 and 1
range(model.data$total.time.in.prop)

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))


full.mod.5<-glmmTMB(total.time.in.prop~dummy*stage+pres.order+z.visibility+z.number.nb+parents.pres+
                      (1|nest.ID),
                    data=model.data, family=beta_family(link="logit"))

mod.5<-lm(total.time.in.prop~dummy+stage+pres.order+z.visibility+z.number.nb+parents.pres, data=model.data)

vif_values <- vif(mod.5)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.041*1.041
max.vif

testDispersion(full.mod.5)
simulationOutput <- simulateResiduals(fittedModel = full.mod.5, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.5)


null.mod.5<-glmmTMB(total.time.in.prop~pres.order+z.visibility+z.number.nb+parents.pres+
                      (1|nest.ID),
                    data=model.data, family=beta_family(link="logit"))


as.data.frame(anova(full.mod.5, null.mod.5, test="Chisq"))
summary(full.mod.5)
#sign. difference between full and null but interaction no sign. effect-> sub.mod

sub.mod.5<-glmmTMB(total.time.in.prop~dummy+stage+pres.order+z.visibility+z.number.nb+parents.pres+
                      (1|nest.ID),
                    data=model.data, family=beta_family(link="logit"))

testDispersion(sub.mod.5)
simulationOutput <- simulateResiduals(fittedModel = sub.mod.5, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(sub.mod.5)

summary(sub.mod.5)


r_sq_p = performance::r2(model = full.mod.5) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.5)
r_sq_p
r_sq_m

plot(allEffects(sub.mod.5))

#PLOT
ggplot(model.data, aes(x=stage, y=total.time.in.prop))+
  geom_boxplot(aes(fill=dummy))

# plot raw data with model estimates
# get model estimates
plot(ggpredict(sub.mod.5, terms = c("dummy", "stage")))
xxp = ggpredict(sub.mod.5, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

p <- ggplot(df1, aes(x = dummy, y = total.time.in.prop)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Total time within radius", fill="dummy", color="dummy", title = "Total time in by at least 1 parent")
ggsave("Total_time_in.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

#

p.1 <- ggplot(df1, aes(x = dummy, y = total.time.in.prop)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3, show.legend = FALSE)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  geom_point(aes(x=dummy, y=predicted, fill=dummy), color = "black",size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE, show.legend = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Total time within radius", fill="dummy", color="dummy", title = "Total time in by at least 1 parent")
print(p.1)
ggsave("Total_time_in_simplified.png", plot = p.1, width = 12, height = 8, units = "in", dpi = 300)

#high - low pattern #### 
#NOT PRESENTED IN THE MANUSCRIPT. FURTHER STUDY INTO THIS IS NEEDED.

#HIGH LOW PATTERN
#(yes/no) - no relationship with dummy and stage, random nestID
#model too complex if nestID was included as covariate

full.mod.6<-glmmTMB(hilo.pattern~dummy*stage+(1|nest.ID),
                    data=data, family=binomial(link="logit"))
summary(full.mod.6)


#individuals present####

xx.fe.re<-fe.re.tab(fe.model="tot.ind.pres~dummy*stage+pres.order+visibility+number.nb", 
                    re="(1|nest.ID)",
                    data=data)

#xx.fe.re<-na.omit(xx.fe.re)
nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

range(model.data$tot.ind.pres)

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))


full.mod.6<-glmmTMB(tot.ind.pres~dummy*stage+pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=gaussian(link="identity"))

mod.6<-lm(tot.ind.pres~dummy+stage+pres.order+z.visibility+z.number.nb, data=model.data)


vif_values <- vif(mod.6)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.016*1.016
max.vif

hist(model.data$tot.ind.pres)
range(model.data$tot.ind.pres)

testDispersion(full.mod.6)
simulationOutput <- simulateResiduals(fittedModel = full.mod.6, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.6)


null.mod.6<-glmmTMB(tot.ind.pres~pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=gaussian(link="identity"))

as.data.frame(anova(full.mod.6, null.mod.6, test="Chisq"))

summary(full.mod.6)
#no effect of interaction

r_sq_p = performance::r2(model = full.mod.6) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.6)
r_sq_p
r_sq_m

sub.mod.6<-glmmTMB(tot.ind.pres~dummy+stage+pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    data=model.data, ziformula=~1, family=gaussian(link="identity"))

testDispersion(sub.mod.6)
simulationOutput <- simulateResiduals(fittedModel = sub.mod.6, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(sub.mod.6)

plot(allEffects(sub.mod.6))

summary(sub.mod.6)

emmeans_results <- emmeans(sub.mod.6, ~ dummy |stage)
pairwise_comparisons <- pairs(emmeans_results)
print(pairwise_comparisons)

#export emmeans summary and pairwise comparisons summaries
emm_df <- as.data.frame(emmeans_results)
pairwise_comparison_df <- as.data.frame(pairwise_comparisons)
write_xlsx(pairwise_comparison_df, path = "Pairwise_comparison_Max.indv.present.xlsx")
write_xlsx(emm_df, path = "Emmeans_Max.indv.present.xlsx")


#PLOT
ggplot(model.data, aes(x=dummy, y=tot.ind.pres))+
  geom_boxplot(aes(fill=dummy))

plot(ggpredict(sub.mod.6, terms = c("dummy", "stage")))
xxp = ggpredict(sub.mod.6, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

p <- ggplot(df1, aes(x = dummy, y = tot.ind.pres)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of individuals", fill="dummy", color="dummy", title = "Total number of individuals present within radius")

ggsave("Max.indv.present.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

p.1 <- ggplot(df1, aes(x = dummy, y = tot.ind.pres)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3, show.legend = FALSE)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy), color="black", size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE, show.legend = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Number of individuals", fill="dummy", color="dummy", title = "Total number of individuals present within radius")
print(p.1)
ggsave("Max.indv.present_simplified.png", plot = p.1, width = 12, height = 8, units = "in", dpi = 300)



#parents alarm####

#manually done in Excel only (0/1)= parents.alarm.1; 0= No parents calling, 1 = Either 1 or 2 parents calling
xx.fe.re<-fe.re.tab(fe.model="parents.alarm.1~dummy*stage+pres.order+visibility+number.nb+parents.pres", 
                    re="(1|nest.ID)",
                    data=data)

xx.fe.re<-na.omit(xx.fe.re)
nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

model.data$z.visibility<-as.vector(scale(model.data$visibility))
model.data$z.number.nb<-as.vector(scale(model.data$number.nb))


full.mod.8<-glmmTMB(parents.alarm.1~dummy*stage+pres.order+z.visibility+z.number.nb+
                      (1|nest.ID),
                    data=model.data, family=binomial(link="logit"))

mod.8<-lm(as.numeric(parents.alarm.1)~dummy+stage+pres.order+z.visibility+number.nb, data=model.data)


vif_values <- vif(mod.8)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.017*1.017
max.vif

testDispersion(full.mod.8)
simulationOutput <- simulateResiduals(fittedModel = full.mod.8, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.8)


null.mod.8<-glmmTMB(parents.alarm.1~pres.order+z.visibility+number.nb+
                      (1|nest.ID),
                    data=model.data, family=binomial(link="logit"))

as.data.frame(anova(full.mod.8, null.mod.8, test="Chisq"))

summary(full.mod.8)

r_sq_p = performance::r2(model = full.mod.8) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.8)
r_sq_p
r_sq_m

emm_df <- emmeans(full.mod.8, ~ dummy * stage, method = "pairwise")
pairwise_comparisons_df <- pairs(emm)
print(pairwise_comparisons)

#export emmeans summary and pairwise comparisons summaries
emm_df <- as.data.frame(emmeans_results)
pairwise_comparison_df <- as.data.frame(pairwise_comparisons)
write_xlsx(pairwise_comparison_df, path = "Pairwise_comparison_Parents_AC_1-2.xlsx")
write_xlsx(emm_df, path = "Emmeans_Parents_AC_1-2.xlsx")

#sign. difference between ani building and owl building and between control building and owl building
#no difference between sgf and owl

model.data$predicted_prob <- predict(full.mod.8, type = "response")
#summarize these predicted probabilities to reflect the proportion of 1s for each combination of dummy and stage.
# Summarize predicted probabilities by dummy and stage
summary_table <- model.data %>%
  group_by(dummy, stage) %>%
  summarise(
    mean_predicted_prob = mean(predicted_prob),
    count = n()
  )

print(summary_table)

all_combinations <- expand.grid(dummy = levels(model.data$dummy),
                                stage = levels(model.data$stage))

model.data_complete <- merge(all_combinations, model.data, by = c("dummy", "stage"), all.x = TRUE)

summary_table <- model.data_complete %>%
  group_by(dummy, stage) %>%
  summarise(
    mean_predicted_prob = mean(predicted_prob, na.rm = TRUE),
    median_predicted_prob = median(predicted_prob, na.rm = TRUE),
    min_predicted_prob = min(predicted_prob, na.rm = TRUE),
    max_predicted_prob = max(predicted_prob, na.rm = TRUE),
    count = sum(!is.na(predicted_prob)),  # Number of observations
    .groups = 'drop'  # Drop the grouping structure
  )

print(summary_table)

###
#plots

#plot of predicted probabilities
ggplot(model.data, aes(x = dummy, y = parents.alarm.1, fill = stage)) +
  geom_boxplot() +
  labs(x = "Dummy", y = "Predicted Probability") +
  theme_classic()

data$parents.alarm.1<-as.factor(data$parents.alarm.1)
data$parents.alarm<-as.factor(data$parents.alarm)


#plot parents.alarm.1
p <- ggplot(data, aes(x = dummy, fill = parents.alarm.1)) +
  geom_bar(position = "dodge") +
  labs(x = "dummy type", y = "Count", title = "Parents Alarm Calling during presentation") +
  facet_wrap(~ stage, scales = "free_x") +  # Facet by stage
  theme_classic()
ggsave("Parents_Alarm_0.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

#plot parents.alarm
p <- ggplot(data, aes(x = dummy, fill = parents.alarm)) +
  geom_bar(position = "dodge") +
  labs(x = "dummy type", y = "Count", title="Number of parents Alarm Calling during the presentation", fill ="Number of parents Alarm Calling") +
  facet_wrap(~ stage, scales = "free_x") +  # Facet by stage
  theme_classic()
ggsave("Parents_Alarm.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)


plot(allEffects(full.mod.8))

#PLOT
ggplot(model.data, aes(x=dummy, y=parents.alarm.1))+
  geom_boxplot(aes(fill=dummy))

plot(ggpredict(full.mod.8, terms = c("dummy", "stage")))
xxp = ggpredict(full.mod.8, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(df1, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(dt1, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

ggplot(df1, aes(x = dummy, y = parents.alarm.1)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = conf.high,ymin = conf.low,x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=predicted, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="parents.alarm (0/1)", fill="dummy", color="dummy")
ggsave("Alarm.parents.png", plot = p, width = 12, height = 8, units = "in", dpi = 300)

#high standard errors and large estimates are indicative for complete- or quasi complete separation in dataset
#(Ani and Owl feeding have no 0 cases!! and 100% of data in this combination are response)
#negative influence on interpretability and trustfullness of model results
#solution for report:show barplot (row 980) and put model output to supplementary material.
#state in text that in all cases when parents were confronted with Ani or Owl during feeding both parents were alarming.
#thus we assume that in feeding stage the birds response in alarming is higher (100%) than in building when these dummies were presented
#calculate non-parametric test (e.g. wilcoxon) - however, this doesn?t include nestID and influence of other covariates (especially pres.order).

# Chi-Square Test for dummy
table_dummy <- table(model.data$dummy, model.data$parents.alarm.1)
print(table_dummy)

chisq_dummy <- chisq.test(table(model.data$dummy, model.data$parents.alarm.1))
print(chisq_dummy)

#expected frequencies
chisq_test_dummy <- chisq.test(table_dummy)
print(chisq_test_dummy$expected)

boxplot(model.data$parents.alarm.1~model.data$dummy)

# Pairwise comparisons for dummy
pairwise_dummy <- pairwiseNominalIndependence(table(model.data$dummy, model.data$parents.alarm.1),
                                              method = "bonferroni")
print(pairwise_dummy)
#differences between all other dummies and owl (look at p. adj.Chisq)


#test dummy at stages independently
#building stage
building.data<-subset(model.data, stage=="Building")
contingency_table <- table(building.data$parents.alarm.1, building.data$dummy)
contingency_table#values below 5, thus fisher exact test better than chi-squared
fisher_test_result <- fisher.test(contingency_table)

levels_dummy <- unique(building.data$dummy)
results <- list()

# Loop through pairs of levels
for (i in 1:(length(levels_dummy) - 1)) {
  for (j in (i + 1):length(levels_dummy)) {
    group1 <- levels_dummy[i]
    group2 <- levels_dummy[j]
    
    # Subset the data for the current pair of groups
    subset_data <- building.data %>%
      filter(dummy %in% c(group1, group2))
    
    # Create a contingency table for the pair
    contingency_table <- table(
      subset_data$parents.alarm.1,
      subset_data$dummy
    )
    
    # Perform Fisher's Exact Test
    fisher_result <- fisher.test(contingency_table)
    
    # Store the results
    results[[paste(group1, "vs", group2)]] <- fisher_result$p.value
  }
}

results_df <- data.frame(
  Comparison = names(results),
  p.value = unlist(results)
)
results_df$p.adj <- p.adjust(results_df$p.value, method = "bonferroni")

print(results_df)
#no difference between all dummies and sgf, but between ani and owl and control and owl

#feeding stage
feeding.data<-subset(model.data, stage=="Feeding")
contingency_table <- table(feeding.data$parents.alarm.1, feeding.data$dummy)
contingency_table#values below 5, thus fisher exact test better than chi-squared
fisher_test_result <- fisher.test(contingency_table)
print(fisher_test_result)

levels_dummy <- unique(feeding.data$dummy)
results <- list()

# Loop through pairs of levels
for (i in 1:(length(levels_dummy) - 1)) {
  for (j in (i + 1):length(levels_dummy)) {
    group1 <- levels_dummy[i]
    group2 <- levels_dummy[j]
    
    # Subset the data for the current pair of groups
    subset_data <- feeding.data %>%
      filter(dummy %in% c(group1, group2))
    
    # Create a contingency table for the pair
    contingency_table <- table(
      subset_data$parents.alarm.1,
      subset_data$dummy
    )
    
    # Perform Fisher's Exact Test
    fisher_result <- fisher.test(contingency_table)
    
    # Store the results
    results[[paste(group1, "vs", group2)]] <- fisher_result$p.value
  }
}

results_df <- data.frame(
  Comparison = names(results),
  p.value = unlist(results)
)
results_df$p.adj <- p.adjust(results_df$p.value, method = "bonferroni")

print(results_df)
#only difference between control and owl and control and ani (p.adj.Fisher)

# Chi-Square Test for stage
table_stage <- table(model.data$stage, model.data$parents.alarm.1)
print(table_stage)

chisq_stage <- chisq.test(table(model.data$stage, model.data$parents.alarm.1))
print(chisq_stage)

chisq_test_stage <- chisq.test(table_stage)
print(chisq_test_stage$expected)

boxplot(model.data$parents.alarm.1~model.data$stage)


#test for the relationship between parents.alarm.1 and the number of parents present
#significant relationship: when both parents were present it was more likely that they were alarming and only if only 1 parent was present no alarming occured
model.data$parents.pres
contingency_table <- table(model.data$parents.alarm.1, model.data$parents.pres)
chi_square_result <- chisq.test(contingency_table)
print(chi_square_result)

boxplot(model.data$parents.alarm.1~model.data$parents.pres)
contingency_table <- table(model.data$parents.alarm.1, model.data$parents.pres)
contingency_df <- as.data.frame(contingency_table)

# Plot a heatmap of the contingency table
ggplot(contingency_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "parents alarming", y = "parents present", fill = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#non-parametric test for the relationship between parents.alarm.1 and visibility and number of neighbours
#no significant effect for both predictors
wilcox_test.1 <- wilcox.test(z.visibility ~ parents.alarm.1, data = model.data)
print(wilcox_test.1)

ggplot(model.data, aes(x = parents.alarm.1, y = z.visibility)) +
  geom_boxplot() +
  labs(x = "Parents Alarm (1/0)", y = "Visibility") +
  theme_minimal()

wilcox_test.2 <- wilcox.test(z.number.nb ~ parents.alarm.1, data = model.data)
print(wilcox_test.2)

ggplot(model.data, aes(x = parents.alarm.1, y = z.number.nb)) +
  geom_boxplot() +
  labs(x = "Parents Alarm (1/0)", y = "Number of neighbours (within 20m radius)") +
  theme_minimal()

#parents present####

xx.fe.re<-fe.re.tab(fe.model="parents.pres~dummy*stage+pres.order+visibility+number.nb+parents.alarm", 
                    re="(1|nest.ID)",
                    data=data)

#xx.fe.re<-na.omit(xx.fe.re)
nrow(data)
nrow(xx.fe.re$data)
#identifying theoretically possible slopes
xx.fe.re$summary
model.data=xx.fe.re$data
table(model.data$dummy, model.data$nest.ID)
str(model.data)

full.mod.9<-glmmTMB(as.factor(parents.pres)~dummy*stage+pres.order+
                      (1|nest.ID),
                    data=model.data, family=binomial(link="logit"))

mod.9<-lm(as.numeric(parents.pres)~dummy+stage+pres.order, data=model.data)


vif_values <- vif(mod.9)
vif_values

#to calculate max. vif = (highest value of right column)?
max.vif<-1.003*1.003
max.vif

testDispersion(full.mod.9)
simulationOutput <- simulateResiduals(fittedModel = full.mod.9, plot = F)
residuals(simulationOutput)
residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
plot(simulationOutput)
testZeroInflation(simulationOutput)

ranef.diagn.plot(full.mod.9)


null.mod.9<-glmmTMB(as.factor(parents.pres)~dummy+stage+pres.order+
                      (1|nest.ID),
                    data=model.data, family=binomial(link="logit"))

as.data.frame(anova(full.mod.9, null.mod.9, test="Chisq"))

summary(full.mod.9)

r_sq_p = performance::r2(model = full.mod.9) 
r_sq_m = MuMIn::r.squaredGLMM(object = full.mod.9)
r_sq_p
r_sq_m

#manually change reference order to compare groups
#model.data$pres.order = relevel(model.data$pres.order, ref = "1")
#significant effect for presentation order, wind, hour

emm <- emmeans(full.mod.9, ~ dummy | stage)

# Conduct pairwise comparisons
pairwise_comparisons <- contrast(emm, method = "pairwise")
summary(pairwise_comparisons)

plot(allEffects(full.mod.9))

emm_df <- as.data.frame(emm)

# Plot the estimated marginal means
ggplot(emm_df, aes(x = dummy, y = emmean, ymin = asymp.LCL, ymax = asymp.UCL, color = stage)) +
  geom_pointrange() +
  labs(title = "Estimated Marginal Means",
       x = "dummy",
       y = "Estimated Marginal Means") +
  theme_minimal()


#PLOT
ggplot(model.data, aes(x=dummy, y=parents.pres))+
  geom_boxplot(aes(fill=dummy))

plot(ggpredict(full.mod.9, terms = c("dummy", "stage")))
xxp = ggpredict(full.mod.9, terms = c("dummy", "stage"))
xxp
dt1 = xxp
head(dt1)
str(dt1)
setDT(dt1)
head(dt1)
setnames(dt1, old = c("x", "group"), new = c("dummy", "stage"))
dt1

df1 = model.data

# set dominant to reference level
df1 <- within(emm_df, stage <- relevel(as.factor(stage), ref = "Building"))
dt1 <- within(emm_df, dummy <- relevel(as.factor(dummy), ref = "Ani"))

# plot

ggplot(emm_df, aes(x = dummy, y = emmean)) +
  geom_boxplot(aes(fill=dummy), outlier.shape = NA, position=position_dodge(.9), alpha=0.3)+
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  new_scale_fill()+ new_scale_color()+
  #error plot using data of post-hoc test dt1 (confidence intervals as whiskers)
  geom_errorbar(mapping = aes(ymax = asymp.UCL, ymin = asymp.LCL, x = dummy,color = dummy), width=0.9, size=1.2,
                data = dt1, inherit.aes = FALSE)+
  #point plot using data of post-hoc test dt1 (estimate as diamond)  
  geom_point(aes(x=dummy, y=emmean, fill=dummy, color=dummy), size=3.5, alpha=1, shape=23,
             data = dt1, inherit.aes = FALSE) +
  #scale_fill_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  #scale_color_manual(name="Top-ranked", values=c("#AA3377", "#4477AA"))+
  new_scale_fill()+
  #point plot using df1 model.data
  geom_point(aes(fill=dummy), position=position_jitterdodge(0.9), size=1.5, alpha=0.8, shape=21) +
  #scale_fill_manual(name="Top-ranked", values=c("#EE6677", "#66CCEE"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(angle = -0))+
  facet_wrap(~stage,  strip.position="bottom")+
  labs(x=NULL, y="Parents present", fill="dummy", color="dummy")

ggplot(filtered_data, aes(x = stage, fill = parents.pres)) +
  geom_bar(position = "dodge") +
  labs(x = "parents.pres", y = "Count") +
  theme_minimal()

ggplot(filtered_data, aes(x = dummy, fill = parents.pres)) +
  geom_bar(position = "dodge") +
  labs(x = "parents.pres", y = "Count") +
  theme_minimal()
