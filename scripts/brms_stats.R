library(BayesFactor)
library(brms)
library(tidyverse)

setwd("D:/human_signaling/bayesian_stats")
set.seed(123)

##### D7 #####

### D7 EPIBLAST NUMBERS ###

D7_EPI <- read_csv("D7_EPI.csv")
D7_EPI$Treatment = factor(D7_EPI$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Count ~ Treatment, data=D7_EPI, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
set.seed(123)
control<-subset(D7_EPI, D7_EPI$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)


### D7 TOTAL HYPOBLAST NUMBERS ###

D7_HypoAll <- read_csv("D7_HypoAll.csv")
D7_HypoAll$Treatment = factor (D7_HypoAll$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Count ~ Treatment, data=D7_HypoAll, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
control<-subset(D7_HypoAll, D7_HypoAll$Treatment=='Control')
sample_size = floor(0.75*nrow(control))

picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)

### D7 CER1+ NUMBER ###

D7_CER1 <- read_csv("D7_CER1.csv")
D7_CER1$Treatment = factor (D7_CER1$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Counts ~ Treatment, data=D7_CER1, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
control<-subset(D7_CER1, D7_CER1$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Counts ~ Treatment, data=control, family=poisson())
summary(fit1)

##### DAY 9 #####
set.seed(1234)

D9_EPI <- read_csv("D9_EPI.csv")
D9_EPI$Treatment = factor(D9_EPI$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Count ~ Treatment, data=D9_EPI, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test

control<-subset(D9_EPI, D9_EPI$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)

### D9 TOTAL HYPOBLAST NUMBERS ###

D9_HypoAll <- read_csv("D9_HypoAll.csv")
D9_HypoAll$Treatment = factor (D9_HypoAll$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Count ~ Treatment, data=D9_HypoAll, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
control<-subset(D9_HypoAll, D9_HypoAll$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)

### D9 CER1+ NUMBER ###

D9_CER1 <- read_csv("D9_CER1.csv")
D9_CER1$Treatment = factor (D9_CER1$Treatment, levels=c("Control", "A83", "ActA-25", "LDN", "BMP6", "DAPT", "CpdE", "MK"))

fit1<-brm(formula=Count ~ Treatment, data=D9_CER1, family=poisson())
summary(fit1)


# randomly split control data for negative control stats test
control<-subset(D9_CER1, D9_CER1$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)

###### TrB Fits #####

### D7 All GATA3+ NUMBERS ###

D7_TrBAll <- read_csv("D7_TrBAll.csv")
D7_TrBAll$Treatment = factor(D7_TrBAll$Treatment, levels=c("Control", "A83", "ActA-25"))

fit1<-brm(formula=Count ~ Treatment, data=D7_TrBAll, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
control<-subset(D7_TrBAll, D7_TrBAll$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)

### D7 GATA3+/GATA6+ NUMBERS ###

D7_TrBG6pos <- read_csv("D7_TrBG6pos.csv")
D7_TrBG6pos$Treatment = factor(D7_TrBG6pos$Treatment, levels=c("Control", "A83", "ActA-25"))

fit1<-brm(formula=Count ~ Treatment, data=D7_TrBG6pos, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
set.seed(9876)
control<-subset(D7_TrBG6pos, D7_TrBG6pos$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)


### D7 GATA3+/GATA6- NUMBERS ###

D7_TrBG6neg <- read_csv("D7_TrBG6neg.csv")
D7_TrBG6neg$Treatment = factor(D7_TrBG6neg$Treatment, levels=c("Control", "A83", "ActA-25"))

fit1<-brm(formula=Count ~ Treatment, data=D7_TrBG6neg, family=poisson())
summary(fit1)

# randomly split control data for negative control stats test
set.seed(123)
control<-subset(D7_TrBG6neg, D7_TrBG6neg$Treatment=='Control')
sample_size = floor(0.75*nrow(control))


picked = sample(seq_len(nrow(control)),size = sample_size)
control_1 =control[picked,]
control_2 =control[-picked,]
control_1$Treatment='control_1'
control_2$Treatment='control_2'
control<-rbind(control_1, control_2)
control$Treatment=factor(control$Treatment, levels=c("control_1", "control_2"))

fit1<-brm(formula=Count ~ Treatment, data=control, family=poisson())
summary(fit1)
