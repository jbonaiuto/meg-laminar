library("Matrix")
library("lme4")
library("ggplot2")
library("lsmeans")
library("car")
library("Rmisc")

data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/behav_data.csv',header=TRUE,sep=',')

model <- glmer(Correct ~ Congruence*Coherence+(1|Subject), data = data, family=binomial, control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 3)
print(results)

twow_pw2<-lsmeans(model, pairwise~Congruence*Coherence|Coherence)
print(summary(twow_pw2)$contrasts)

model <- lmer(RT ~ Congruence*Coherence+(1|Subject), data = data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
results<-Anova(model, type = 3, test = "F")
print(results)

twow_pw2<-lsmeans(model, pairwise~Congruence*Coherence|Coherence)
print(summary(twow_pw2)$contrasts)

sink()
