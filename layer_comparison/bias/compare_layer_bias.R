library("Matrix")
library("lme4")
library("ggplot2")
library("lsmeans")
library("car")
library("Rmisc")

data<-read.csv(file='/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/bias.csv',header=TRUE,sep=',')

sink('/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/bias_stats.txt')


print('Whole Brain - LFN diff')
wb_lfn_data<-data[data$ROI=='whole_brain' & data$Measure=='LFNDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = wb_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

pw<-lsmeans(model, pairwise~FOI*Type|FOI)
print(summary(pw)$contrasts)

pw<-lsmeans(model, pairwise~FOI*Type|Type)
print(summary(pw)$contrasts)

model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = wb_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

pw<-lsmeans(model, pairwise~FOI*Type|FOI)
print(summary(pw)$contrasts)

pw<-lsmeans(model, pairwise~FOI*Type|Type)
print(summary(pw)$contrasts)


print('')
print('')
print('')
print('Whole Brain - Depth diff')
wb_depth_data<-data[data$ROI=='whole_brain' & data$Measure=='DepthDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = wb_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

pw<-lsmeans(model, pairwise~FOI*Type|FOI)
print(summary(pw)$contrasts)

pw<-lsmeans(model, pairwise~FOI*Type|Type)
print(summary(pw)$contrasts)


model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = wb_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

pw<-lsmeans(model, pairwise~FOI*Type|FOI)
print(summary(pw)$contrasts)

pw<-lsmeans(model, pairwise~FOI*Type|Type)
print(summary(pw)$contrasts)



print('')
print('')
print('')
print('Functional ROI - LFN diff')
func_lfn_data<-data[data$ROI=='func' & data$Measure=='LFNDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = func_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = func_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)


print('')
print('')
print('')
print('Functional ROI - Depth diff')
func_depth_data<-data[data$ROI=='func' & data$Measure=='DepthDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = func_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = func_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

pw<-lsmeans(model, pairwise~FOI*Type|FOI)
print(summary(pw)$contrasts)

pw<-lsmeans(model, pairwise~FOI*Type|Type)
print(summary(pw)$contrasts)





print('')
print('')
print('')
print('Anatomical ROI - LFN diff')
anat_lfn_data<-data[data$ROI=='anat' & data$Measure=='LFNDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = anat_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = anat_lfn_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)


print('')
print('')
print('')
print('Anatomical ROI - Depth diff')
anat_depth_data<-data[data$ROI=='anat' & data$Measure=='DepthDiff',]

model <- lmer(Z ~ FOI*Type+ (1|Subject), data = anat_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

model <- lmer(PartialZ ~ FOI*Type+ (1|Subject), data = anat_depth_data, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
print(summary(model))
results<-Anova(model, type = 3, test = "F")
print(results)

sink()
