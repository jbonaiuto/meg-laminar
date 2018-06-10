library("doBy")

data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/dots_alpha/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('Dots-Alpha')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)


data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/dots_gamma/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('Dots-Gamma')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)


data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/instr_gamma/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('Instr-Gamma')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)


data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/dots_beta_erd/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('Beta-ERD')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)


data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/resp_beta_rebound/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('Beta-Rebound')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)


data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/against_shuffled/resp_mrgs/subjects_data.csv',header=TRUE,sep=',')
data$Subject<-as.factor(data$Subject)
agg_data<-summaryBy(Pial.White ~ Subject + ROIType, FUN=mean, data=data, keep.names=TRUE)

print('MRGS')
print('Global')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Global'],mu=0,alternative='two.sided')
print(results)

print('Functional')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Func'],mu=0,alternative='two.sided')
print(results)

print('Anatomical')
results<-wilcox.test(agg_data$Pial.White[agg_data$ROIType=='Anat'],mu=0,alternative='two.sided')
print(results)

sink()
