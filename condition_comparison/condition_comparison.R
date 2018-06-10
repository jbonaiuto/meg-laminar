beta_rebound_data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/resp_beta_rebound_data.csv',header=TRUE,sep=',')

beta_rebound_data_agg <- aggregate(beta_rebound_data$Power, by=list(Subject=beta_rebound_data$Subject, Congruence=beta_rebound_data$Congruence), FUN=mean)

wilcox.test(beta_rebound_data_agg$x[beta_rebound_data_agg$Congruence=='congruent'], beta_rebound_data_agg$x[beta_rebound_data_agg$Congruence=='incongruent'],, paired=TRUE)



dots_beta_erd_data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/dots_beta_erd_data.csv',header=TRUE,sep=',')

dots_beta_erd_data_agg <- aggregate(dots_beta_erd_data$Power, by=list(Subject=dots_beta_erd_data$Subject, Coherence=dots_beta_erd_data$Coherence), FUN=mean)

dots_beta_erd_data_agg$Coherence<-factor(dots_beta_erd_data_agg$Coherence, levels=c("low","med","high"))

friedman.test(x ~ Coherence|Subject, data = dots_beta_erd_data_agg)



instr_gamma_data<-read.csv(file='/home/bonaiuto/meg_laminar/derivatives/spm12/instr_gamma_data.csv',header=TRUE,sep=',')

instr_gamma_data_agg <- aggregate(instr_gamma_data$Power, by=list(Subject=instr_gamma_data$Subject, Congruence=instr_gamma_data$Congruence), FUN=mean)

wilcox.test(instr_gamma_data_agg$x[instr_gamma_data_agg$Congruence=='congruent'], instr_gamma_data_agg$x[instr_gamma_data_agg$Congruence=='incongruent'],, paired=TRUE)

