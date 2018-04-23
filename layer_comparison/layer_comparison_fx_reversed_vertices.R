library("doBy")

sink('/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/layer_comparison_fx_reversed_vertices_stats.txt')

comparisons<-c('dots_alpha','dots_gamma','instr_gamma','dots_beta_erd','resp_beta_rebound','resp_mrgs')
for(i in 1:length(comparisons)) {
	comparison<-comparisons[i]
	data<-read.csv(file=paste0('/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/',comparison,'/bias/reversed_depth/whole_brain/subjects.csv'),header=TRUE,sep=',')
	data$Subject<-as.factor(data$Subject)
	agg_data<-summaryBy(Pial.White ~ Subject, FUN=mean, data=data, keep.names=TRUE)

	print(comparison)
	print('Global')
	m<-mean(agg_data$Pial.White)
	print(m)
	results<-wilcox.test(agg_data$Pial.White,mu=0,alternative='two.sided')
	print(results)

	data<-read.csv(file=paste0('/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/',comparison,'/bias/reversed_depth/func/subjects.csv'),header=TRUE,sep=',')
	data$Subject<-as.factor(data$Subject)
	agg_data<-summaryBy(Pial.White ~ Subject, FUN=mean, data=data, keep.names=TRUE)

	print('Functional')
	m<-mean(agg_data$Pial.White)
	print(m)
	results<-wilcox.test(agg_data$Pial.White,mu=0,alternative='two.sided')
	print(results)

	data<-read.csv(file=paste0('/home/bonaiuto/Dropbox/meg/pred_coding/plots/layer_comparison/',comparison,'/bias/reversed_depth/anat/subjects.csv'),header=TRUE,sep=',')
	data$Subject<-as.factor(data$Subject)
	agg_data<-summaryBy(Pial.White ~ Subject, FUN=mean, data=data, keep.names=TRUE)

	print('Anatomical')
	m<-mean(agg_data$Pial.White)
	print(m)
	results<-wilcox.test(agg_data$Pial.White,mu=0,alternative='two.sided')
	print(results)
}

sink()
