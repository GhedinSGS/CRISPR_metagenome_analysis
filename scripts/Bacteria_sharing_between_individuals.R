#load library
library("ggplot2")
library("RColorBrewer")
library("data.table")
#bacterial MAGs have spacers
spacers_analysis=read.table("~/Desktop/Microbiome_resubmission/crispr_array_hits.tsv",header=T,sep="\t")
metadata=read.table("~/Desktop/Microbiome_resubmission/Metadata.txt",header=T,sep="\t")
spacers_analysis$sample=gsub("MG_","",spacers_analysis$sample)
colnames(spacers_analysis)=gsub("sample","Sample_ID",colnames(spacers_analysis))
spacers_analysis=merge(spacers_analysis,metadata,by=c("Sample_ID"))
spacers_analysis$bacteria=gsub(".*g__","g__",spacers_analysis$bacteria)
#spacers clustered together at 90% similarities. 
Data=read.table("~/Desktop/Draft\ paper\ 2021/analysis\ data/spacers_pool_rename_cluster_results",header=FALSE,sep="\t")
Data$V4=gsub(" >","",Data$V4)
colnames(Data)=gsub("V4","spacer",colnames(Data))
######merge data tables 
spacer_output=merge(Data,spacers_analysis,by=c("spacer"))
cluster_shared=data.frame(table(spacer_output$V1))
cluster_shared=as.vector(cluster_shared[which(cluster_shared$Freq >1),]$Var1)
spacer_output=spacer_output[which(spacer_output$V1 %in% cluster_shared),]
new_d=list()
spacer_output$V2<-NULL
spacer_output$V3<-NULL

spacer_output$New_Sample_ID<-NULL
spacer_output$Age<-NULL
spacer_output$Influenza_positive<-NULL
spacer_output$Gender<-NULL
spacer_output$Subject_ID<-NULL

########identify the shared bacterial MAGs with shared spacers between samples
for(i in 1:length(cluster_shared)){
  d=spacer_output[which(spacer_output$V1==cluster_shared[i]),]
  a=dim(d)[1]
  data=expand.grid((1:a),(1:a))
  data=data[which(data$Var1!=data$Var2),]
  data=unique(t(apply(data, 1, sort)))
  data=as.data.frame(data)
  new_d[[i]]=cbind(d[data[,1],],d[data[,2],])
  colnames(new_d[[i]])=c("spacer_1","cluster_1","sample_ID_1","bacteria_1","virus_1","HH_1","Subject_1","Infection_HH_1","spacer_2","cluster_2","sample_ID_2","bacteria_2","virus_2","HH_2","subject_2","Infection_HH_2")
}

new_d=data.frame(rbindlist(new_d))

#####filter the spp by number of spacers mapped to the bacterial MAGs
spacers_analysis$genus=gsub(";s__.*","", spacers_analysis$bacteria)
spacers_analysis$genus=gsub("g__Haemophilus_.*","g__Haemophilus", spacers_analysis$genus)
spacers_analysis$count="1"
spacers_analysis$count=as.numeric(spacers_analysis$count)
spp_prevalence=aggregate(count ~ Sample_ID+bacteria, spacers_analysis, sum)
spp_filter=spp_prevalence[which(spp_prevalence$count >=3),]
spp_filter$count<-NULL
######plot for prevalence using spp filter
colnames(spp_filter)=c("sample_ID_1","bacteria_1")
new_d_spp=merge(new_d,spp_filter,by=c("sample_ID_1","bacteria_1"))
colnames(spp_filter)=c("sample_ID_2","bacteria_2")
new_d_spp=merge(new_d_spp,spp_filter,by=c("sample_ID_2","bacteria_2"))
new_d_spp$genus=gsub(";s__.*","", new_d_spp$bacteria_1)
new_d_spp$genus=gsub("g__Haemophilus_.*","g__Haemophilus", new_d_spp$genus)
new_d_spp=new_d_spp[-which(new_d_spp$bacteria_1 %in% c("Unclassified","Unclassified Bacteria") | new_d_spp$bacteria_2 %in% c("Unclassified","Unclassified Bacteria")),]
new_d_spp$cat="same_individual"
new_d_spp$cat[which(new_d_spp$HH_1 == new_d_spp$HH_2 & new_d_spp$Subject_1 != new_d_spp$subject_2)]<-"same_household"
new_d_spp$cat[which(new_d_spp$HH_1 != new_d_spp$HH_2)]<-"diff_household"
new_d_spp=new_d_spp[which(new_d_spp$bacteria_1 ==new_d_spp$bacteria_2),]
#####normalize the data by number of samples in each household 
samples=unique(gsub(".*\\|","",Data$spacer))
Metadata_select=metadata[which(metadata$Sample_ID %in% samples),]
Household_counts=data.frame(table(Metadata_select$Household))
Individual_counts=unique(data.frame(Metadata_select$Household,Metadata_select$New_Subject_ID))
Individual_counts=data.frame(table(Individual_counts$Metadata_select.Household))
Total=sum(Household_counts$Freq)
timepoint_counts=data.frame(table(Metadata_select$New_Subject_ID))
average_samples_per_household=15.9
average_person_per_household=5.4
average_timepoint=2.9
prevalence=unique(data.frame(new_d_spp$Subject_1,new_d_spp$subject_2,new_d_spp$bacteria_1,new_d_spp$cat,new_d_spp$genus))
prevalence$count="1"
prevalence$count=as.numeric(prevalence$count)
####remove the samples counted twice
prevalence=prevalence[!duplicated(data.frame(t(apply(prevalence[1:2], 1, sort)), prevalence[3:6])),]

prevalence_plot=aggregate(count ~ new_d_spp.bacteria_1 + new_d_spp.cat +new_d_spp.genus,data=prevalence,sum)
prevalence_plot$normalize_counts="0"
#for(i in 1:nrow(prevalence_plot)){
#if (prevalence_plot$new_d_spp.cat[i] =="diff_household"){#prevalence_plot$normalize_counts[i]=prevalence_plot$count[i]/(10*average_samples_per_household*9*average_samples_per_household)}
#  else
#    {#prevalence_plot$normalize_counts[i]=prevalence_plot$count[i]/(10*average_samples_per_household*(average_person_per_household-1)*average_timepoint)}
#}
for(i in 1:nrow(prevalence_plot)){
  if (prevalence_plot$new_d_spp.cat[i] =="diff_household"){prevalence_plot$normalize_counts[i]=prevalence_plot$count[i]/(10*average_person_per_household*9*average_person_per_household)}
  else
  {prevalence_plot$normalize_counts[i]=prevalence_plot$count[i]/(10*average_person_per_household*(average_person_per_household-1))}
}


prevalence_plot=prevalence_plot[which(prevalence_plot$new_d_spp.cat !="same_individual"),]
prevalence_plot$new_d_spp.genus=gsub("g__Leptotrichia_A","g__Leptotrichia",prevalence_plot$new_d_spp.genus)
prevalence_plot$new_d_spp.genus=factor(prevalence_plot$new_d_spp.genus,levels=c("g__Neisseria","g__Leptotrichia","g__Anaeroglobus","g__Fusobacterium","g__Veillonella","g__Pauljensenia","g__Prevotella","g__Rothia"))
prevalence_plot$normalize_counts=as.numeric(prevalence_plot$normalize_counts)
color_code=brewer.pal(n = 8, name = "Dark2")
color_match=data.frame(c("g__Neisseria","g__Leptotrichia","g__Anaeroglobus","g__Fusobacterium","g__Pauljensenia","g__Prevotella","g__Rothia","g__Veillonella"),color_code)
colnames(color_match)=c("new_d_spp.genus","color")
prevalence_plot=merge(prevalence_plot,color_match,by=c("new_d_spp.genus"))
#prevalence_plot$color=reorder(prevalence_plot$color,prevalence_plot$new_d_spp.genus)
prevalence_plot$spp=gsub(".*s__","",prevalence_plot$new_d_spp.bacteria_1)
p=ggplot(data=prevalence_plot,aes(x=normalize_counts,y=reorder(spp, normalize_counts)))+geom_bar(stat="identity",aes(fill=new_d_spp.cat))+facet_grid(new_d_spp.genus~new_d_spp.cat,scale="free",space="free")+xlim(c(0,0.03))+
  theme_bw()+theme(text = element_text(size = 15),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())+
  labs(x="Normalized number of connections",y="Bacterial species")+
  scale_fill_manual(values=c("orange","purple"))

ggsave("~/Desktop/Microbiome_resubmission/Shared_of_bacteria.pdf",p,width = 10,height = 5)



