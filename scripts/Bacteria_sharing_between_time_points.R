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
new_d_spp=new_d_spp[which(new_d_spp$cat!="same_individual"),]
########for each bacteria, we plot out the individuals and time points shared spacers
bacteria=unique(new_d_spp$bacteria_1)
#analzye the time points and spacers
spacers_analysis=read.table("~/Desktop/Microbiome_resubmission/crispr_array_hits.tsv",header=T,sep="\t")
metadata=read.table("~/Desktop/Microbiome_resubmission/Metadata.txt",header=T,sep="\t")
spacers_analysis$sample=gsub("MG_","",spacers_analysis$sample)
colnames(spacers_analysis)=gsub("sample","Sample_ID",colnames(spacers_analysis))
spacers_analysis=merge(spacers_analysis,metadata,by=c("Sample_ID"))
spacers_analysis$bacteria=gsub(".*g__","g__",spacers_analysis$bacteria)
Data=read.table("~/Desktop/Draft\ paper\ 2021/analysis\ data/spacers_pool_rename_cluster_results",header=FALSE,sep="\t")
Data$V4=gsub(" >","",Data$V4)
colnames(Data)=gsub("V4","spacer",colnames(Data))

spacer_output=merge(Data,spacers_analysis,by=c("spacer"))
spacer_output$V1=gsub(">","",spacer_output$V1)
spacer_output$timepoints=gsub(".*\\.","",spacer_output$Sample_ID)
###output the spacers
p=list()
plot.new()
d=list()

for(i in 1:length(bacteria)){
  samples=unique(c(new_d_spp[which(new_d_spp$bacteria_1==bacteria[i]),]$sample_ID_2,new_d_spp[which(new_d_spp$bacteria_1==bacteria[i]),]$sample_ID_1))
  d[[i]]=spacer_output[which(spacer_output$bacteria==bacteria[i] & spacer_output$Sample_ID %in% samples),]
}
p[[1]]=ggplot(d[[1]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[1])
#p[[i]]=ggplot(d,aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[1],"_TP.pdf"),p[[1]],width = 5,height = 2.5)
p[[2]]=ggplot(d[[2]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[2])
#p[[i]]=ggplot(d,aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[2],"_TP.pdf"),p[[2]],width = 12,height = 7)
p[[3]]=ggplot(d[[3]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[3])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[3],"_TP.pdf"),p[[3]],width = 6,height = 3)
p[[4]]=ggplot(d[[4]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[4])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[4],"_TP.pdf"),p[[4]],width = 6,height = 3)
p[[5]]=ggplot(d[[5]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[5])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[5],"_TP.pdf"),p[[5]],width = 4,height = 3)
p[[6]]=ggplot(d[[6]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[6])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[6],"_TP.pdf"),p[[6]],width = 3,height = 2.5)
p[[7]]=ggplot(d[[7]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[7])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[7],"_TP.pdf"),p[[7]],width = 2,height = 1)
p[[8]]=ggplot(d[[8]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[8])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[8],"_TP.pdf"),p[[8]],width = 4,height = 1.5)
p[[9]]=ggplot(d[[9]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[9])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[9],"_TP.pdf"),p[[9]],width = 2.5,height = 1.5)
p[[10]]=ggplot(d[[10]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[10])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[10],"_TP.pdf"),p[[10]],width = 3,height = 1.5)
p[[11]]=ggplot(d[[11]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[11])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[11],"_TP.pdf"),p[[11]],width = 3,height = 1)
p[[12]]=ggplot(d[[12]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[12])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[12],"_TP.pdf"),p[[12]],width = 2,height = 1)
p[[13]]=ggplot(d[[13]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[13])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[13],"_TP.pdf"),p[[13]],width = 2,height = 1)
p[[14]]=ggplot(d[[14]],aes(x=V1,y=Sample_ID))+geom_tile()+facet_grid(V3 ~ ., space="free",scale="free")+ theme_bw()+theme(panel.spacing = unit(0.1, "lines"),axis.title.y=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill = NA, color = NA), panel.border = element_rect(colour = "black", fill = NA),axis.text.x = element_blank(),strip.text.y = element_text(angle = 0))+
  labs(y="samples",x="spacers",title=bacteria[14])
ggsave(paste("~/Desktop/Microbiome_resubmission/",bacteria[14],"_TP.pdf"),p[[14]],width = 2,height = 1)

