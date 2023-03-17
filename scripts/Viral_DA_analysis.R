####load in library
library(vegan)
library("DAtest")
####load in data
vcontact2_phage=read.table("~/Desktop/Microbiome_resubmission/counts.vcontact2.tsv",header=T,sep="\t")
colnames(vcontact2_phage)=gsub("MG_","",colnames(vcontact2_phage))
samples=intersect(colnames(vcontact2_phage),metadata$Sample_ID)
vcontact2_phage=vcontact2_phage[,samples]
vcontact2_phage=preDA(vcontact2_phage,min.samples = 10,min.reads = 10)
###transmission
rownames(metadata)=metadata$Sample_ID
transmission_metadata=metadata[samples,]
flu.transmission=as.factor(transmission_metadata$Infection_Household)
flu.transmission=factor(flu.transmission,levels=c("Control","Low","High "))
Household=as.factor(transmission_metadata$Household)
Age=as.numeric(transmission_metadata$Age)
Age=ifelse(Age < 19,"0","1")
Age=as.factor(Age)
Individuals=as.factor(transmission_metadata$Subject_ID)
final <- DA.lia(vcontact2_phage, paired=Individuals,predictor = flu.transmission, covars = list(Age = Age))
final_output=final[which(final$pval.adj<=0.05),]
features=final_output$Feature
final_output=melt(final_output,id.vars=c("Feature","pval.adj"),measure.vars=c("predictorLow","predictorHigh."))
final_output=final_output[which(abs(final_output$value)>=1),]
final_output$Feature=gsub("g__","",final_output$Feature)
#####barplots
p=ggplot(data=final_output,aes(x=value,y=reorder(Feature,value),fill=variable))+geom_bar(stat="identity")+facet_grid(.~variable)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic"),text = element_text(size = 20),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())+
  labs(x="logFC",y="Viral genus")+
  scale_fill_manual(values=c("#0000FF","turquoise"))
ggsave("~/Desktop/Microbiome_resubmission/viral_DA.pdf",p,width = 6,height = 3)
