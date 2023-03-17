#load library
library(vegan)
library("DAtest")
#load in data
bacteria=read.table("~/Desktop/Microbiome_resubmission/counts.gtdbtk.tsv",header=T,sep="\t")
#bacteria=preDA(bacteria, min.samples = 16, min.reads = 20,min.abundance = 0.01)
metadata=read.table("~/Desktop/Microbiome_resubmission/Metadata.txt",header=T,sep="\t")
colnames(bacteria)=gsub("MG_","",colnames(bacteria))
samples=intersect(colnames(bacteria),metadata$Sample_ID)
bacteria=bacteria[,samples]
bacteria=bacteria[!(rownames(bacteria) %in% c("Unclassified","Unclassified Bacteria")),]
rownames(metadata)=metadata$Sample_ID
metadata=metadata[samples,]
#####Beta diversity analysis
transmission.dist<-vegan::vegdist(t(bacteria), method='bray')
bacteria_data=as.data.frame(t(bacteria))
bacteria_data$Sample_ID=rownames(bacteria_data)
combine=merge(bacteria_data,metadata,by=c("Sample_ID"))
combine$Age=as.numeric(combine$Age)
combine$age_group<-"None"
combine[which(combine$Age>=18),]$age_group<-"Adults"
combine[which(combine$Age<18),]$age_group<-"Children"
combine$age_group=as.factor(combine$age_group)
transmission.div<-adonis2(transmission.dist~ Infection_Household + age_group, data=combine, permutations = 9999, method="bray",block=Subject_ID)
combine$Infection_Household=factor(combine$Infection_Household,levels=c("High ","Low","Control"))
dispersion<-betadisper(transmission.dist, group=combine$Infection_Household)
labs <- paste("Dimension", 1:2, "(", 
              round(100*dispersion$eig / sum(dispersion$eig), 2), "%)")
pdf("~/Desktop/Microbiome_resubmission/transmission_dispersion.pdf")
plot(dispersion, hull=FALSE,xlab=labs[1],ylab=labs[2], ellipse=TRUE,cex.lab=1,col=c("blue","turquoise","black"))
dev.off()
#####differential abundance
library("DAtest")
########MG analysis
transmission_metadata=metadata[colnames(bacteria),]
flu.transmission=as.factor(transmission_metadata$Infection_Household)
flu.transmission=factor(flu.transmission,levels=c("Control","Low","High "))
Household=as.factor(transmission_metadata$Household)
Age=as.numeric(transmission_metadata$Age)
Age=ifelse(Age < 19,"0","1")
Age=as.factor(Age)
Individuals=as.factor(transmission_metadata$Subject_ID)
test <- testDA(bacteria, predictor = flu.transmission, covars = list(Age = Age))
final <- DA.lic(bacteria, paired=Individuals,predictor = flu.transmission, covars = list(Age = Age))
final_output=final[which(final$pval.adj<=0.05),]
features=final_output$Feature
final <- DA.lic(bacteria, paired=Individuals,predictor = flu.transmission, covars = list(Age = Age),allResults=TRUE)
p_values=final$p.value
p_values[,c("predictorLow","predictorHigh ")]=apply(p_values[,c("predictorLow","predictorHigh ")],2,function(m){p.adjust(m,method="fdr")})
p_select=as.data.frame(p_values[features,])
p_select$features=rownames(p_select)
p_select=p_select[which(p_select$predictorLow<=0.05 | p_select$predictorHigh <=0.05),]
p_select=melt(p_select,id.vars=c("features"),measure.vars=c("predictorLow","predictorHigh "))
colnames(p_select)=gsub("value","FDR",colnames(p_select))
p_select=p_select[which(p_select$FDR<=0.05),]

coef=final$coefficients
coef_select=as.data.frame(coef[features,])
coef_select$features=rownames(coef_select)
coef_select=melt(coef_select,id.vars=c("features"),measure.vars=c("predictorLow","predictorHigh "))
colnames(coef_select)=gsub("value","coefficients",colnames(coef_select))

output=merge(p_select,coef_select,by=c("features","variable"))
output$variable=factor(output$variable,levels=c("predictorHigh ","predictorLow"))
output$spp=gsub(".*g__","",output$features)
spp=unique(as.vector(output$spp))

p=ggplot(data=output,aes(x=coefficients,y=reorder(spp,coefficients),fill=variable))+geom_bar(stat="identity")+facet_grid(.~variable)+
  theme_bw()+theme(axis.text.y = element_text(face = "italic"),text = element_text(size = 20),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())+
  labs(x="logFC",y="Bacterial species")+
  scale_fill_manual(values=c("#0000FF","turquoise"))
ggsave("~/Desktop/Microbiome_resubmission/households_MG_comparions.pdf",p,width = 12,height = 12)

######compare between negative samples between H/L and N households
#####beta diversity
metadata_negative=metadata[which(metadata$Influenza_positive=="n"),]
bacteria_negative=bacteria[,metadata_negative$Sample_ID]
bacteria_data=as.data.frame(t(bacteria_negative))
bacteria_data$Sample_ID=rownames(bacteria_data)
combine=merge(bacteria_data,metadata_negative,by=c("Sample_ID"))
combine$Age=as.numeric(combine$Age)
combine$age_group<-"None"
combine[which(combine$Age>=18),]$age_group<-"Adults"
combine[which(combine$Age<18),]$age_group<-"Children"
combine$age_group=as.factor(combine$age_group)
combine$group=combine$Infection_Household
combine$group=gsub("High |Low","Y",combine$group)
combine$group=gsub("Control","N",combine$group)

transmission.dist<-vegan::vegdist(t(bacteria_negative), method='bray')
transmission.div<-adonis2(transmission.dist~ group + age_group, data=combine, permutations = 9999, method="bray",block=Subject_ID)
combine$group=factor(combine$group,levels=c("N","Y"))
dispersion<-betadisper(transmission.dist, group=combine$group)
labs <- paste("Dimension", 1:2, "(", 
              round(100*dispersion$eig / sum(dispersion$eig), 2), "%)")
pdf("~/Desktop/Microbiome_resubmission/transmission_negative_dispersion.pdf")
plot(dispersion, hull=FALSE,xlab=labs[1],ylab=labs[2], ellipse=TRUE,cex.lab=1,col=c("black","red"))
dev.off()
#####DA analysis
metadata_negative=metadata_negative[colnames(bacteria_negative),]
flu.transmission=as.factor(metadata_negative$Infection_Household)
Household=as.factor(metadata_negative$Household)
Age=as.numeric(metadata_negative$Age)
Age=ifelse(Age < 19,"0","1")
Age=as.factor(Age)
Individuals=as.factor(metadata_negative$New_Subject_ID)
transmission=gsub("High |Low","Y",flu.transmission)
transmission=factor(transmission,levels=c("Control","Y"))
test <- testDA(bacteria_negative, predictor = transmission, covars = list(Age = Age))

final <- DA.lic(bacteria_negative, paired=Individuals,predictor = transmission, covars = list(Age = Age))
final_output=final[which(final$pval.adj<=0.05),]
final_output$spp=gsub(".*g__","",final_output$Feature)
final_output=final_output[which(abs(final_output$logFC) >=1),]
final_output$cat="control"
final_output$cat[which(final_output$logFC >0)]<-"Y"
final_output$cat=factor(final_output$cat,levels=c("control","Y"))
p=ggplot(data=final_output,aes(x=logFC,y=reorder(spp,logFC),fill=cat))+geom_bar(stat="identity")+
  theme_bw()+theme(axis.text.y = element_text(face = "italic"),text = element_text(size = 20),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.text.x = element_blank(),strip.text.y = element_blank(),strip.background = element_blank())+
  labs(x="logFC",y="Bacterial species")+scale_fill_manual(values=c("grey","red"))
ggsave("~/Desktop/Microbiome_resubmission/households_negative_samples_MG.pdf",p,width = 12,height = 12)
