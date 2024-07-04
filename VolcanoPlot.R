library(ggplot2)
library(ggrepel)

setwd("D:/OneDrive - University of Nebraska Medical Center/NAS/VolcanoPlot")
###########################################################################
# Diff Methylation volcano plot
# Volcano plot OPIOIDNAS_vc_OPIOIDNoNAS
meth <- read.csv("20-New Groupdata_OPIOIDNAS_vc_OPIOIDNoNAS Nov 5 2020 for Volcono Plot to Nitish.txt",header = TRUE, sep = "\t")
meth.OPIOIDNAS_vc_OPIOIDNoNAS <- meth
Gene_Type <- "methylation"
#genes <- meth
meth1 <- meth[meth$FDR.p.Val < 0.00001,]
genes <-meth1
#list_cpg <- c("cg03544320","cg03692651","cg07915921","cg15811515","cg19717586","cg22620090","cg22674699","cg22784954",
#              "cg22797031","cg00446123","cg02772121","cg07388969","cg11201447","cg13860281","cg14931884", 
#              "cg16389901","cg20765408","cg20852851")

#genes_subset <- genes[genes$X %in% list_cpg,] 
x.cut=5;y.cut=0.01
Significance <- ifelse(meth$Methylation.diff >= x.cut & meth$FDR.p.Val < y.cut, "Hypermethylated", ifelse(meth$Methylation.diff <= -x.cut & meth$FDR.p.Val < y.cut, "Hypomethylated", "Not significant"))

ggplot(meth, aes(x = Methylation.diff, y = -log10(FDR.p.Val))) +
  geom_point(aes(color = Significance), size = 1.2) +
  scale_color_manual(values = c( "blue4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  #geom_text_repel(data = subset(meth, ( FDR.p.Val< 0.00001 & abs(Methylation.diff) >=12)),
   #               aes(label = TargetID)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("Delta Beta Value in Percentage") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = paste0(Gene_Type,"_VolcanoPlot_Groupdata_OPIOIDNAS_vc_OPIOIDNoNAS.pdf"), width = 12, height = 8, dpi = 800)

#####################################################################################
### Volcano plot for OPIOIDusedNAS_OPIOIDuedNONAS-vs_NormalControls
meth <- read.csv("10_Group_DataOPIOIDusedNAS_OPIOIDuedNONAS-vs_NormalControlsNov52020 Nitish volcono.txt",header = TRUE, sep = "\t")
meth$Methylation.diff <- meth$X..diffreence
meth.OPIOIDusedNAS_OPIOIDuedNONAS_vs_NormalControls <- meth
#Gene_Type <- "methylation"
#meth1 <- meth[meth$FDR.p.Val < 0.00001,]
#genes <-meth1
x.cut=5;y.cut=0.01

Significance <- ifelse(meth$Methylation.diff >= x.cut & meth$FDR.p.Val < y.cut, "Hypermethylated", ifelse(meth$Methylation.diff <= -x.cut & meth$FDR.p.Val < y.cut, "Hypomethylated", "Not significant"))

ggplot(meth, aes(x = Methylation.diff, y = -log10(FDR.p.Val))) +
  geom_point(aes(color = Significance), size = 1.2) +
  scale_color_manual(values = c( "red4","blue4" ,"gray41")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  #geom_text_repel(data = subset(meth, ( FDR.p.Val< 0.00001 & abs(Methylation.diff) >=12)),
  #               aes(label = TargetID)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("Delta Beta Value in Percentage") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = paste0(Gene_Type,"_VolcanoPlot_Groupdata_OPIOIDusedNAS_OPIOIDuedNONAS-vs_NormalControls.pdf"), width = 12, height = 8, dpi = 800)

#####################################################################################
### Volcano plot for OPIOIDusedNAS_vs_NormalControls

meth <- read.csv("10_NAS-AUC-Opioid Used NAS vs Normal ControlsMay 19 2020 80249 text volcono to Nitish.txt",header = TRUE, sep = "\t")
meth$Methylation.diff <- meth$X..methylation.Diff
meth.OPIOIDusedNAS_vs_NormalControls <- meth
#Gene_Type <- "methylation"
#meth1 <- meth[meth$FDR.p.Val < 0.00001,]
#genes <-meth1
x.cut=5;y.cut=0.01

Significance <- ifelse(meth$Methylation.diff >= x.cut & meth$FDR.p.Val < y.cut, "Hypermethylated", ifelse(meth$Methylation.diff <= -x.cut & meth$FDR.p.Val < y.cut, "Hypomethylated", "Not significant"))

ggplot(meth, aes(x = Methylation.diff, y = -log10(FDR.p.Val))) +
  geom_point(aes(color = Significance), size = 1.2) +
  scale_color_manual(values = c( "red4","blue4" ,"gray41")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  #geom_text_repel(data = subset(meth, ( FDR.p.Val< 0.00001 & abs(Methylation.diff) >=12)),
  #               aes(label = TargetID)) +
  
  theme(legend.position="right")+
  ggtitle("Volcano plot")+
  xlab("Delta Beta Value in Percentage") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = paste0(Gene_Type,"_VolcanoPlot_Opioid Used NAS vs Normal Control.pdf"), width = 12, height = 8, dpi = 800)

save.image("VolcanoPlot.RData")
