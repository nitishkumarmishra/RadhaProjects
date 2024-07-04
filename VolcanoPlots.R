setwd("C:/Users/siddesh.southekal/University of Nebraska Medical Center/Mishra, Nitish K - NAS/HTSeqCount/OBBM")


set.seed(8456)
library(ComplexHeatmap)
library(circlize)
library(genefilter)
library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)


Group1_vs_4 <- read.csv("Group1_Vs_Group4_DEG.txt",header = TRUE,stringsAsFactors = FALSE)
  
  

x.cut=1.5;y.cut=0.05

Significance <- ifelse(Group1_vs_4$log2FoldChange >= x.cut & Group1_vs_4$padj < y.cut, "Upregulated", ifelse(Group1_vs_4$log2FoldChange <= -x.cut & Group1_vs_4$padj < y.cut, "Downregulated", "Not significant"))

ggplot(Group1_vs_4, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 1) +
  scale_color_manual(values = c( "green4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(Group1_vs_4, (padj < 0.01 & abs(log2FoldChange) >=2)),
                  aes(label = Symbol)) +
  
  theme(legend.position="right")+
  ggtitle("Group1 Vs Group4 Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = "Volcano_Group1_vs_Group4.pdf", width = 12, height = 8, dpi = 800)



#####



Group2_vs_4 <- read.csv("Group2_Vs_Group4_DEG.txt",header = TRUE,stringsAsFactors = FALSE)



x.cut=1.5;y.cut=0.05

Significance <- ifelse(Group2_vs_4$log2FoldChange >= x.cut & Group2_vs_4$padj < y.cut, "Upregulated", ifelse(Group2_vs_4$log2FoldChange <= -x.cut & Group2_vs_4$padj < y.cut, "Downregulated", "Not significant"))

ggplot(Group2_vs_4, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 1) +
  scale_color_manual(values = c( "green4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(Group2_vs_4, (padj < 0.01 & abs(log2FoldChange) >=2)),
                  aes(label = Symbol)) +
  
  theme(legend.position="right")+
  ggtitle("Group2 Vs Group4 Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = "Volcano_Group2_vs_Group4.pdf", width = 12, height = 8, dpi = 800)

####


Group3_vs_4 <- read.csv("Group3_Vs_Group4_DEG.txt",header = TRUE,stringsAsFactors = FALSE)



x.cut=1.5;y.cut=0.05

Significance <- ifelse(Group3_vs_4$log2FoldChange >= x.cut & Group3_vs_4$padj < y.cut, "Upregulated", ifelse(Group3_vs_4$log2FoldChange <= -x.cut & Group3_vs_4$padj < y.cut, "Downregulated", "Not significant"))

ggplot(Group3_vs_4, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Significance), size = 1) +
  scale_color_manual(values = c( "green4","gray41" ,"red4")) +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_vline(aes(xintercept = -x.cut), colour = "black", linetype = "dashed") + 
  geom_vline(aes(xintercept = x.cut),  colour = "black", linetype = "dashed") + 
  geom_hline(aes(yintercept = -1 * log10(y.cut)), colour = "black", linetype = "dashed") +
  
  geom_text_repel(data = subset(Group3_vs_4, (padj < 0.01 & abs(log2FoldChange) >=2)),
                  aes(label = Symbol)) +
  
  theme(legend.position="right")+
  ggtitle("Group3 Vs Group4 Volcano plot")+
  xlab("log2(FoldChange)") + ylab("-log10(padj)") +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.5))+
  theme(legend.key = element_rect(colour = "white", fill = "white"))+ ### remove box from legends
  theme(legend.title=element_text(colour="black", size = 10, face = "bold"))+
  theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+
  theme(axis.title = element_text(color="black", face="bold", size=14)) 

ggsave(filename = "Volcano_Group3_vs_Group4.pdf", width = 12, height = 8, dpi = 800)

