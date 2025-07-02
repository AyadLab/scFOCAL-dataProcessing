# 200nM gene signature characterization

## Load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(msigdbr)
library(ggrepel)
## Load Bioconductor packages
# library(GEOquery)
library(DESeq2)
library(Biobase)
library(fgsea)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(singscore)
library(stringr)
## load data
gse <- read.delim(file = "CTL vs 200nM CT179-rawdata.txt", row.names = 1, check.names = F)

## experiment design
pheno <- data.frame(row.names = colnames(gse))
rownames(pheno)

for (i in 1:length(rownames(pheno))){
  if ((unlist(strsplit(rownames(pheno)[i], split = "-", fixed = T))[2]) == "ctl"){
    pheno$group[i] <- "ctl"
  } else {
    pheno$group[i] <- "CT179"
  }
  }
pheno$replicate <- c("1", "2", "3", "1", "2", "3")

pheno$group <- as.factor(pheno$group)
pheno$replicate <- as.factor(pheno$replicate)
pheno <- pheno[order(match(rownames(pheno), colnames(gse))), ]

pheno

## DEG analysis
dds <- DESeqDataSetFromMatrix(countData = gse, colData = pheno, design = ~ replicate + group)
keep <- rowSums(counts(dds)) >= 10
dim(dds) # before filtering: 35236
dds <- dds[keep,]
dim(dds) # after filtering: 18302

res <- DESeq(dds)
degs <- results(res, contrast = c("group", "CT179", "ctl"))
degs2 <- as.data.frame(degs) #degs before lfc

pdf("Olig2_200nM_volc_before_LFC.pdf")
volc <- EnhancedVolcano::EnhancedVolcano(toptable = degs2, lab = rownames(degs2), x = "log2FoldChange", y = "pvalue", pCutoff = .05, ylim = c(0, max(-log10(degs2[["pvalue"]]), na.rm = TRUE) + 0.5))
volc
dev.off()

## shrinkLFC
resLFC <- lfcShrink(res, contrast = c("group", "CT179", "ctl"), type = "normal")
resLFC_DF <- as.data.frame(resLFC)
resLFC_DF

pdf(file = "200nM_CT179_vs_ctl_volcano.pdf")
volc <- EnhancedVolcano::EnhancedVolcano(toptable = resLFC_DF, 
                                         lab = rownames(resLFC_DF), 
                                         x = "log2FoldChange", 
                                         y = "padj", 
                                         pCutoff = .05, FCcutoff = 0.25,
                                         ylim = c(0, 170), 
                                         xlim = c(min(resLFC_DF[["log2FoldChange"]], na.rm = TRUE) - .5, 
                                                  max(resLFC_DF[["log2FoldChange"]], na.rm = TRUE) + .5))

volc
dev.off()

## Revised volcano plot version
pdf(file = "Revised_200nM_CT179_vs_ctl_volcano.pdf")

logFC_t <- with(resLFC_DF,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange))) 
## 这里使用的是动态阈值，也可以自定义例如logFC=1或2这种静态阈值
logFC_t <- round(logFC_t, 3) # 取前三位小数，这步也可以不运行
logFC_t #看一下动态阈值

resLFC_DF$Change = as.factor(ifelse(resLFC_DF$padj < 0.05 & abs(resLFC_DF$log2FoldChange) > logFC_t,
                                    ifelse(resLFC_DF$log2FoldChange > logFC_t ,'UP','DOWN'),'STABLE'))
## 定义校正p值<0.05和差异倍数大于动态阈值的结果
table(resLFC_DF$Change) #看一下上下调基因的数量

resLFC_DF$label <- ifelse(resLFC_DF$padj< 0.0005& abs(resLFC_DF$log2FoldChange) >= 1.0,rownames(resLFC_DF),"")

ggplot(resLFC_DF, aes(x=log2FoldChange, y=-log10(pvalue),color=Change)) + 
  geom_point(alpha=0.4, size=2) +  # 设置点的透明度和大小
  theme_bw(base_size = 12) +  #设置一个主题背景
  xlab("Log2(Fold change)") + # x轴名字
  ylab("-Log10(P.adj)") + # y轴名字
  theme(plot.title = element_text(size=15,hjust = 0.5)) + 
  scale_colour_manual(values = c('steelblue','gray','brown')) + # 各自的颜色
  geom_hline(yintercept = -log10(0.05), lty = 4) + #定义p值和线形
  geom_vline(xintercept = c(-logFC_t, logFC_t), lty = 4)+ #定义差异倍数和线形
  labs(title = 'volcano')+ #加上题目
  geom_label_repel(data = resLFC_DF, aes(label = label),
                 size = 3,box.padding = unit(0.5, "lines"),
                 point.padding = unit(0.8, "lines"),
                 segment.color = "black",
                 show.legend = FALSE, max.overlaps = 10000)

dev.off()

resLFC_DF <- resLFC_DF[,-8]

## select the significant genes
class(counts(dds))
select <- subset(resLFC_DF, abs(resLFC_DF$log2FoldChange) > logFC_t & resLFC_DF$pvalue < 0.05)
degs_names <- rownames(select)

pheno$ID_1
duplicated <- rownames(pheno[duplicated(pheno$ID_1),])
pheno2 <- pheno
colanno <- data.frame(row.names = rownames(pheno2), group = pheno2$group)

## heatmap
pdf(file = "200nM_CT179_vs_ctl_degs_heatmap.pdf")
pheatmap(log(counts(dds)[degs_names,] + 1), 
         scale = "row", annotation_col = colanno, 
         color = colorRampPalette(c('steelblue','white','brown'))(100), # color = viridis(256, option = "turbo"),
         show_rownames = F, show_colnames = F)
dev.off()

## gene signature
sigres2 <- subset(resLFC_DF, resLFC_DF$pval < .05 & abs(resLFC_DF$log2FoldChange) > logFC_t) # filter for significance
sigres2 <- sigres2 %>% arrange(desc(log2FoldChange))
dim(sigres2)

sigres3 <- sigres2$log2FoldChange # replace estimate with column name for FC values
names(sigres3) <- rownames(sigres2)

## pathway analysis
msigdbr_df_all <- msigdbr(species = "human")

msigdbr_df_hallmarks <- msigdbr::msigdbr(species = "human", category = "H")

# fixing format to work with fgsea
pathwaysH = split(x = msigdbr_df_hallmarks$gene_symbol, f = msigdbr_df_hallmarks$gs_name)
pathwaysH

pathways <- split(x = msigdbr_df_all$gene_symbol, f = msigdbr_df_all$gs_name)

hallmarks <- fgsea(pathways=pathwaysH, sigres3)
all <- fgsea(pathways=pathways, sigres3)

hallmarks
all

topPathwaysUp <- hallmarks[NES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDown <- hallmarks[NES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathwaysH[topPathways], sigres3, hallmarks, 
              gseaParam=0.5)

topPathwaysUpAll <- all[NES > 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysDownAll <- all[NES < 0 & padj < 0.05][head(order(padj), n=10), pathway]
topPathwaysAll <- c(topPathwaysUpAll, rev(topPathwaysDownAll))

pdf("pathway200nM.pdf")
plotGseaTable(pathways[topPathwaysAll], sigres3, all, 
              gseaParam=0.5)
dev.off()

## GO Enrichment
A <- rownames(sigres2)

genelist <- mapIds(org.Hs.eg.db,A,"ENTREZID","SYMBOL")
genelist <- na.omit(genelist)
go <- enrichGO(gene = genelist,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "ALL", # 可选,BP(生物学过程)/CC(细胞组分)/MF(分子功能)/ALL(同时指定)
               pAdjustMethod = "fdr", # P值校正方法,还可以是fdr
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q值阈值
)

go.res <- data.frame(go) # 将GO结果转为数据框，方便后续分析（不转也可以，看个人习惯）
#write.csv(go.res,"Table_GO_result.csv",quote = F) # 输出GO富集分析结果
# 绘制GO富集分析条形图，结果默认按qvalue升序，分别选出前十的term进行绘图即可
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:10,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:10,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:10,]
go.df <- rbind(goBP,goCC,goMF)
# 使画出的GO term的顺序与输入一致
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图
go_bar <- ggplot(data = go.df, # 绘图使用的数据
                 aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
  geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
  coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
  labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
  theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
ggsave(go_bar,filename = "200nM_GO_Enrichment.pdf",width = 9,height = 7)

## Output for future use
sigres2$gene <- rownames(sigres2)
sigres2 <- sigres2[,-7]
sigres200 <- sigres2
save(sigres200,file = "sigres200.RData")

sigres3 <- subset(sigres2,abs(sigres2$log2FoldChange) > logFC_t)
sigres4 <- subset(sigres2,abs(sigres2$log2FoldChange) > 0.75)

sigres3$label <- NULL
sigres3$gene <- rownames(sigres3)

write.csv(sigres3, file = "lfc_t_sigres.csv")

save(sigres3,file = "sigres200_552.RData")
save(sigres4,file = "sigres200_167.RData")


nm100 <- read.delim(file = "CTL vs 100nM CT179-rawdata.txt", row.names = 1, check.names = F)
gse

comb <- cbind(nm100, gse)[-(1:3)]
colnames(comb)

write.csv(comb, file = "raw_counts_GBM8_CT179.csv")

library(singscore)
colnames(sigres3)

up_genes <- sigres3$gene[sigres3$log2FoldChange > 0]
down_genes <- sigres3$gene[sigres3$log2FoldChange < 0]

gene_ranks <- rankGenes(comb)

score <- simpleScore(rankData = gene_ranks, upSet = up_genes, downSet = down_genes, centerScore = T)
score_df <- as.data.frame(score)
score_df$Sample <- rownames(score_df)
score_df$treatment <- str_extract(rownames(score_df), "CT179|ctl")
score_df$label <- str_extract(rownames(score_df), "100nM|200nM|ctl")

score_df

library(dplyr)

score_df$label <- factor(score_df$label, levels = c("200nM", "100nM", "ctl"))

# Compute summary statistics (mean and standard error)
summary_df <- score_df %>%
  group_by(label) %>%
  summarize(meanScore = mean(TotalScore),
            sdScore = sd(TotalScore),
            n = n(),
            se = sdScore / sqrt(n))

# Quick normality test using the Shapiro-Wilk test
shapiro_result <- shapiro.test(score_df$TotalScore)
print(shapiro_result)

# Quick normality test using the Shapiro-Wilk test
shapiro_result <- shapiro.test(summary_df$meanScore)
print(shapiro_result)

# Visual inspection with a histogram and Q-Q plot
# par(mfrow = c(1, 2))
# hist(score_df$TotalScore, 
#      main = "Histogram of TotalScore", 
#      xlab = "TotalScore", 
#      col = "lightblue", 
#      border = "white")
# qqnorm(score_df$TotalScore, 
#        main = "Q-Q Plot of TotalScore")
# qqline(score_df$TotalScore, col = "red")
# par(mfrow = c(1, 1))

# Create the plot with a global mapping for stat_compare_means to work
enrich <- ggplot(score_df, aes(x = label, y = TotalScore, fill = label)) +
  # Plot aggregated bars using summary data
  geom_bar(data = summary_df, 
           aes(x = label, y = meanScore), 
           stat = "identity", width = 0.7, alpha = 0.6) +
  # Add error bars for standard error
  geom_errorbar(data = summary_df, 
                aes(x = label, y = meanScore, ymin = meanScore - se, ymax = meanScore + se),
                width = 0.2, size = 1) +
  # Overlay individual data points with jitter
  geom_jitter(width = 0.2, size = 4, alpha = 0.8, color = "black") +
  # Add statistical comparison between "ctl" and "CT179"
  stat_compare_means(comparisons = list(c("100nM", "200nM"), c("100nM", "ctl"), c("200nM", "ctl")), 
                     method = "t.test", label = "p.format", bracket.size = 1, size = 10) +
  theme_minimal() +
  labs(title = NULL, 
       x = "Treatment", 
       y = "CT-179 Response Enrichment Score") +
  theme(legend.position = "none",
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 18),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        )

pdf(file = "enrichmentByDose.pdf", height = 8, width = 8)
enrich
dev.off()
