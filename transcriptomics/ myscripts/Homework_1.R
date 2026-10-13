#### Homework 1: Changes in Gene Expression Between F0 and F4 ####
setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(vsn)  
library(patchwork)


#Import counts matrix 
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)

#Round counts matrix for DESeq2
countsTableRound <- round(countsTable)

#Import sample descriptions table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)

#normalize counts table data
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
dim(dds)
# 67916 genes and 38 samples total

#filter so that there are at least 15 counts per gene in 75% of samples
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]

dim(dds)
## Subset the DESeqDataSet to FO
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# 25260 genes and 12 samples in F0

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

#Check DESeq2 results for F0
resultsNames(dds_F0)

#Name FO OWA vs AM result
F0_res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
#9.3% are up-regulated and 6.2% are down-regulated

#Name F0 OW vs AM
F0_res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
#14% are up-regulated and 8.1% are down-regulated

#Name F0 OA vs AM
F0_res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
#1.7% are up-regulated and 0.7% are down-regulated


## Subset the DESeqDataSet to F4
dds_F4 <- subset(dds, select = generation == 'F4')
dim(dds_F4)
# 25260 genes and 12 samples in F4

# Perform DESeq2 analysis on the subset
dds_F4 <- DESeq(dds_F4)

#Check DESeq2 results for F4
resultsNames(dds_F4)

#Name F4 OWA vs AM result
F4_res_OWAvsAM <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
summary(F4_res_OWAvsAM)
#0.41% are up-regulated and 0.55% are down-regulated

#Name F4 OW vs AM
F4_res_OWvsAM <- results(dds_F4, name="treatment_OW_vs_AM", alpha=0.05)
summary(F4_res_OWvsAM)
#0.44% are up-regulated and 0.17% are down-regulated

#Name F4 OA vs AM
F4_res_OAvsAM <- results(dds_F4, name="treatment_OA_vs_AM", alpha=0.05)
summary(F4_res_OAvsAM)
#0.42% are up-regulated and 0.2% are down-regulated

#plot individual genes to see if the statistics are behaving as we expect for F0
d <-plotCounts(dds_F0, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment")), returnData=TRUE)
d
p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
#plot individual genes to see if the statistics are behaving as we expect for F4
j <-plotCounts(dds_F4, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment")), returnData=TRUE)
j
k <-ggplot(j, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
k <- k + geom_point(position=position_jitter(w=0.2,h=0), size=3)
k <- k + stat_summary(fun = mean, geom = "line")
k <- k + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
k

###PCAplot###
#transform the data for plotting using variance stabilization
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup = c("treatment", "generation"), returnData = TRUE) %>%
  filter(generation %in% c("F0", "F4"))
percentVar <- round(100 * attr(pcaData,"percentVar"))

PCAplot <- ggplot(pcaData, aes(PC1, PC2, color=treatment, shape=generation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

PCAplot

### Euler plot ###
# Comparing DEGs between generations of each treatment

# For OW vs AM
# make a list of significant differentially expressed genes for this contrast

# pull out the results for the contrast of interest
res_OWvsAM_F0 <- results(
  dds_F0, name="treatment_OW_vs_AM", 
  alpha=0.05
)
res_OWvsAM_F0 <- res_OWvsAM_F0[order(res_OWvsAM_F0$padj),] # order them by significance
res_OWvsAM_F0 <- res_OWvsAM_F0[!is.na(res_OWvsAM_F0$padj),] # get rid of any NAs

# make a list of significant differentially expressed genes for this contrast
degs_OWvsAM_F0 <- row.names(res_OWvsAM_F0[res_OWvsAM_F0$padj < 0.05,]) 

#F4
res_OWvsAM_F4 <- results(dds_F4, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM_F4 <- res_OWvsAM_F4[order(res_OWvsAM_F4$padj),]
res_OWvsAM_F4 <- res_OWvsAM_F4[!is.na(res_OWvsAM_F4$padj),]
degs_OWvsAM_F4 <- row.names(res_OWvsAM_F4[res_OWvsAM_F4$padj < 0.05,])

### For OA vs AM
# F0
res_OAvsAM_F0 <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F0 <- res_OAvsAM_F0[order(res_OAvsAM_F0$padj),]
res_OAvsAM_F0 <- res_OAvsAM_F0[!is.na(res_OAvsAM_F0$padj),]
degs_OAvsAM_F0 <- row.names(res_OAvsAM_F0[res_OAvsAM_F0$padj < 0.05,])
#F4
res_OAvsAM_F4 <- results(dds_F4, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM_F4 <- res_OAvsAM_F4[order(res_OAvsAM_F4$padj),]
res_OAvsAM_F4 <- res_OAvsAM_F4[!is.na(res_OAvsAM_F4$padj),]
degs_OAvsAM_F4 <- row.names(res_OAvsAM_F4[res_OAvsAM_F4$padj < 0.05,])

# For OWA vs AM
#F0
res_OWAvsAM_F0 <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F0 <- res_OWAvsAM_F0[order(res_OWAvsAM_F0$padj),]
res_OWAvsAM_F0 <- res_OWAvsAM_F0[!is.na(res_OWAvsAM_F0$padj),]
degs_OWAvsAM_F0 <- row.names(res_OWAvsAM_F0[res_OWAvsAM_F0$padj < 0.05,])
#F4
res_OWAvsAM_F4 <- results(dds_F4, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM_F4 <- res_OWAvsAM_F4[order(res_OWAvsAM_F4$padj),]
res_OWAvsAM_F4 <- res_OWAvsAM_F4[!is.na(res_OWAvsAM_F4$padj),]
degs_OWAvsAM_F4 <- row.names(res_OWAvsAM_F4[res_OWAvsAM_F4$padj < 0.05,])

library(eulerr)

# Total
length(degs_OAvsAM_F0)  # 602
length(degs_OAvsAM_F4)  # 157
length(degs_OWvsAM_F0)  # 5517 
length(degs_OWvsAM_F4)  # 153
length(degs_OWAvsAM_F0)  # 3918
length(degs_OWAvsAM_F4)  # 241

# Intersections
length(intersect(degs_OAvsAM_F0,degs_OWvsAM_F0))  # 444
length(intersect(degs_OAvsAM_F4,degs_OWvsAM_F4))  # 58
length(intersect(degs_OAvsAM_F0,degs_OWAvsAM_F0))  # 380
length(intersect(degs_OAvsAM_F4,degs_OWAvsAM_F4))  # 41
length(intersect(degs_OWAvsAM_F0,degs_OWvsAM_F0))  # 2743
length(intersect(degs_OWAvsAM_F4,degs_OWvsAM_F4))  # 44

# Shared across all
intWA_F0 <- intersect(degs_OAvsAM_F0,degs_OWvsAM_F0)
length(intersect(degs_OWAvsAM_F0,intWA_F0)) # 338
intWA_F4 <- intersect(degs_OAvsAM_F4,degs_OWvsAM_F4)
length(intersect(degs_OWAvsAM_F4,intWA_F4)) # 18

# Number unique to each treatment
# OA F0
602-444-380+338 # 116 
# OW F0
5517-444-2743+338 # 2668
# OWA F0
3918-380-2743+338 # 1133 

# OA F4
157-58-41+18 # 76 
# OW F4
153-58-44+18 # 69
# OWA F4
241-41-44+18 # 174 

# Number shared in pairs of treatments
# OA & OW F0
444-338 # 106 
# OA & OWA F0
380-338 # 42 
# OWA & OW F0
2743-338 # 2405

# OA & OW F4
58-18 # 40
# OA & OWA F4
41-18 # 23
# OWA & OW F4
44-18 # 26

# Making the plot

fit1 <- euler(
  c(
    "OA" = 116, 
    "OW" = 2668, 
    "OWA" = 1133, 
    "OA&OW" = 106, 
    "OA&OWA" = 42, 
    "OW&OWA" = 2405, 
    "OA&OW&OWA" = 338
  )
)

F0_euler <- plot(
  fit1,  
  quantities = TRUE,  
  fill = list(fill = c("gold","orangered", "sienna1", alpha = 0.5)),
  legend = list(labels = c("OA", "OW", "OWA"))
)

fit2 <- euler(
  c(
    "OA" = 76, 
    "OW" = 69, 
    "OWA" = 174, 
    "OA&OW" = 40, 
    "OA&OWA" = 23, 
    "OW&OWA" = 26, 
    "OA&OW&OWA" = 18
  )
)

F4_euler <- plot(
  fit2,  
  quantities = TRUE,  
  fill = list(fill = c("magenta","darkorchid1", "violetred", alpha = 0.5)),
  legend = list(labels = c("OA", "OW", "OWA"))
)

F0_euler
F4_euler

### Volcano plot ###
# Upregulated and downregulated genes for F4 vs. F0 of each treatment

#OW F0

volcano_df <- as.data.frame(res_OWvsAM_F0)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )
#Number of down/up-regulated genes for OW F0
table(volcano_df$sig) 
#Down    NS    Up 
#258 22317  2173 

F0_OW_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "turquoise4",
    "NS"   = "grey70",
    "Up"   = "orangered"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

#OA F0
volcano_df <- as.data.frame(res_OAvsAM_F0)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )
#Number of down/up-regulated genes for OA F0
table(volcano_df$sig) 
#Down    NS    Up 
#99 21391   320 

F0_OA_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "turquoise4",
    "NS"   = "grey70",
    "Up"   = "orangered"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

#OWA F0

volcano_df <- as.data.frame(res_OWAvsAM_F0)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

#Number of down/up-regulated genes for OWA F0
table(volcano_df$sig) 
# Down    NS    Up 
#280 22643  1336 

F0_OWA_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "turquoise4",
    "NS"   = "grey70",
    "Up"   = "orangered"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

#OW F4

volcano_df <- as.data.frame(res_OWvsAM_F4)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

#Number of down/up-regulated genes for OW F4
table(volcano_df$sig) 
#Down    NS    Up 
#18 22731    50

F4_OW_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "palegreen4",
    "NS"   = "grey70",
    "Up"   = "magenta"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

#OA F4

volcano_df <- as.data.frame(res_OAvsAM_F4)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

#Number of down/up-regulated genes for OA F4
table(volcano_df$sig)
# Down    NS    Up 
#24 21739    56 

F4_OA_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "palegreen4",
    "NS"   = "grey70",
    "Up"   = "magenta"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

#OWA F4

volcano_df <- as.data.frame(res_OWAvsAM_F4)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

#Number of down/up-regulated genes for OWA F4
table(volcano_df$sig)
# Down    NS    Up 
#59 24139    70 

F4_OWA_Volcano <- ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "palegreen4",
    "NS"   = "grey70",
    "Up"   = "magenta"
  )) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  theme_classic(base_size = 14) +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = NULL
  )

F0_OW_Volcano
F0_OA_Volcano
F0_OWA_Volcano

F4_OW_Volcano
F4_OA_Volcano
F4_OWA_Volcano


#Combine using patchwork
combined_volcano <-(F0_OW_Volcano | F4_OW_Volcano) /
  (F0_OA_Volcano | F4_OA_Volcano) /
  (F0_OWA_Volcano | F4_OWA_Volcano) +
  plot_annotation(
    title = "Differential Gene Expression",
    subtitle = "F0 vs. F4"
  ) &
  theme(
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

combined_volcano

#Combine volcano plots into one plot
combined_volcano <- ggarrange(
  F0_OW_Volcano, F4_OW_Volcano,
  F0_OA_Volcano, F4_OA_Volcano,
  F0_OWA_Volcano, F4_OWA_Volcano,
  ncol = 2,
  nrow = 3,
  align = "hv",
  common.legend = FALSE
) 

combined_volcano

