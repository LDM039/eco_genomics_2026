print ("hello")
setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")
## Import the libraries that we're likely to need in this session

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample descriptions table (gets metadata of samples)
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

#normalizes data in counts table
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 130580 genes     38 samples- (same as counts table)

# Filter 15 counts per gene in 75% of samples
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# went from 67916 genes to 25260 genes

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# [1] 25260 genes    12 samples

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)
#looks at impact of treatment on only generation 0

### Check on the DE results from the DESeq 

####################################################

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

#Name results so we pull correct names from DESeq
res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
#sort by significance
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
#show most significant results
head(res_OWAvsAM)  
summary(res_OWAvsAM)
#out of significant results, 9.3% are up-regulated in OWA and 6.2% are down-regulated in OWA

#lets do the same for ow vs am
res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),]
head(res_OWvsAM) 
summary(res_OWvsAM)
#out of significant results, 14% are up-regulated in OW and 8.1% are down-regulated in OWA

#aaand the same for oa vs am
res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj),]
head(res_OAvsAM) 
summary(res_OAvsAM)
#out of significant results, 1.7% are up-regulated in OA and 0.7% are down-regulated in OA

### Plot Individual genes ### 

# Counts of specific top interaction gene! (important validation that the normalization, model is working)
d <-plotCounts(dds_F0, gene="TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130", intgroup = (c("treatment")), returnData=TRUE)
d

p <-ggplot(d, aes(x=treatment, y=count, color=treatment)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
#jitter makes sure things dont overlap
p <- p + geom_point(position=position_jitter(w=0.2,h=0), size=3)
p <- p + stat_summary(fun = mean, geom = "line")
p <- p + stat_summary(fun = mean, geom = "point", size=5, alpha=0.7) 
p
#plot shows count of a single specific gene (large dot signifies mean)

#Make an MA plot
plotMA(res_OWvsAM, ylim=c(-5,5))
#shows average number of normalized counts of that gene in the 
#ambient (AM) is 0 line, shows 0-fold change from ambient
#majority up-regulation of genes in OW
#blue=significant, grey=non-significant

######
#Make volcano plot
######
volcano_df <- as.data.frame(res_OWvsAM)

volcano_df <- volcano_df %>%
  mutate(
    sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

ggplot(volcano_df,
       aes(x = log2FoldChange,
           y = -log10(padj),
           color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Down" = "steelblue",
    "NS"   = "grey70",
    "Up"   = "firebrick"
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
#y axis shows the significance factor
#x axis shows change from AM with OW
#blue is down-regulation, red is up-regulation, grey is insignificant results

########################################

# Heatmap of top 20 genes sorted by pvalue

library(pheatmap)

# By environment
vsd <- vst(dds_F0, blind=FALSE)

topgenes <- head(rownames(res_OWvsAM),100)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds_F0)[,c("treatment", "generation")])
pheatmap(mat, annotation_col=df)
pheatmap(mat, annotation_col=df, cluster_cols = F)

# Can you read that? Try this... How/why is it better?
pheatmap(mat, annotation_col=df, cluster_cols = F, show_rownames = F)
#shows genes being expressed similarly

#################################################################

#### PLOT OVERLAPPING DEGS IN VENN EULER DIAGRAM

#################################################################

# For OW vs AM
res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05) # pull out the results for the contrast of interest
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),] # order them by significance
res_OWvsAM <- res_OWvsAM[!is.na(res_OWvsAM$padj),] # get rid of any NAs (high insignificant p vals)
degs_OWvsAM <- row.names(res_OWvsAM[res_OWvsAM$padj < 0.05,]) # make a list of significant differentially expressed genes for this contrast

length(degs_OWvsAM) #5517 significantly different genes
# For OA vs AM
res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj),]
res_OAvsAM <- res_OAvsAM[!is.na(res_OAvsAM$padj),]
degs_OAvsAM <- row.names(res_OAvsAM[res_OAvsAM$padj < 0.05,])

length(degs_OAvsAM) #602 significantly different genes
# For OWA vs AM
res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
res_OWAvsAM <- res_OWAvsAM[!is.na(res_OWAvsAM$padj),]
degs_OWAvsAM <- row.names(res_OWAvsAM[res_OWAvsAM$padj < 0.05,])

length(degs_OWAvsAM) #3918 significantly different genes

library(eulerr)

# Total in each circle
length(degs_OAvsAM)  # 602
length(degs_OWvsAM)  # 5517 
length(degs_OWAvsAM)  # 3918

# Intersections between circles
length(intersect(degs_OAvsAM,degs_OWvsAM))  # 444
length(intersect(degs_OAvsAM,degs_OWAvsAM))  # 380
length(intersect(degs_OWAvsAM,degs_OWvsAM))  # 2743

# Shared across all circles
intWA <- intersect(degs_OAvsAM,degs_OWvsAM)
length(intersect(degs_OWAvsAM,intWA)) # 338

# Number unique to each treatment

602-444-380+338 # 116 OA
5517-444-2743+338 # 2668 OW 
3918-380-2743+338 # 1133 OWA

# Number shared in pairs of treatments

444-338 # 106 OA & OW
380-338 # 42 OA & OWA
2743-338 # 2405 OWA & OW

# Now assemble the results
# Note that the names are important and have to be specific to line up the diagram
fit1 <- euler(c("OA" = 116, "OW" = 2668, "OWA" = 1133, "OA&OW" = 106, "OA&OWA" = 42, "OW&OWA" = 2405, "OA&OW&OWA" = 338))

# And make the plot!
plot(fit1,  lty = 1:3, quantities = TRUE)
# lty changes the lines

plot(fit1, quantities = TRUE, fill = "transparent",
     lty = 1:3,
     labels = list(font = 4))


#cross check with above lengths of DEGS: the four values, unique, shared with one other, shared with the second other, shared across all treatments, should sum to the length of DEGs for each treatment contrast to AM
2668+2405+338+106 # 5517 total OW
1133+2405+338+42  # 3918 total OWA
116+42+106+338    # 602  total OA

#######
#Lets make an upset plot
######

install.packages("UpSetR")
library(UpSetR)

#organize numbers into format for upset plot
all_genes <- unique(c(
  degs_OAvsAM,
  degs_OWvsAM,
  degs_OWAvsAM
))

upset_df <- data.frame(
  gene = all_genes,
  OA = all_genes %in% degs_OAvsAM,
  OW = all_genes %in% degs_OWvsAM,
  OWA = all_genes %in% degs_OWAvsAM
)

head(upset_df)
#True false values assigns membership for one gene to other genes

deg.list <- list(
  OA  = degs_OAvsAM,
  OW  = degs_OWvsAM,
  OWA = degs_OWAvsAM
)
#Pull differential expressed genes, order based on T or F, label them
upset(
  fromList(deg.list),
  order.by = "freq",
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)

####################### A bit prettier
data.
upset(
  fromList(deg.list),
  order.by = "freq",
  main.bar.color = "grey30",
  sets.bar.color = c("#00A08A", "#CC3333", "#F2AD00"), # had to manually adjust the order
  mainbar.y.label = "Number of DEGs",
  sets.x.label = "Total DEGs"
)

#This shows same info from Euler plot but bars show # of differentially expressed genes, sorted by combinations between groups
#OW and OWA have second biggest effect