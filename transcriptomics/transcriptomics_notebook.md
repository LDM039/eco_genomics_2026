# Transcriptomics Notebook

**Course**: Intro to Ecological Genomics - Fall 2026

**Name**: Lisa Mathews

------------------------------------------------------------------------

## 9.15.2026 - Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to github

**Working Directory**

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics`

**Input Files**:

`None`

**Output Files**:

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**:

`none`

**Code**:

``` r

print("Hello World")
```

**table**:

#/table

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

**Image**:

#/image

![](images/markdown-syntax-cheatsheet.webp)

**Notes/Observations**:

-   interpretation blahblahblah

**Next steps:**

-   next project!

------------------------------------------------------------------------

## 9.17.2026 - Introduce Study System

-   Learning how to locate .fq files in biol3990 folder through the VACC shell access

**Output Files**:

`None`

**Programs and dependencies**:

-   `VACC Shell access`

**Code:**

prints directory for transcriptomics folder

```{r}
ll | cd CleanData/ 
```

unzips big file and shows first 4 lines (avoids overloading VACC)

```{r}
zcat filename | head -n 4/
```

counts number of lines in the file

```{r}
zcat filename | wc -l
```

**Next steps:**

test for differential expression with fastq data

------------------------------------------------------------------------

## 9.22.2026 - Day 3 of Transcriptomics

Today we set up our R working file to look at A. hutsonica DESeq data

-   Got set up with your Rstudio working environment, repo, data files, and script

-   Continued working in .rmd file to keep your differential gene expression analysis notes together and annotated

-   Imported the counts matrix into DESeq2

-   Visualized reads and variation

-   Visualized global variation in gene expression using Principal Component Analysis (PCA)

**Working Directory**

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/myscripts`

**Input Files**:

`salmon.isoform.counts.matrix.filteredAssembly`

`ahud_samples_R.txt`

**Output Files**:

`/myresults/PCA_allGens.png` `ahud_DESeq2 inclass.R`

**Programs and dependencies**:

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**:

`ahud_DESeq2 inclass.R`

**Code**:

Remove all genes with counts \< 15 in more than 75% of samples (genes w/ too few reads)

```{r}
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
```

Runs DESeq function

```{r}
dds <- DESeq(dds)

Log 2 (n+1) and variance stabalizing transformation graphs

ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

vsd <- vst(dds, blind=FALSE)
meanSdPlot(assay(vsd))
```

Heatmap of sample distance/dissimilarity

```{r}
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$line, vsd$generation, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

Cluster tree that looks for outliers

```{r}
sampleTree <- hclust(dist(sampleDists), method="average")}
plot(sampleTree, main="Sample clustering to detect outliers", sub="", xlab="",cex.lab=1.5, cex.axis=1.5, cex.main=2)
```

Transform the data for plotting using variance stabilization

```{r}
vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("line","generation"), returnData=TRUE) percentVar <- round(100 * attr(pcaData,"percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=line, shape=generation)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed(
```

... Then made lots of PCA plots using ggplot

**PCA Plots:**

![](myresults/PCA_allGens.png)

**Notes/Observations**:

-   PCA:
    -   After 4 generations, all treatments basically synced physiology
    -   After 11 generations, OWA shifted physiology

**Next steps:**

-   Explore the data with more visualizations

-   Plot individual genes

-   Run another model to focus within generation F0 between treatments

-   Make a heat map of the top differentially expressed genes \## 9.24.2026 - Reviewing R coding and Bash Basics

-   Going over bash commands used last week

-   Reviewing A.hudsonica code

**Working Directory**

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/myscripts`

**Input Files**:

`None`

**Output Files**:

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies**:

-   `R version Tidyverse 4.5.1`

-   `R-Studio`

**Scripts**:

`ahud_DESeq2 inclass.R`

**Code**:

**bash** Print working directory `pwd`

Change working directory `cd`

Move back directory `..`

Home directory shortcut `~`

List long (includes file info) `ll`

Print all code entered/changed during session `history`

Copy `cp`

Remove `rm` \*This is permanent!

Unzip file (prints whole file) `zcat`

**ahud working script**

Shows output dimensions

``` r
dim()
```

![](images/Linux-bas-cheatsheet-pg1.webp)

**Notes**:

-   Much of A.hud code included in 9-22 notebook entry
-   All bash commands pertain to terminal coding

**Next steps:**

-   Continue processing + visualizing A.hud data

------------------------------------------------------------------------

## 9.24.2026 - Reviewing R coding and Bash Basics

-   Going over bash commands used last week
-   Reviewing A.hudsonica code

**Working Directory**

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/myscripts`

**Input Files**:

`None`

**Output Files**:

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies**:

-   `R version Tidyverse 4.5.1`

-   `R-Studio`

**Scripts**:

`ahud_DESeq2 inclass.R`

**Code**:

**bash** Print working directory `pwd`

Change working directory `cd`

Move back directory `..`

Home directory shortcut `~`

List long (includes file info) `ll`

Print all code entered/changed during session `history`

Copy `cp`

Remove `rm` \*This is permanent!

Unzip file (prints whole file) `zcat`

**ahud working script**

Shows output dimensions

``` r
dim()
```

![](images/Linux-bas-cheatsheet-pg1.webp)

**Notes**:

-   Much of A.hud code included in 9-22 notebook entry
-   All bash commands pertain to terminal coding

**Next steps:**

-   Continue processing + visualizing A.hud data

------------------------------------------------------------------------

## 9.29.2026 - Day 4: Differential gene expression analysis

-   Analyzed and visualized the counts matrix using a simplied data set (just generation F0)

-   Understanding a contrast and up- versus down-regulation?

-   Learn how to make the various common types of differential gene expression visualizations (ie. Volcano plot, Euler diagram, Heatmap, etc.)

**Working Directory**

`/gpfs1/home/l/d/ldmathew/projects/eco_genomics_2026/transcriptomics/mydata`

**Input Files**:

`salmon.isoform.counts.matrix.filteredAssembly`

`ahud_samples_R.txt`

**Output Files**:

`9.29.26_ahud_DESeqpt2.R`

**Programs and dependencies**:

-   `R version Tidyverse 4.5.1`

-   `R-Studio`

**Scripts**:

`9.26.26_ahud_DESeqpt2.R`

**Code**:

`setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")`

Sets our working directory

`class(volcano_df)`

Tells r what kind of data "volcano_df" holds

`%in%`

This asks, “is this member of that group?”

**Plots**:

Figure 1. Plot of specific top interaction gene TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130

![](myresults/9.29.26_TRINITY_DN30_c0_g2::TRINITY_DN30_c0_g2_i1::g.130::m.130_plot.png){width="445"}

Figure 2. MA plot of log fold change in OW vs AM. Blue indicates significance whereas grey is insignificant.

![](myresults/9.29.26_MAplot.png){width="388"}

Figure 3. Volcano plot of OW vs AM. Red indicates up-regulations while blue indicates down-regulation.

![](myresults/9.29.26_Volcanoplot.png){width="430"}

Figure 4. Heat map of top 20 expressed genes sorted by p-value

![](myresults/9.29.26_heatmap.png){width="437"}

Figure 5. Venn Euler diagram for OWA, OA, and OW.

![](myresults/9.29.26_Euler_plot.png){width="525"}

Figure 6. Upset Plot for OWA, OA, and OW.

![](myresults/9.29.26_Upset_plot.png){width="418"}

**Observations**:

-   Out of significant results:

    -   9.3% are up-regulated in OWA, 14% in OW, and 1.7% in OA

    -   6.2% are down-regulated in OWA, 8.1% in OW, and 0.7% in OA

-   MA plot showed a majority up-regulation of genes in OW vs AM

-   Volcano plot of OW vs AM showed same results as MA plot, but also greater significance in up-regulated genes than down-regulated genes

    -   Positive log fold change is up-regulation, opposite is down-regulation

-   Upset plot shows \# of deferentially expressed genes, sorted by combinations between groups

**Next steps:**

-   DGEA wrap up

-   GO enrichment analysis

-   WGCNA analyses

------------------------------------------------------------------------

## 10.01.2026 - DGEA wrap up and Scatterplotting

-   Scatter plot was used to compare expression responses to OW relative to OWA (each vs. AM control)

-   Color coding plot based on significance

**Working Directory:**

`~/projects/eco_genomics_2026/transcriptomics/mydata`

**Input Files**:

`salmon.isoform.counts.matrix.filteredAssembly`

`ahud_samples_R.txt`

**Output Files**:

`10.01.26_ahud_GO_WGCNA.R`

**Programs and dependencies**:

-   `R version tidyverse 4.5.1`

-   `R-Studio`

**Scripts**:

`10.01.26_ahud_GO_WGCNA.R`

**Code**:

`alpha ()` changes opacity of ggplot dots

`annotate ()` adds text to a figure

*Tidyverse functions*:

`filter()` to remove rows

`mutate()` to add a new variable

`case_when()` to classify genes into categories

`arrange()` to sort the rows

**Images:**

![](myresults/OWA vs AM contrast scatterplot.png){width="455"}

**Notes/Observations**:

-   Log2 Fold Change is the change in gene regulation (up or down) between groups

**Next steps:**

-   Go enrichment and WGCNA analyses
