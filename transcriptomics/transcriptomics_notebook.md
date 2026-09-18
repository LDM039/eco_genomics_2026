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
