
###COMMAND LINE CODE ###############################################################
##Pruning with 0.5 treshold with plink
#cd /users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/mydata
#DATA_DIR=/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/mydata
#VCF=$DATA_DIR/Lottia_gigantea_pl_flt3.renamed.vcf.gz
#module load plink
#plink --version   #version 1.9.0-b7.7

#plink --vcf $VCF --double-id --allow-extra-chr \
#--set-missing-var-ids @:# \
#  --indep-pairwise 50 10 0.5 \
#--out Lg_LD_0.5
#see how many SNPs were kept and removed
#wc -l Lg_LD_0.5.prune.in Lg_LD_0.5.prune.out

  ## how many SNPs did we keep and remove?
  #1821750 total
  #983923 kept
  #837827 removed
#run the PCA
#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
#--extract Lg_LD_0.5.prune.in \
#--make-bed --pca --out Lg_LD_0.5
#see how much variation each PC captures
## first 3 individuals, first 4 PCs (cut -f keeps only the columns we ask for)
#head -n 3 Lg_LD_0.5.eigenvec | cut -d " " -f 1-6
#cat Lg_LD_0.5.eigenval
  #PC1 = 6.24863
  #PC2 = 4.15376
  #PC3 = 2.05948
  #PC4 = 1.55805

##Repeat with 0.1 threshold

#plink --vcf $VCF --double-id --allow-extra-chr \
#--set-missing-var-ids @:# \
#  --indep-pairwise 50 10 0.1 \
#--out Lg_LD_0.1

#see how many SNPs were kept and removed
#wc -l Lg_LD_0.1.prune.in Lg_LD_0.1.prune.out

## how many SNPs did we keep and remove?
# 1821750 total
# 333605 kept
# 1488145 removed

#run the PCA
#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
#--extract Lg_LD_0.1.prune.in \
#--make-bed --pca --out Lg_LD_0.1

#see how much variation each PC captures
## first 3 individuals, first 4 PCs (cut -f keeps only the columns we ask for)
#head -n 3 Lg_LD_0.1.eigenvec | cut -d " " -f 1-6
#cat Lg_LD_0.1.eigenval
  #PC1 = 6.97725
  #PC2 = 3.65953
  #PC3 = 1.8639
  #PC4 = 1.72634

#0.2 Threshold (done in class)

#plink --vcf $VCF --double-id --allow-extra-chr \
#--set-missing-var-ids @:# \
#  --indep-pairwise 50 10 0.2 \
#--out Lg_LD_0.2

## how many SNPs did we keep and remove?
#wc -l Lg_LD_0.2.prune.in Lg_LD_0.2.prune.out
  # 17011 kept
  # 1105 removed
  # 18116 total

#plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# \
#--extract Lg_LD_0.2.prune.in \
#--make-bed --pca --out Lg_LD_0.2

#head -n 3 Lg_LD_0.2.eigenvec | cut -d " " -f 1-6
#cat Lg_LD_0.2.eigenval
  #PC1 = 6.11983
  #PC2 = 5.07366
  #PC3 = 4.48116
  #PC4 = 4.10018


##Calculating Nucleotide Diversity##

#mkdir -p pi ##make subdirectory for pi data
#vcftools --gzvcf $VCF --window-pi 10000 --out pi/all
#head pi/all.windowed.pi
##
3rm pi/all*

  #Calculate pi for each population
  #WINDOW=10000 ###we can change this but keep the code teh same to test multiple window sizes
  
  #for KEEP in pops/*.txt; do
  #POP=$(basename $KEEP .txt)   # e.g. pops/BOD.txt -> BOD
  #echo "== $POP =="
  #vcftools --gzvcf $VCF \
  #--keep $KEEP \
  #--window-pi $WINDOW \
  #--out pi/$POP
  #done

###check the files exist and its how many to expect
  #ls pi
  #ls -l pi/*.pi | wc -l 


#R CODE#################################################################################
library(tidyverse)
setwd ("/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/mydata")
####read in the PCA results from plink
# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.5.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.5.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so drop column 2 (remember: [rows, columns])
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# paste0() joins text and numbers, e.g. paste0("PC", 1:3) gives "PC1" "PC2" "PC3"
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

####see how much of the variation each PC explains
# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# make a folder for our figures

pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.5.png"), width = 6, height = 4)


###PCA PLOT

meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (18.7%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot
OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.5_pca.png"), width = 8, height = 6)

###LATITUDINAL PLOT
lat_plot <- pca %>%
  select(ID, Site, Latitude, PC1, PC2) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.5_pca_latitude.png"), width = 8, height = 6)

###Map of our sites

library(patchwork) ##combined plots

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.5_pca_map.png"), width = 8, height = 6)


##Repeat for 0.1 Threshold##
# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.1.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.1.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so drop column 2 (remember: [rows, columns])
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# paste0() joins text and numbers, e.g. paste0("PC", 1:3) gives "PC1" "PC2" "PC3"
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

####see how much of the variation each PC explains
# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# make a folder for our figures

pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.1.png"), width = 6, height = 4)


###PCA PLOT

meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (18.7%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()
pca_plot
OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.1_pca.png"), width = 8, height = 6)

###LATITUDINAL PLOT
lat_plot <- pca %>%
  select(ID, Site, Latitude, PC1, PC2) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.1_pca_latitude.png"), width = 8, height = 6)

###Map of our sites

library(patchwork) ##combined plots

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults/homework2"

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.1_pca_map.png"), width = 8, height = 6)

###Plotting Pi####
OUTDIR="/users/l/d/ldmathew/projects/eco_genomics_2026/population_genomics/myresults"

meta <- read_csv("sample_metadata.csv")

# one latitude per site, so we can order the sites north to south
sites <- meta %>%
  group_by(Population) %>%
  summarise(Latitude = mean(Latitude))

# the files to read: pi/BOD.windowed.pi, pi/PGR.windowed.pi, ...
pi_files <- list.files("pi", pattern = "windowed.pi$", full.names = TRUE)

pi <- pi_files %>%
  # name each file by its site, e.g. "pi/BOD.windowed.pi" -> "BOD"
  set_names(str_remove(basename(pi_files), fixed(".windowed.pi"))) %>%
  # read every file and stack them, putting the names in a new column, Site
  map_dfr(read_tsv, .id = "Population", show_col_types = FALSE) %>%
  left_join(sites, by = "Population") %>%
  mutate(Site = fct_reorder(Population, Latitude))

# summarise pi for each site grouped by population and latitude
pi_summary <- pi %>%
  group_by(Population, Latitude) %>%
  summarise(n_windows = n(),
            mean_pi   = mean(PI),
            median_pi = median(PI),
            .groups = "drop")%>%
  arrange(desc(mean_pi))

pi_summary
  #PUN, WIL, and SOB sites with greatest mean pi

# pi by lat plot
pi_lat_plot <- ggplot(pi_summary, aes(Latitude, mean_pi, fill = Latitude)) +
  geom_point(size = 4, shape =21, stroke =1.5) +
  labs(x = "Latitude (°N)", y = "Nucleotide diversity") +
  theme_classic()+
  theme(legend.position = "none")+ scale_fill_viridis_c()+
  coord_flip()

ggsave(paste0(OUTDIR, "/pi_vs_latitude.png"), plot = pi_lat_plot, width = 7, height = 7)

