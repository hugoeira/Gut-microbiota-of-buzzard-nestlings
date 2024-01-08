# Diferential Abundant Analysis with ANOCOM-BC2



Pipeline adapted from https://bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html



## A) DAA 16s rRNA



### 1. Read in the data

```R
# Load libraries
library(qiime2R)
library(phyloseq)
library(tidyverse)
library(janitor)
library(microbiome)
library(ANCOMBC)
library(ggrepel)

# Import data from qiime2 and create a phyloseq object

ps <- qza_to_phyloseq(
  features="beta-table.qza",
  tree="rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv")


# Edit metadata file

metadata <- sample_data(ps)
metadata <- clean_names(metadata)
taxonomy <- as.data.frame(tax_table(ps))
taxonomy$Kingdom <- gsub("d__","",as.character(taxonomy$Kingdom))
tree <- phy_tree(ps)
asv <- otu_table(ps)

metadata$identifier <- as.factor(metadata$identifier)
metadata$ring_number <- as.factor(metadata$ring_number)
metadata$habitat <- as.factor(metadata$habitat)
metadata$nest <- as.factor(metadata$nest)
metadata$rank <- as.factor(metadata$rank)
metadata$year <- as.factor(metadata$year)
metadata$lbinom <- as.factor(metadata$lbinom)
metadata$sex <- as.factor(metadata$sex)
metadata$age_days <- as.numeric(metadata$age_days)
metadata$std_age <- as.numeric(metadata$std_age)
metadata$bci_two <- as.numeric(metadata$bci_two)
metadata$std_bci_two <- as.numeric(metadata$std_bci_two)
metadata$faith_pd <- as.numeric(metadata$faith_pd)
metadata$shannon_entropy <- as.numeric(metadata$shannon_entropy)

# New phyloseq object
ps <- phyloseq(asv, taxonomy, metadata, tree)

# Save edited ps object as rds
saveRDS(ps,"16S_phyloseq.rds")


```



### 2. Run ANCOM-BC2

#### 2.1 ANCOM-BC2 at ASV level

```R
# Diferential abundace analysis on final model. Genus level
# age + bci + rank + sex + year + habitat + lbinom

ps <- readRDS("16S_phyloseq.rds")

set.seed(123)
output_final = ancombc2(data = ps, assay_name = "counts", tax_level = NULL,
                        fix_formula = " std_age + std_bci_two + rank + sex  + year + habitat+lbinom", # input model
                        rand_formula = "(1|nest/ring_number)", # input random effects
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE, # multicomparison p adjustment
                        prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                        #group = NULL, struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE, # p values cut off
                        #global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                        iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE), # number of iterations
                        em_control = list(tol = 1e-5, max_iter = 100),
                        lme_control = lme4::lmerControl(optimizer ="Nelder_Mead"),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                        nrow = 2, byrow = TRUE), matrix(c(-1, 0, 1, -1), nrow = 2, 
                        byrow = TRUE)), node = list(2, 2), solver = "ECOS", B = 100))
```



##### 2.1.1 ANCOM-BC2 primary analysis

```R
res_prim_asv = output_final_asv$res
```

**Only found asvs that co-vary with age** 



##### 2.1.2 Sensitivity scores

ANCOM-BC2 uses a sensitivity analysis to assess the impact of different pseudo-counts on zero counts for each taxon. The sensitivity score is determined by performing linear regression models on the bias-corrected log abundance table using various pseudo-counts and calculating the proportion of times the p-value exceeds the significance level (alpha). This helps identify taxa that are not sensitive to the pseudo-count addition, ensuring robustness in the analysis.

```
tab_sens_asv = output_asv$pseudo_sens_tab
```



##### 2.2.3 Plot ANCOM-BC2 results

```R
# Volcano plots
volc_asv <- ggplot(data=res_prim_asv, aes(x=lfc_std_age, y=-log10(p_std_age), col=diff_std_age)) + geom_point() + 
  geom_text_repel(aes(label = ifelse(diff_std_age, taxon, ""))) + 
  theme_bw()
volc_asv
```



<img src="/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/volcano-asv-16s.svg" alt="volcano-asv-16s" style="zoom:80%;" />



```R
# Plot log fold changes with unit of age

## Subset a dataset only for age results

df_age_asv = res_prim_asv %>% dplyr::select(taxon, ends_with("age")) # create a dataframe with values only for age

df_fig_age_asv = df_age_asv %>%
  filter(diff_std_age == TRUE) %>% 
  arrange(desc(lfc_std_age)) %>%
  mutate(direct = ifelse(lfc_std_age > 0, "Positive LFC", "Negative LFC")) # prepare the dataframe for plotting

df_fig_age_asv$taxon = factor(df_fig_age_asv$taxon, levels = df_fig_age_asv$taxon)
df_fig_age_asv$direct = factor(df_fig_age_asv$direct, 
                                 levels = c("Positive LFC", "Negative LFC"))

## Make the plot

fig_age_asv = df_fig_age_asv %>%
  ggplot(aes(x = taxon, y = lfc_std_age, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_std_age - se_std_age, ymax = lfc_std_age + se_std_age), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))

fig_age_asv
```



![log-fold-change-age-asv-16s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/log-fold-change-age-asv-16s.svg)

```R
# Pseudo-count sensitivity analysis for age

sens_asv_age = tab_sens_asv %>%
  transmute(taxon, sens_asv_age = std_age) %>%
  left_join(df_age_asv, by = "taxon")

sens_asv_age$diff_std_age = recode(sens_asv_age$diff_std_age * 1, 
                               `1` = "Significant",
                               `0` = "Non-significant")

fig_sens_asv_age = sens_asv_age %>%
  ggplot(aes(x = taxon, y = sens_asv_age, color = diff_std_age)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  labs(x = "ASVs", y = "Sensitivity Score") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())  # Hide the axis text on the x-axis

fig_sens_asv_age

```



![sensitivity-age-asv-16s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/sensitivity-age-asv-16s.svg)

**For the co-variate of age, no outlying sensitivity scores are observed. All significant taxa have low sensitivity scores.**



#### 2.2  ANCOM-BC2  at Genus level

```R
# Diferential abundace analysis on final model. Genus level
# age + bci + rank + sex + year + habitat + lbinom

# Agglomerate at Genus level
ps_genus <- ps %>% tax_glom("Genus") # agglomerate at ge

set.seed(123)

# Run ANCOM BC2

output_genus = ancombc2(data = ps_genus, assay_name = "counts", tax_level = "Genus",
                        fix_formula = " std_age + std_bci_two + rank + sex  + year + habitat + lbinom", # input model
                        rand_formula = "(1|nest/ring_number)", # input random effect
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE, # multicomparison p adjustment if needed
                        prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                        #group = NULL, struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE, # p values cut off
                        #global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                        iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE), # number of iterations
                        em_control = list(tol = 1e-5, max_iter = 100),
                        lme_control = lme4::lmerControl(optimizer ="Nelder_Mead"),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                        nrow = 2, byrow = TRUE), matrix(c(-1, 0, 1, -1), nrow = 2, 
                        byrow = TRUE)), node = list(2, 2), solver = "ECOS", B = 100))

```



##### 2.2.1 ANCOM BC2 primary analysis

```
res_prim_genus = output_genus$res 
```



##### 2.2.2 Sensitivity scores

ANCOM-BC2 uses a sensitivity analysis to assess the impact of different pseudo-counts on zero counts for each taxon. The sensitivity score is determined by performing linear regression models on the bias-corrected log abundance table using various pseudo-counts and calculating the proportion of times the p-value exceeds the significance level (alpha). This helps identify taxa that are not sensitive to the pseudo-count addition, ensuring robustness in the analysis.

```
tab_sens_genus = output_genus$pseudo_sens_tab
```



##### 2.2.3 Plot ANCOM-BC2 results

```R
# Volcano plots

volc_genus <- ggplot(data=res_prim_genus, aes(x=lfc_std_age, y=-log10(p_std_age), col=diff_std_age)) + geom_point() + 
  geom_text_repel(aes(label = ifelse(diff_std_age, taxon, ""))) + theme_bw()
volc_genus

```



```R
# Plot log fold changes with unit of age

df_age_genus = res_prim_genus %>% dplyr::select(taxon, ends_with("age")) # create a dataframe with values only for age

df_fig_age_genus = df_age_genus %>%
  filter(diff_std_age == TRUE) %>% 
  arrange(desc(lfc_std_age)) %>%
  mutate(direct = ifelse(lfc_std_age > 0, "Positive LFC", "Negative LFC")) # prepare the dataframe for plotting

df_fig_age_genus$taxon = factor(df_fig_age_genus$taxon, levels = df_fig_age_genus$taxon)
df_fig_age_genus$direct = factor(df_fig_age_genus$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

fig_age_genus = df_fig_age_genus %>%
  ggplot(aes(x = taxon, y = lfc_std_age, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_std_age - se_std_age, ymax = lfc_std_age + se_std_age), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_age_genus


```



```R
# Pseudo count sensitivity analysis for age

sens_genus_age = tab_sens_asv %>%
  transmute(taxon, sens_asv_age = std_age) %>%
  left_join(df_age_asv, by = "taxon")

sens_genus_age$diff_std_age = recode(sens_asv_age$diff_std_age * 1, 
                               `1` = "Significant",
                               `0` = "Non-significant")

fig_genus_asv_age = sens_asv_age %>%
  ggplot(aes(x = taxon, y = sens_genus_age, color = diff_std_age)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  labs(x = "ASVs", y = "Sensitivity Score") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
theme(axis.text.x = element_blank())  # Hide the axis text on the x-axis

fig_sens_genus_age
```



### 3. Prepare a biom file for functional prediction of the differential abundant genus



#### 3.1 Filter the 9 DAA ASVs
```R

# Create a vector with ASV ID
DAA_asv <- c("51ad0deb6331cb77f17cde5970059539", "eaf0ffd15882c0410c527effbe35f17e", "c287db728494cd6027afb0494a103927", "799f6ae8d556107ff50a68b663e4b663", "af9942270500e33a4c495ffa5dbc42eb", "2f08c5a0c0fecd2166cafb8e15e50202", "92e1de7ddb7eed4bd11841dac5c32fb1", "fb26c90c0b2afc91efe361f689f6b1ca", "6428019f5f4085ec42c0dc69e2cb7d47")

# Filter phyloseq object by ASV ID 
DAA_asv_ps <- prune_taxa(taxa_names(ps) %in% DAA_asv, ps)

# DAA taxa table
DAA_asv_ps@tax_table
Taxonomy Table:     [9 taxa by 7 taxonomic ranks]:
                                 Kingdom    Phylum             Class                 Order               Family              
51ad0deb6331cb77f17cde5970059539 "Bacteria" "Firmicutes"       "Bacilli"             "Lactobacillales"   "Enterococcaceae"   
eaf0ffd15882c0410c527effbe35f17e "Bacteria" "Firmicutes"       "Bacilli"             "Bacillales"        "Bacillaceae"       
c287db728494cd6027afb0494a103927 "Bacteria" "Actinobacteriota" "Actinobacteria"      "Bifidobacteriales" "Bifidobacteriaceae"
799f6ae8d556107ff50a68b663e4b663 "Bacteria" "Actinobacteriota" "Actinobacteria"      "Micrococcales"     "Microbacteriaceae" 
af9942270500e33a4c495ffa5dbc42eb "Bacteria" "Actinobacteriota" "Actinobacteria"      "Actinomycetales"   "Actinomycetaceae"  
2f08c5a0c0fecd2166cafb8e15e50202 "Bacteria" "Actinobacteriota" "Actinobacteria"      "Actinomycetales"   "Actinomycetaceae"  
92e1de7ddb7eed4bd11841dac5c32fb1 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Burkholderiales"   "Comamonadaceae"    
fb26c90c0b2afc91efe361f689f6b1ca "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Xanthomonadales"   "Xanthomonadaceae"  
6428019f5f4085ec42c0dc69e2cb7d47 "Bacteria" "Proteobacteria"   "Gammaproteobacteria" "Xanthomonadales"   "Xanthomonadaceae"  
                                 Genus               Species                        
51ad0deb6331cb77f17cde5970059539 "Enterococcus"      NA                             
eaf0ffd15882c0410c527effbe35f17e NA                  NA                             
c287db728494cd6027afb0494a103927 "Bifidobacterium"   NA                             
799f6ae8d556107ff50a68b663e4b663 "Microbacterium"    NA                             
af9942270500e33a4c495ffa5dbc42eb "Varibaculum"       "uncultured_bacterium"         
2f08c5a0c0fecd2166cafb8e15e50202 "Varibaculum"       "uncultured_bacterium"         
92e1de7ddb7eed4bd11841dac5c32fb1 "Tepidimonas"       NA                             
fb26c90c0b2afc91efe361f689f6b1ca "Pseudoxanthomonas" "Pseudoxanthomonas_taiwanensis"
6428019f5f4085ec42c0dc69e2cb7d47 "Pseudoxanthomonas" "Pseudoxanthomonas_taiwanensis"


# Build new phyloseq object with only DAA genera
DAA_asv_ps <- prune_taxa(taxa_names(ps_genus) %in% DAA_asv, ps_genus)
```




#### 3.2 Convert the filtered phyloseq object to a biom table
```R
library(biomformat)

biom_table <- otu_table(DAA_asv_ps, taxa_are_rows = TRUE)
biom_table<- make_biom(biom_table)

write_biom(biom_table, "DAA-asv.biom")
```



#### 3.3 Filter ASV sequences from a biom table in Qiime2

```python
# Activate qiime
conda activate qiime2-2022.11

#Import biom table into qiime2
qiime tools import --input-path DAA-asv.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path DAA-asv.qza

# Filter DAA representative sequences 
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table DAA-asv.qza --o-filtered-data DAA-seqs.qza

# Deactivate conda environment
conda deactivate

```



### 4. Functional prediction with picrust2

We followed the pipeline here : https://github.com/picrust/picrust2/wiki/Full-pipeline-script

```python
# Activate picrust2
conda activate picrust2

# Run the pipeline
picrust2_pipeline.py -s DAA-seqs.fasta -i DAA-table.biom -o picrust2_out_pipeline -p 1

All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.

All ASVs were below the max NSTI cut-off of 2.0 and so all were retained for downstream analyses.

```



### 5. DAA of functional prediction (ggpicrust2)

Workflow adapted from https://cran.r-project.org/web/packages/ggpicrust2/readme/README.html#workflow

```R
# Load libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(ggprism)
library(patchwork)

# Load KEGG pathway abundance
kegg_abundance_buzzard <- ko2kegg_abundance("pred_metagenome_unstrat.tsv")
saveRDS(kegg_abundance_buzzard,"kegg_abundance_buzzard.rds")

# Perform pathway differential abundance analysis (DAA) using ALDEx2 method
daa_results_df <- pathway_daa(abundance = kegg_abundance_buzzard, metadata = metadata, group = "sampling_point", daa_method = "ALDEx2", p.adjust = "holm", select = NULL, reference = NULL)

# Keep only results p<0.01
daa_sub_sig <- daa_results_df[daa_results_df$p_adjust < 0.01, ]

# Subset by method
daa_sig_results_wilcoxon <- daa_sub_sig[daa_sub_sig$method == "ALDEx2_Wilcoxon rank test", ]
daa_sig_results_wilcoxon <- daa_sig_results_wilcoxon[order(daa_sig_results_wilcoxon$p_adjust), ][1:70, ] # further subset to 70 with lowest p values

daa_sig_results_welch <-daa_results_df[daa_results_df$method == "ALDEx2_Welch's t test", ]
daa_sig_results_welch <- daa_sig_results_welch[order(daa_sig_results_welch$p_adjust), ][1:70, ]


# Annotate pathway results using KO to KEGG conversion
daa_annotated_results_wilcoxon <-pathway_annotation(pathway = "KO", daa_results_df = daa_sig_results_wilcoxon, ko_to_kegg = TRUE)
daa_annotated_results_wilco <-pathway_annotation(pathway = "KO", daa_results_df = daa_sig_results_welch, ko_to_kegg = TRUE)


# Filter annotated results for wilcoxon test
daa_annotated_filtered_wilcoxon <- daa_annotated_results_wilcoxon[!is.na(daa_annotated_results_wilcoxon$pathway_name), ] #remove NA
daa_annotated_filtered_wilcoxon <- daa_annotated_filtered_wilcoxon[!grepl("Human", daa_annotated_filtered_wilcoxon$pathway_class), ]#filter pathway related to humans
daa_annotated_filtered_wilcoxon <- daa_annotated_filtered_wilcoxon[order(daa_annotated_filtered_kruskal$p_adjust), ][1:30, ] # subset to lowest p values

# Filter annotated results for Welch's t test
daa_annotated_filtered_welch <- daa_annotated_results_welch[!is.na(daa_annotated_results_welch$pathway_name), ] #remove NA
daa_annotated_filtered_welch <- daa_annotated_filtered_welch[!grepl("Human", daa_annotated_filtered_welch$pathway_class), ]#filter pathway related to humans
daa_annotated_filtered_welch <- daa_annotated_filtered_welch[order(daa_annotated_filtered_welch$p_adjust), ][1:30, ]


# Generate pathway error bar plot

kegg_abundance_buzzard <- kegg_abundance_buzzard[, -which(names(kegg_abundance_buzzard) == "S213")] # sample 213 has zero kegg abundances
metadata <- metadata %>% filter(sample_id != "S213") # sample 213 has zero kegg abundances 

p_wilcoxon <-  pathway_errorbar(
  abundance = kegg_abundance_buzzard,
  daa_results_df = daa_annotated_filtered_wilcoxon,
  Group = metadata$sampling_point,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "pathway_name")

p_wilcoxon

p2<- pathway_errorbar(
  abundance = kegg_abundance_buzzard,
  daa_results_df = daa_annotated_filtered_welch,
  Group = metadata$sampling_point,
  p_values_threshold = 0.05,
  order = "pathway_class",
  select = NULL,
  ko_to_kegg = TRUE,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "pathway_name")

p2
```



#### 5.1 DAA functional prediction Wilcoxon rank test

<img src="/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/functional-pred-wilcoxon1-16s.svg" alt="functional-pred-wilcoxon1-16s" style="zoom: 200%;" />

#### 5.2. DAA functional prediction Welch t test

<img src="/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/functional-pred-welch-16s.svg" alt="functional-pred-welch-16s" style="zoom: 200%;" />

------



## B) DAA 28s rRNA



### 1. Read in the data

```R
# Import data from qiime2 and create a phyloseq object

ps <- qza_to_phyloseq(
  features="beta-table.qza",
  tree=".rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv")


# Edit metadata file

metadata <- sample_data(ps)
metadata <- clean_names(metadata)
taxonomy <- as.data.frame(tax_table(ps))
taxonomy$Kingdom <- gsub("d__","",as.character(taxonomy$Kingdom))
tree <- phy_tree(ps)
asv <- otu_table(ps)

metadata$identifier <- as.factor(metadata$identifier)
metadata$ring_number <- as.factor(metadata$ring_number)
metadata$habitat <- as.factor(metadata$habitat)
metadata$nest <- as.factor(metadata$nest)
metadata$rank <- as.factor(metadata$rank)
metadata$year <- as.factor(metadata$year)
metadata$lbinom <- as.factor(metadata$lbinom)
metadata$sex <- as.factor(metadata$sex)
metadata$age_days <- as.numeric(metadata$age_days)
metadata$std_age <- as.numeric(metadata$std_age)
metadata$bci_two <- as.numeric(metadata$bci_two)
metadata$std_bci_two <- as.numeric(metadata$std_bci_two)
metadata$faith_pd <- as.numeric(metadata$faith_pd)
metadata$shannon_entropy <- as.numeric(metadata$shannon_entropy)

# New phyloseq object
ps <- phyloseq(asv, taxonomy, metadata, tree)

# Save edited ps object as rds
saveRDS(ps,"28S_phyloseq.rds")


```



### 2. Run ANCOM-BC2



#### 2.1 ANCOM-BC2 at ASV level

```R
# Diferential abundace analysis on final model. Genus level
# age + bci + rank + sex + year + habitat + lbinom

ps <- readRDS("28S_phyloseq.rds")

set.seed(123)
output_final = ancombc2(data = ps, assay_name = "counts", tax_level = NULL,
                        fix_formula = " std_age + std_bci_two + rank + sex  + year + habitat+lbinom", # input model
                        rand_formula = "(1|nest/ring_number)", # input random effects
                        p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE, # multicomparison p adjustment
                        prv_cut = 0.10, lib_cut = 0, s0_perc = 0.05,
                        #group = NULL, struc_zero = TRUE, neg_lb = TRUE,
                        alpha = 0.05, n_cl = 2, verbose = TRUE, # p values cut off
                        #global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                        iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE), # number of iterations
                        em_control = list(tol = 1e-5, max_iter = 100),
                        lme_control = lme4::lmerControl(optimizer ="Nelder_Mead"),
                        mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                        trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                        nrow = 2, byrow = TRUE), matrix(c(-1, 0, 1, -1), nrow = 2, 
                        byrow = TRUE)), node = list(2, 2), solver = "ECOS", B = 100))
```



##### 2.1.1 ANCOM-BC2 primary analysis

```R
res_prim_asv = output_final_asv$res
```

**Only found asvs that co-vary with age** 



##### 2.1.2 Sensitivity scores

```
tab_sens_asv = output_asv$pseudo_sens_tab
```



##### 2.2.3 Plot ANCOM-BC2 results



```R
# Volcano plots
volc_asv <- ggplot(data=res_prim_asv, aes(x=lfc_habitatsouth, y=-log10(p_habitatsouth), col=diff_habitatsouth)) + geom_point() + 
  geom_text_repel(aes(label = ifelse(diff_habitatsouth, taxon, ""))) + 
  theme_bw()
volc_asv

```



![volcano-asv-28s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/volcano-asv-28s.svg)



```R
## Pseudo count sensitivity analysis for habitat_south

df_habitat_asv = res_prim_asv %>% dplyr::select(taxon, ends_with("south")) # create a dataframe with values only for age

sens_asv_habitat = tab_sens_asv %>%
  transmute(taxon, sens_asv_habitat = habitatsouth) %>%
  left_join(df_habitat_asv, by = "taxon")

sens_asv_habitat$diff_habitatsouth = recode(sens_asv_habitat$diff_habitatsouth * 1, 
                                   `1` = "Significant",
                                   `0` = "Non-significant")

fig_sens_asv_habitat = sens_asv_habitat %>%
  ggplot(aes(x = taxon, y = sens_asv_habitat, color = diff_habitatsouth)) +
  geom_point() +
  scale_color_brewer(palette = "Dark2", name = NULL) +
  labs(x = "ASVs", y = "Sensitivity Score") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5))
  theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())  # Hide the axis text on the x-axis

fig_sens_asv_habitat

#For the covariate of habitat no outlying sensitivity scores are observed. All significant taxa have low sensitivity scores.
```



![sensitivity-age-asv-28s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/sensitivity-habitat-asv-28s.svg)

**For the co-variate of habitat (south), no outlying sensitivity scores are observed. All significant taxa have low sensitivity scores.**



```R
# Plot log fold changes. Pairwise comparisson with North

res_pair = output_final_asv2$res_pair

df_habitat = res_prim %>%
  dplyr::select(taxon, contains("habitat")) 


df_fig_habitat = df_habitat %>%
  dplyr::filter((diff_habitatsouth == 1 | diff_habitatteuto == 1)) %>%
  dplyr::mutate(lfc_habitatteuto = ifelse(diff_habitatteuto == 1, 
                                        lfc_habitatteuto, 0),
                lfc_habitatsouth = ifelse(diff_habitatsouth == 1, 
                                  lfc_habitatsouth, 0)) %>%
  dplyr::transmute(taxon, 
                   `Teuto vs. North` = round(lfc_habitatteuto, 2),
                   `South vs. North` = round(lfc_habitatsouth, 2)) %>%
  tidyr::pivot_longer(cols = `Teuto vs. North`:`South vs. North`, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

lo = floor(min(df_fig_habitat$value))
up = ceiling(max(df_fig_habitat$value))
mid = (lo + up)/2
fig_habitat = df_fig_habitat %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to North") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

fig_habitat
```



![log-fold-change-age-asv-28s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/log-fold-change-habitat-asv-28s.svg)



```R
# Plot log fold changes. Multiple pairwise comparisons

res_pair = output_final_asv1$res_pair

df_fig_pair = res_pair %>%
  dplyr::filter((diff_habitatsouth == 1 |diff_habitatteuto == 1 |diff_habitatteuto_habitatsouth == 1)) %>%
  dplyr::mutate(lfc_habitatsouth = ifelse(diff_habitatsouth == 1, 
                                  lfc_habitatsouth, 0),
                lfc_habitatteuto = ifelse(diff_habitatteuto == 1, 
                                        lfc_habitatteuto, 0),
                lfc_habitatteuto_habitatsouth = ifelse(diff_habitatteuto_habitatsouth == 1, 
                                             lfc_habitatteuto_habitatsouth, 0)) %>%
  dplyr::transmute(taxon, 
                   `South vs. North` = round(lfc_habitatsouth, 2),
                   `Teuto vs. North` = round(lfc_habitatteuto, 2),
                   `South vs. Teuto` = round(lfc_habitatteuto_habitatsouth, 2)
  ) %>%
  tidyr::pivot_longer(cols = `South vs. North`:`South vs. Teuto`, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)
df_fig_pair$group = factor(df_fig_pair$group, 
                           levels = c("South vs. North",
                                      "Teuto vs. North",
                                      "South vs. Teuto"))

lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
fig_pair = df_fig_pair %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold change of pairwise comparisons") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_pair
```

![log-fold-pairwie-change-habitat-asv-28s](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/DAA/results/log-fold-pairwise-change-habitat-asv-28s.svg)





### 3. Prepare a biom file for functional prediction of the differential abundant genus



#### 3.1 Filter the 2 DAA ASVs

```R
# Create a vector with ASV ID
DAA_asv <- c("ASV_118", "ASV_352")

# Filter phyloseq object by ASV ID 
DAA_asv_ps <- prune_taxa(taxa_names(ps) %in% DAA_asv, ps)

# DAA taxa table
DAA_asv_ps@tax_table
Taxonomy Table:     [2 taxa by 7 taxonomic ranks]:
        Kingdom        Phylum       Class             Order         Family        Genus         Species               
ASV_118 "d__Eukaryota" "Ascomycota" "Dothideomycetes" "Dothideales" "Dothideales" "Dothideales" "Hormonema_carpetanum"
ASV_352 "d__Eukaryota" "Ascomycota" "Dothideomycetes" "Dothideales" "Dothideales" "Dothideales" "Hormonema_carpetanum"

# Build new phyloseq object with only DAA genera
DAA_asv_ps <- prune_taxa(taxa_names(ps_genus) %in% DAA_asv, ps_genus)
```



#### 3.2 Convert the filtered phyloseq object to a biom table

```R
library(biomformat)

biom_table <- otu_table(DAA_asv_ps, taxa_are_rows = TRUE)
biom_table<- make_biom(biom_table)

write_biom(biom_table, "DAA-asv.biom")
```



#### 3.3 Filter ASV sequences from a biom table in Qiime2

```python
# Activate qiime
conda activate qiime2-2022.11

#Import biom table into qiime2
qiime tools import --input-path DAA-asv.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path DAA-asv.qza

# Filter DAA representative sequences 
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table DAA-asv.qza --o-filtered-data DAA-seqs.qza

# Deactivate conda environment
conda deactivate

```



### 4. Functional prediction with picrust2

Not able to do functional prediction as picrust2 does not have a reference 28S rRNA database.
