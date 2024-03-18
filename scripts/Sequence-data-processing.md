# Sequence data processing

- [A) 16s rRNA sequence data processing](#a--16s-rrna-sequence-data-processing)
- [1.  Activate Qiime2](#1--activate-qiime2)
- [2. Import sequences](#2-import-sequences)
- [3. Visualize quality plots](#3-visualize-quality-plots)
  * [3.1. Summary of raw reads](#31-summary-of-raw-reads)
  * [3.2. Quality plots](#32-quality-plots)
- [4. Run dada2](#4-run-dada2)
  * [4.1. DADA2 results](#41-dada2-results)
    + [4.1.1. Table summary](#411-table-summary)
    + [4.1.2. Summary reads per samples](#412-summary-reads-per-samples)
    + [4.1.2. Sequence Length Statistics](#412-sequence-length-statistics)
    + [4.1.3. Denoising statistics](#413-denoising-statistics)
- [5. Taxonomy assignment](#5-taxonomy-assignment)
  * [5.1. Taxonomy visualisation](#51-taxonomy-visualisation)
- [6. Exit qiime2](#6-exit-qiime2)
- [7. In R run decontam](#7-in-r-run-decontam)
  * [7.1. Read in the data](#71-read-in-the-data)
  * [7.2. Run decontam](#72-run-decontam)
  * [7.3. Export feature table as biom file](#73-export-feature-table-as-biom-file)
- [8. Import biom table from R to qiime2](#8-import-biom-table-from-r-to-qiime2)
- [9. Remove control samples from dataset](#9-remove-control-samples-from-dataset)
- [10. Taxonomy based filtering](#10-taxonomy-based-filtering)
- [11. Filter unique features](#11-filter-unique-features)
- [12. Filter samples with less than 500 reads](#12-filter-samples-with-less-than-500-reads)
- [13. Filter representative sequences](#13-filter-representative-sequences)
- [14. Building a phylogenetic tree](#14-building-a-phylogenetic-tree)
- [15. Rarefaction curves](#15-rarefaction-curves)
- [16. Taxa bar plots](#16-taxa-bar-plots)
  * [16.1. Plot core phylum (present in 70% of samples)](#161-plot-core-phylum--present-in-70--of-samples-)
  * [16.2. Plot core families (present in 70% of samples)](#162-plot-core-families--present-in-70--of-samples-)
- [17. Calculate alpha diversity metrics and rarefy the data-set](#17-calculate-alpha-diversity-metrics-and-rarefy-the-data-set)
- [18. Filter samples with less than 4000 reads (unrarefied table for beta analysis)](#18-filter-samples-with-less-than-4000-reads--unrarefied-table-for-beta-analysis-)
  * [18.1. Table Summary](#181-table-summary)
  * [18.2. Summary of final reads per sample](#182-summary-of-final-reads-per-sample)
  * [18.3. Sequence length statistics](#183-sequence-length-statistics)
- [B) 28s rRNA sequence data processing](#b--28s-rrna-sequence-data-processing)
- [1.  Activate Qiime2](#1--activate-qiime2-1)
- [2. Import sequences](#2-import-sequences-1)
- [3. Visualize quality plots](#3-visualize-quality-plots-1)
  * [3.1 Summary of raw reads](#31-summary-of-raw-reads)
  * [3.2 Quality plots](#32-quality-plots)
- [4. Exit qiime](#4-exit-qiime)
- [5. Run R pipeline to concatenate reads](#5-run-r-pipeline-to-concatenate-reads)
- [6. Import data from R to qiime2](#6-import-data-from-r-to-qiime2)
  * [6.1. Visualising results from R pipeline](#61-visualising-results-from-r-pipeline)
    + [6.1.1 Summary of primer trimmed reads](#611-summary-of-primer-trimmed-reads)
    + [6.1.2. Primer trimmed quality plots](#612-primer-trimmed-quality-plots)
    + [6.1.3 Denoising statistics](#613-denoising-statistics)
    + [6.1.4. Table summary (after denoising)](#614-table-summary--after-denoising-)
    + [6.1.5 Summary reads per sample](#615-summary-reads-per-sample)
    + [6.1.6. Sequence Length Statistics](#616-sequence-length-statistics)
- [7. Taxonomy assignment](#7-taxonomy-assignment)
- [8. Exit qiime](#8-exit-qiime)
- [8. In R run decontam](#8-in-r-run-decontam)
  * [8.1. Read the data](#81-read-the-data)
  * [8.2. Run decontam](#82-run-decontam)
  * [8.3. Export feature table as biom file](#83-export-feature-table-as-biom-file)
- [9. Import biom table from R to qiime2](#9-import-biom-table-from-r-to-qiime2)
- [10. Remove control samples from dataset](#10-remove-control-samples-from-dataset)
- [16. Taxonomy based filtering](#16-taxonomy-based-filtering)
- [17. Filter unique features](#17-filter-unique-features)
- [18. Filter samples with less than 500 reads](#18-filter-samples-with-less-than-500-reads)
- [19. Filter representative sequences](#19-filter-representative-sequences)
- [20. Build a phylogenetic tree](#20-build-a-phylogenetic-tree)
- [21. Rarefaction curves](#21-rarefaction-curves)
    + [21.1 Rarefaction plots (sampling depth 20,000)](#211-rarefaction-plots--sampling-depth-20-000-)
    + [21.1 Rarefaction plots (sampling depth 2000)](#211-rarefaction-plots--sampling-depth-2000-)
  * [22. Taxa bar plots](#22-taxa-bar-plots)
  * [22.1. Plot core phylum (present in 70% of samples)](#221-plot-core-phylum--present-in-70--of-samples-)
  * [22.2. Plot core families (present in 70% of samples)](#222-plot-core-families--present-in-70--of-samples-)
- [23. Calculate alpha diversity metrics and rarefy the data-set](#23-calculate-alpha-diversity-metrics-and-rarefy-the-data-set)
- [24. Filter samples with less than 2000 reads (unrarefied table for beta analysis)](#24-filter-samples-with-less-than-2000-reads--unrarefied-table-for-beta-analysis-)
  * [24.1. Table Summary](#241-table-summary)
  * [24.2. Summary of final reads per sample](#242-summary-of-final-reads-per-sample)
  * [24.2. Sequence length statistics](#242-sequence-length-statistics)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>



## A) 16s rRNA sequence data processing



## 1.  Activate Qiime2

```bash
# Activate base env
. ~/.bashrc

# Activate qiime2
conda activate qiime2-2022.11
```



## 2. Import sequences

```python
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path seqs --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza 
```



## 3. Visualize quality plots

```python
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-quality-plots.qzv 
```



### 3.1. Summary of raw reads

|             | Forward reads | Reverse reads |
| ----------- | :-----------: | :-----------: |
| **Minimum** |      77       |      77       |
| **Median**  |     59119     |     59119     |
| **Mean**    |  57956.0754   |  57956.0754   |
| **Maximum** |     96300     |     96300     |
| **Total**   |   14604931    |   14604931    |



### 3.2. Quality plots

![16s-quality-plots](/pics/16s-quality-plots.svg)



## 4. Run dada2 

```python
qiime dada2 denoise-paired --i-demultiplexed-seqs demux-paired-end.qza --p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 253 --p-trunc-len-r 185  --p-trunc-q 2 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza
```



### 4.1. DADA2 results

```python
qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv | 

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv | 

qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file buzzard_meta.tsv
```



#### 4.1.1. Table summary

|                               |  Sample   |
| ----------------------------- | :-------: |
| **Number of samples**         |    252    |
| **Number of features (ASVs)** |  10,065   |
| **Total frequency**           | 6,319,664 |



#### 4.1.2. Summary reads per samples

|                       | Frequency |
| --------------------- | :-------: |
| **Minimum frequency** |     0     |
| **1st quartile**      | 19,110.25 |
| **Median frequency**  | 25,091.50 |
| **3rd quartile**      | 31,849.25 |
| **Maximum frequency** | 50,886.00 |
| **Mean frequency**    | 25,078.03 |



#### 4.1.2. Sequence Length Statistics

| Sequence count | Min Lenght | Max Lenght | Mean Lenght | Range | Standard Deviation |
| :------------: | :--------: | :--------: | :---------: | :---: | :----------------: |
|     10065      |    233     |    383     |   254.814   |  150  |       17.408       |



#### 4.1.3. Denoising statistics
Not available here.


## 5. Taxonomy assignment

```python
qiime feature-classifier classify-sklearn --i-classifier silva-138.1-SSU-nr99-515F-806R-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza 
```



### 5.1. Taxonomy visualisation

```python
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
```



## 6. Exit qiime2

```python
conda deactivate 
```



## 7. In R run decontam

### 7.1. Read in the data

```R
# Load Libraries
library(decontam)
library(qiime2R)
library(phyloseq)
library(biomformat)

# Make a phyloseq object
ps <- qza_to_phyloseq(features = "table.qza", taxonomy = "taxonomy.qza", metadata = "buzzard_meta.tsv")

# Choose which samples are the negative controls 
sample_data(ps)$is.neg <- sample_data(ps)$type == "negative"
```



### 7.2. Run decontam 

```R
# Identify contaminants based on prevalence method (treshold 0.1 is the standard)
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Remove contaminants from the phyloseq object
ps.nocontam <- prune_taxa(!contamdf.prev$contaminant,ps)
```

decontam results



### 7.3. Export feature table as biom file

```R
# Extract asv table from the phyloseq object
table_nocontam <- as(otu_table(ps.nocontam),"matrix",) 

#'t' to transform if taxa_are_rows=FALSE table_nocontam<- t(as(otu_table(ps.nocontam),"matrix",))#if taxa_are_rows=TRUE

# Make a biom table
table_nocontam_biom <- make_biom(data=table_nocontam)
write_biom(table_nocontam_biom,"table-nocontam.biom")
```



## 8. Import biom table from R to qiime2

```python
conda activate qiime2-2022.11

qiime tools import --input-path table-nocontam.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path table-nocontam.qza
```



## 9. Remove control samples from dataset
```python
qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file buzzard_meta.tsv --p-where "type='positive'" --p-exclude-ids --o-filtered-table table-nocontam.qza

qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file buzzard_meta.tsv --p-where "type='negative'" --p-exclude-ids --o-filtered-table table-nocontam.qza
```



## 10. Taxonomy based filtering

Filter out mitochondrial, chloroplast, unassigned, Vertebrata, Eukaryote and taxa not assigned to phylum

```python
qiime taxa filter-table --i-table table-nocontam.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast,Unassigned,Vertebrata,Eukaryota --p-include p_ --o-filtered-table table-taxa-filter.qza

qiime feature-table summarize --i-table table-taxa-filter.qza --o-visualization table-taxa-filter.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 11. Filter unique features

```python
qiime feature-table filter-features --i-table table-taxa-filter.qza --p-min-samples 2 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 12. Filter samples with less than 500 reads

```python
qiime feature-table filter-samples --i-table table-taxa-filter-no_singles.qza --p-min-frequency 500 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 13. Filter representative sequences

```python
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table filtered-table.qza --o-filtered-data filter-seqs.qza
qiime feature-table tabulate-seqs --i-data filter-seqs.qza --o-visualization filter-seqs.qzv 
```



## 14. Building a phylogenetic tree

```python
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences filter-seqs.qza  --o-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-seqs.qza --o-tree unrooted-tree.qza  --o-rooted-tree rooted-tree.qza 
```



## 15. Rarefaction curves

```python
qiime diversity alpha-rarefaction --i-table filtered-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 20000 --m-metadata-file buzzard_meta.tsv  --o-visualization alpha-rarefaction.qzv

qiime diversity alpha-rarefaction --i-table filtered-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 4000 --m-metadata-file buzzard_meta.tsv  --o-visualization alpha-rarefaction.qzv
```



![16s-rarefaction-plots](/pics/16s-rarefaction-plots.svg)



## 16. Taxa bar plots

```
qiime taxa barplot --i-table filtered-table.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization taxa-bar-plots.qzv 
```



### 16.1. Plot core phylum (present in 70% of samples)

```python
qiime taxa collapse --i-table filtered-table.qza --i-taxonomy taxonomy.qza --o-collapsed-table phylum-table.qza  --p-level 2

qiime feature-table summarize --i-table phylum-table.qza --o-visualization phylum-table.qzv --m-sample-metadata-file buzzard_meta.tsv

qiime feature-table core-features --i-table phylum-table.qza --o-visualization phylum-core-table.qzv 

```

From the visualization file download the tsv relative to features present in 70% of the samples (or any other % that I want)

```python
qiime metadata tabulate --m-input-file core-features-0.700.tsv --o-visualization core-features-0.700.qzv

qiime feature-table filter-features --i-table filtered-table.qza --o-filtered-table table-phylum-core-0.700.qza --m-metadata-file core-features-0.700.tsv

qiime taxa barplot --i-table table-phylum-core-0.700.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization core-phylum-bar-plots.qzv
```



| Feature ID                   | 2%      | 9%      | 25%     | 50%    | 75%     | 91%      | 98%      |
| ---------------------------- | ------- | ------- | ------- | ------ | ------- | -------- | -------- |
| dBacteria;p_Firmicutes       | 2184.9  | 3847.59 | 6083.25 | 8120.5 | 11548.5 | 14100.62 | 17146.1  |
| dBacteria;p_Actinobacteriota | 1985.22 | 3404.05 | 4962.5  | 7219   | 9077.5  | 10892.29 | 13493.08 |
| dBacteria;p_Proteobacteria   | 757.18  | 1443.83 | 2297.75 | 3928.5 | 5697.5  | 7565.78  | 10199.84 |
| dBacteria;p_Bacteroidota     | 0       | 15.61   | 79      | 260    | 674     | 2138.07  | 8270.38  |
| dBacteria;p_Campylobacterota | 0       | 0       | 17.5    | 177.5  | 526.25  | 1400.73  | 3014.1   |



### 16.2. Plot core families (present in 70% of samples)

```python
qiime taxa collapse --i-table filtered-table.qza --i-taxonomy taxonomy.qza --o-collapsed-table family-table.qza  --p-level 5

qiime feature-table summarize --i-table -table.qza --o-visualization family-table.qzv --m-sample-metadata-file buzzard_meta.tsv

qiime feature-table core-features --i-table family-table.qza --o-visualization family-core-table.qzv 
```

From the visualization file download the tsv relative to features present in 70% of the samples (or any other % that I want)

```python
qiime metadata tabulate --m-input-file core-features-0.700.tsv --o-visualization core-features-0.700.qzv

qiime feature-table filter-features --i-table filtered-table.qza --o-filtered-table table-family-core-0.700.qza --m-metadata-file core-features-0.700.tsv

qiime taxa barplot --i-table table-family-core-0.700.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization core-family-bar-plots.qzv
```



| Feature ID                                                   | 2%     | 9%     | 25%     | 50%    | 75%     | 91%     | 98%     |
| ------------------------------------------------------------ | ------ | ------ | ------- | ------ | ------- | ------- | ------- |
| dBacteria;pActinobacteriota;cActinobacteria;oCorynebacteriales;f__Corynebacteriaceae | 350.74 | 872.57 | 1722.75 | 3061   | 4669.25 | 6330.39 | 8077.36 |
| dBacteria;pFirmicutes;cClostridia;oPeptostreptococcales-Tissierellales;f__Peptostreptococcaceae | 44.64  | 221.27 | 977.5   | 2371.5 | 3711.5  | 5851.84 | 8581.56 |
| dBacteria;pActinobacteriota;cActinobacteria;oActinomycetales;f__Actinomycetaceae | 0      | 9.61   | 128.5   | 1542   | 3403.25 | 5059.41 | 6521.82 |
| dBacteria;pProteobacteria;cGammaproteobacteria;oEnterobacterales;f__Enterobacteriaceae | 17.16  | 85.37  | 290.5   | 747    | 1892    | 3955.99 | 7019.76 |
| dBacteria;pActinobacteriota;cActinobacteria;oPropionibacteriales;f__Propionibacteriaceae | 0      | 39.61  | 116.75  | 344.5  | 803     | 1854.51 | 3679.56 |
| dBacteria;pFirmicutes;cBacilli;oStaphylococcales;f__Gemellaceae | 0      | 0      | 49.25   | 302    | 1325    | 3811.18 | 7817.5  |
| dBacteria;pFirmicutes;cBacilli;oStaphylococcales;f__Staphylococcaceae | 0      | 25.61  | 92.25   | 246.5  | 596.25  | 1271.9  | 2761.72 |
| dBacteria;pFirmicutes;cBacilli;oLactobacillales;f__Lactobacillaceae | 0      | 18     | 67      | 232    | 684     | 1414.26 | 3788.36 |
| dBacteria;pFirmicutes;cBacilli;oLactobacillales;f__Enterococcaceae | 0      | 0      | 64.5    | 208    | 569     | 1411.95 | 3548.1  |
| dBacteria;pProteobacteria;cGammaproteobacteria;oXanthomonadales;f__Xanthomonadaceae | 0      | 0      | 49.25   | 199    | 445     | 821.73  | 1358.34 |
| dBacteria;pFirmicutes;cBacilli;oBacillales;f__Bacillaceae    | 0      | 0      | 50.75   | 182    | 399.75  | 837.68  | 1659.94 |
| dBacteria;pCampylobacterota;cCampylobacteria;oCampylobacterales;f__Campylobacteraceae | 0      | 0      | 12.25   | 173.5  | 505     | 1294.68 | 2993.46 |
| dBacteria;pProteobacteria;cGammaproteobacteria;oBurkholderiales;f__Comamonadaceae | 0      | 13.61  | 61.25   | 171.5  | 411     | 827     | 1608.6  |
| dBacteria;pFirmicutes;cClostridia;oClostridiales;f__Clostridiaceae | 0      | 0      | 6       | 87.5   | 482     | 1778.88 | 3636.3  |
| dBacteria;pProteobacteria;cGammaproteobacteria;oPseudomonadales;f__Moraxellaceae | 0      | 0      | 16.5    | 75     | 272.5   | 758.23  | 2257.12 |
| dBacteria;pActinobacteriota;cActinobacteria;oBifidobacteriales;f__Bifidobacteriaceae | 0      | 0      | 8.25    | 74.5   | 185.5   | 372.34  | 1114.06 |
| dBacteria;pProteobacteria;cAlphaproteobacteria;oRhizobiales;f__Rhizobiaceae | 0      | 0      | 14      | 50.5   | 174     | 424.24  | 827.8   |
| dBacteria;pBacteroidota;cBacteroidia;oFlavobacteriales;f__Weeksellaceae | 0      | 0      | 5       | 31.5   | 117.75  | 467.9   | 1931.38 |
| dBacteria;pProteobacteria;cAlphaproteobacteria;oSphingomonadales;f__Sphingomonadaceae | 0      | 0      | 4.25    | 27     | 79.75   | 189.07  | 461.82  |



## 17. Calculate alpha diversity metrics and rarefy the data-set

```python
# Calculates observed features, shannon diversity, Faith PD and retrieves rarefied table
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table filtered-table.qza --p-sampling-depth 4000 --m-metadata-file buzzard_meta.tsv --output-dir alpha-metrics-results

#Add alpha diversity metrics to the metadata
qiime metadata tabulate --m-input-file buzzard_metadata.tsv --m-input-file shannon_vector.qza --m-input-file observed_features.qza --m-input-file faith_pd_vector.qza --o-visualization buzzard_meta_alpha.qzv


```

From the "buzzard_meta_alpha.qzv" one can extract the tsv file with the metadata + alpha metrics.



## 18. Filter samples with less than 4000 reads (unrarefied table for beta analysis)

```
qiime feature-table filter-samples --i-table filtered-table.qza --p-min-frequency 4000 --o-filtered-table beta-table.qza

qiime feature-table summarize --i-table beta-table.qza --o-visualization beta-table.qzv --m-sample-metadata-file buzzard_meta_alpha.tsv
```



### 18.1. Table Summary

|                               |  Sample   |
| ----------------------------- | :-------: |
| **Number of samples**         |    230    |
| **Number of features (ASVs)** |   2,078   |
| **Total frequency**           | 5,121,868 |



### 18.2. Summary of final reads per sample

|                       | Frequency |
| --------------------- | :-------: |
| **Minimum frequency** | 4,295.00  |
| **1st quartile**      | 16,586.25 |
| **Median frequency**  | 22,162.00 |
| **3rd quartile**      | 27,903.50 |
| **Maximum frequency** | 43,074.00 |
| **Mean frequency**    | 22,268.99 |



### 18.3. Sequence length statistics

| Sequence Count | Min Length | Max Length | Mean Length | Range | Standard Deviation |
| :------------: | :--------: | :--------: | :---------: | :---: | ------------------ |
|      2078      |    247     |    291     |   251.96    |  44   | 1.34               |

------





## B) 28s rRNA sequence data processing



## 1.  Activate Qiime2

```bash
# Activate base env
. ~/.bashrc

# Activate qiime2
conda activate qiime2-2022.11
```



## 2. Import sequences

```python
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path seqs --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza
```



## 3. Visualize quality plots

```python
qiime demux summarize --i-data demux-paired-end.qza --o-visualization demux-quality-plots.qzv
```



### 3.1 Summary of raw reads

|             | Forward reads | Reverse reads |
| ----------- | ------------- | ------------- |
| **Minimum** | 311           | 311           |
| **Median**  | 55034.5       | 55034.5       |
| **Mean**    | 83582.93254   | 83582.93254   |
| **Maximum** | 475700        | 475700        |
| **Total**   | 21062899      | 21062899      |



### 3.2 Quality plots

![28s-quality-plots](/pics/28s-quality-plots.svg)



## 4. Exit qiime

```python
conda deactivate
```



## 5. Run R pipeline to concatenate reads

```bash
Rscript 28s-concatenate-reads.R
```



## 6. Import data from R to qiime2

```python
#Activate qiime2
conda activate qiime2-2022.11

# Import biom table (final result from DADA2 pipeline)
qiime tools import --input-path table.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path table.qza

# Import ASV sequences (fasta file)
qiime tools import --input-path ASV_nochim.fasta --output-path rep-seqs.qza --type 'FeatureData[Sequence]'
```



### 6.1. Visualising results from R pipeline

```python
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file buzzard_meta.tsv
```



#### 6.1.1 Summary of primer trimmed reads

|             | Forward reads | Reverse reads |
| ----------- | :-----------: | :-----------: |
| **Minimum** |      173      |      173      |
| **Median**  |    30195.5    |    30195.5    |
| **Mean**    |  34482.53571  |  34482.53571  |
| **Maximum** |    251247     |    251247     |
| **Total**   |    8689599    |    8689599    |

#### 6.1.2. Primer trimmed quality plots

![28s-primer-trimmed-quality-plots](/pics/28s-primer-trimmed-quality-plots.svg)



#### 6.1.3 Denoising statistics

Not available here.


#### 6.1.4. Table summary (after denoising)

|                               |  Sample   |
| ----------------------------- | :-------: |
| **Number of samples**         |    252    |
| **Number of features (ASVs)** |  14,143   |
| **Total frequency**           | 2,710,821 |



#### 6.1.5 Summary reads per sample

|                       | Frequency |
| --------------------- | :-------: |
| **Minimum frequency** |     0     |
| **1st quartile**      | 6,773.50  |
| **Median frequency**  | 9,894.00  |
| **3rd quartile**      | 15,024.00 |
| **Maximum frequency** | 33,905.00 |
| **Mean frequency**    | 10,757.23 |
|                       |           |

#### 6.1.6. Sequence Length Statistics

| Sequence Count | Min Length | Max Length | Mean Length | Range | Standard Deviation |
| :------------: | :--------: | :--------: | :---------: | :---: | :----------------: |
|     14143      |    230     |    415     |   410.82    |  185  |        20.7        |



## 7. Taxonomy assignment

```python
qiime feature-classifier classify-sklearn --i-classifier silva-138.1-LSU-nr99-GA20F-RM9R-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

# Taxonomy visualisation
qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv
```



## 8. Exit qiime

```bash
conda deactivate
```



## 8. In R run decontam

### 8.1. Read the data

```R
# Load Libraries
library(decontam)
library(qiime2R)
library(phyloseq)
library(biomformat)

# Make a phyloseq object
ps <- qza_to_phyloseq(features = "table.qza", taxonomy = "taxonomy.qza", metadata = "buzzard_meta.tsv")

# Choose which samples are the negative controls 
sample_data(ps)$is.neg <- sample_data(ps)$type == "negative"
```



### 8.2. Run decontam 

```R
# Identify contaminants based on prevalence method (treshold 0.1 is the standard)
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1)

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Remove contaminants from the phyloseq object
ps.nocontam <- prune_taxa(!contamdf.prev$contaminant,ps)
```



### 8.3. Export feature table as biom file

```R
# Extract asv table from the phyloseq object
table_nocontam <- as(otu_table(ps.nocontam),"matrix",) 

#'t' to transform if taxa_are_rows=FALSE table_nocontam<- t(as(otu_table(ps.nocontam),"matrix",))#if taxa_are_rows=TRUE

# Make a biom table
table_nocontam_biom <- make_biom(data=table_nocontam)
write_biom(table_nocontam_biom,"table-nocontam.biom")
```



## 9. Import biom table from R to qiime2

```python
conda activate qiime2-2022.11

qiime tools import --input-path table-nocontam.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path table-nocontam.qza
```



## 10. Remove control samples from dataset

```python
qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file buzzard_meta.tsv --p-where "type='positive'" --p-exclude-ids --o-filtered-table table-nocontam.qza

qiime feature-table filter-samples --i-table table-nocontam.qza --m-metadata-file buzzard_meta.tsv --p-where "type='negative'" --p-exclude-ids --o-filtered-table table-nocontam.qza
```



## 16. Taxonomy based filtering

Filter out mitochondrial, Chloroplast, Unassigned, Vertebrata, Bacteria and taxa not assigned to phylum

```python
qiime taxa filter-table --i-table table-nocontam.qza --i-taxonomy taxonomy.qza --p-exclude mitochondria,chloroplast,Unassigned,Vertebrata --p-include p_ --o-filtered-table table-taxa-filter.qza

qiime feature-table summarize --i-table table-taxa-filter.qza --o-visualization table-taxa-filter.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 17. Filter unique features

```python
qiime feature-table filter-features --i-table table-taxa-filter.qza --p-min-samples 2 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 18. Filter samples with less than 500 reads

```python
qiime feature-table filter-samples --i-table table-taxa-filter-no_singles.qza --p-min-frequency 500 --o-filtered-table table-taxa-filter-no_singles.qza

qiime feature-table summarize --i-table table-taxa-filter-no_singles.qza --o-visualization table-taxa-filter-no_singles.qzv --m-sample-metadata-file buzzard_meta.tsv 
```



## 19. Filter representative sequences

```python
qiime feature-table filter-seqs --i-data rep-seqs.qza --i-table filtered-table.qza --o-filtered-data filter-seqs.qza
qiime feature-table tabulate-seqs --i-data filter-seqs.qza --o-visualization filter-seqs.qzv
```



## 20. Build a phylogenetic tree

```python
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences filter-seqs.qza  --o-alignment aligned-seqs.qza --o-masked-alignment masked-aligned-seqs.qza --o-tree unrooted-tree.qza  --o-rooted-tree rooted-tree.qza
```



## 21. Rarefaction curves

```python
qiime diversity alpha-rarefaction --i-table filtered-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 20000 --m-metadata-file buzzard_meta.tsv  --o-visualization alpha-rarefaction.qzv

qiime diversity alpha-rarefaction --i-table filtered-table.qza --i-phylogeny rooted-tree.qza --p-max-depth 2000 --m-metadata-file buzzard_meta.tsv  --o-visualization alpha-rarefaction.qzv
```



#### 21.1 Rarefaction plots (sampling depth 20,000)

![28s-rarefaction-plots](/pics/28s-rarefaction-plots.svg)

#### 21.1 Rarefaction plots (sampling depth 2000)

![28s-rarefaction-plots-2000](/pics/28s-rarefaction-plots-2000.svg)



### 22. Taxa bar plots

```
qiime taxa barplot --i-table filtered-table.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization taxa-bar-plots.qzv 
```



### 22.1. Plot core phylum (present in 70% of samples)

```python
qiime taxa collapse --i-table filtered-table.qza --i-taxonomy taxonomy.qza --o-collapsed-table phylum-table.qza  --p-level 2

qiime feature-table summarize --i-table phylum-table.qza --o-visualization phylum-table.qzv --m-sample-metadata-file buzzard_meta.tsv

qiime feature-table core-features --i-table phylum-table.qza --o-visualization phylum-core-table.qzv 

```

From the visualization file download the tsv relative to features present in 70% of the samples (or any other % that I want)

```python
qiime metadata tabulate --m-input-file core-features-0.700.tsv --o-visualization core-features-0.700.qzv

qiime feature-table filter-features --i-table filtered-table.qza --o-filtered-table table-phylum-core-0.700.qza --m-metadata-file core-features-0.700.tsv

qiime taxa barplot --i-table table-phylum-core-0.700.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization core-phylum-bar-plots.qzv
```



| Feature ID                     | 2%     | 9%     | 25%   | 50%    | 75%     | 91%      | 98%      |
| ------------------------------ | ------ | ------ | ----- | ------ | ------- | -------- | -------- |
| dEukaryota;pAscomycota         | 194.02 | 1249.7 | 2510  | 3538.5 | 5763.75 | 10004.59 | 13126.14 |
| dEukaryota;pBasidiomycota      | 19.08  | 250.31 | 771.5 | 1338.5 | 2443    | 5386.57  | 8451.72  |
| dEukaryota;pPhragmoplastophyta | 0      | 0      | 1.5   | 168    | 549     | 1439.67  | 4006.36  |



### 22.2. Plot core families (present in 70% of samples)

```python
qiime taxa collapse --i-table filtered-table.qza --i-taxonomy taxonomy.qza --o-collapsed-table family-table.qza  --p-level 5

qiime feature-table summarize --i-table -table.qza --o-visualization family-table.qzv --m-sample-metadata-file buzzard_meta.tsv

qiime feature-table core-features --i-table family-table.qza --o-visualization family-core-table.qzv 
```

From the visualization file download the tsv relative to features present in 70% of the samples (or any other % that I want)

```python
qiime metadata tabulate --m-input-file core-features-0.700.tsv --o-visualization core-features-0.700.qzv

qiime feature-table filter-features --i-table filtered-table.qza --o-filtered-table table-family-core-0.700.qza --m-metadata-file core-features-0.700.tsv

qiime taxa barplot --i-table table-family-core-0.700.qza --i-taxonomy taxonomy.qza --m-metadata-file buzzard_meta.tsv --o-visualization core-family-bar-plots.qzv
```



| Feature ID                                                   | 2%   | 9%     | 25%   | 50%   | 75%    | 91%     | 98%     |
| ------------------------------------------------------------ | ---- | ------ | ----- | ----- | ------ | ------- | ------- |
| dEukaryota;pBasidiomycota;cExobasidiomycetes;oExobasidiomycetes;f__Exobasidiomycetes | 0    | 105.44 | 293   | 775   | 1621.5 | 3861.9  | 7110.5  |
| dEukaryota;pAscomycota;cDothideomycetes;oCapnodiales;f__Capnodiales | 0    | 17.66  | 220   | 457   | 797    | 1134.9  | 1368.08 |
| dEukaryota;pAscomycota;cDothideomycetes;oCapnodiales;__      | 0    | 14.44  | 111.5 | 236.5 | 520.25 | 1088.58 | 2432.76 |
| dEukaryota;pAscomycota;cDothideomycetes;oCapnodiales;f__Cladosporiaceae | 0    | 0      | 29.75 | 161   | 466.25 | 1123.67 | 3475.32 |
| dEukaryota;pAscomycota;;;__                                  | 0    | 0      | 10.75 | 144   | 354.5  | 587.78  | 1235.62 |
| dEukaryota;pPhragmoplastophyta;cPhragmoplastophyta;oPhragmoplastophyta;f__Phragmoplastophyta | 0    | 0      | 0     | 106.5 | 416.75 | 1247.84 | 3978.52 |



## 23. Calculate alpha diversity metrics and rarefy the data-set

```python
# Calculates observed features, shannon diversity, Faith PD and retrieves rarefied table
qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table filtered-table.qza --p-sampling-depth 2000 --m-metadata-file buzzard_meta.tsv --output-dir alpha-metrics-results

#Add alpha diversity metrics to the metadata
qiime metadata tabulate --m-input-file buzzard_metadata.tsv --m-input-file shannon_vector.qza --m-input-file observed_features.qza --m-input-file faith_pd_vector.qza --o-visualization buzzard_meta_alpha.qzv

```

From the "buzzard_meta_alpha.qzv" one can extract the tsv file with the metadata + alpha metrics.



## 24. Filter samples with less than 2000 reads (unrarefied table for beta analysis)

```python
qiime feature-table filter-samples --i-table filtered-table.qza --p-min-frequency 4000 --o-filtered-table beta-table.qza

qiime feature-table summarize --i-table beta-table.qza --o-visualization beta-table.qzv --m-sample-metadata-file buzzard_meta_alpha.tsv
```



### 24.1. Table Summary

|                               |  Sample   |
| ----------------------------- | :-------: |
| **Number of samples**         |    180    |
| **Number of features (ASVs)** |   1,770   |
| **Total frequency**           | 1,478,279 |



### 24.2. Summary of final reads per sample

|                       | **Frequency** |
| --------------------- | :-----------: |
| **Minimum frequency** |   2,036.00    |
| **1st quartile**      |   5,388.25    |
| **Median frequency**  |   7,138.50    |
| **3rd quartile**      |   9,985.75    |
| **Maximum frequency** |   33,249.00   |
| **Mean frequency**    |   8,212.66    |



### 24.2. Sequence length statistics

| Sequence Count | Min Length | Max Length | Mean Length | Range | Standard Deviation |
| :------------: | :--------: | :--------: | :---------: | :---: | :----------------: |
|      1770      |    230     |    415     |   413.84    |  185  |       11.08        |
