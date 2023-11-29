# Alpha diversity statistical analysis





## Table of Contents

[TOC]

------





# A) 16S rRNA alpha diversity statistical analysis



### 1. Shannon diversity Index

```R
#Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library (performance)
library(datawizard)
library(car)
library(effects)
library(ggpubr)
library(jtools)

#Load dataset
metadata <- readRDS("metadata-rarefied.rds")

# Calculate age and body condition index
metadata$bci_two<-resid(glm(weight~ log10(wing) + sex, gaussian, metadata, na.action="na.exclude")) #calculate body condition
metadata$std_bci <- scale(metadata$bci_two) # scale bci values

metadata$age_days <- buteo_age(df = metadata, wing = "wing", sex = TRUE, unit = c("cm"), .plot = F, decimals = 2,.show_model = T)$fit
metadata$std_age <- scale(metadata$age_days) # scale age values

saveRDS(metadata,"metadata-rarefied.rds")
```



#### 1.1 Model Shannon

```R
model_shannon <- lmer(shannon_entropy ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)
```



#### 1.2 Check Normality

```R
> check_normality(model_shannon)
OK: residuals appear as normally distributed (p = 0.394).
```



#### 1.3  Model Diagnostics

```R
check_model(model_shannon)
```

![shannon-model-diagnostics](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/shannon-model-diagnostics.svg)



#### 1.4 Model Summary

```R
> summary(model_shannon_final)

Linear mixed model fit by REML ['lmerMod']
Formula: shannon_entropy ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1 | nest/ring_number)
   Data: metadata

REML criterion at convergence: 475.7

Scaled residuals: 
     Min       1Q        Median       3Q       Max 
-2.38433    -0.63636   -0.03356    0.65990   2.21304 

Random effects:
 Groups           Name        Variance   Std.Dev.
 ring_number:nest (Intercept) 0.00000     0.0000  
 nest             (Intercept) 0.04932     0.2221  
 Residual                     0.40710     0.6380  
Number of obs: 226, groups:  ring_number:nest, 117; nest, 54

Fixed effects:
               Estimate    Std. Error    t value
(Intercept)   4.3068106     0.1058171     40.701
std_age      -0.1652697     0.0529112     -3.124
std_bci       0.0001679     0.0511950      0.003
rank2        -0.0925053     0.0964370     -0.959
rank3        -0.1613338     0.1423461     -1.133
rank4         0.2139437     0.5077520      0.421
sexM         -0.1144732     0.0963586     -1.188
year2021      0.0683771     0.1981118      0.345
habitatsouth -0.1657539     0.1744792     -0.950
habitatteuto -0.2976098     0.2869742     -1.037
lbinom1       0.0787340     0.1036844      0.759

Correlation of Fixed Effects:
            (Intr) std_ag std_b_ rank2  rank3  rank4  sexM   yr2021 hbttst hbtttt
std_age      0.143                                                               
std_bci     -0.087  0.023                                                        
rank2       -0.404  0.196  0.101                                                 
rank3       -0.246  0.216  0.210  0.378                                          
rank4       -0.025  0.112  0.240  0.132  0.167                                   
sexM        -0.388 -0.043 -0.088 -0.035 -0.146 -0.138                            
year2021    -0.177 -0.106  0.100 -0.028  0.029  0.041 -0.138                     
habitatsoth -0.143 -0.013  0.156  0.048  0.168  0.094 -0.207  0.134              
habitatteut -0.087 -0.137  0.239  0.074  0.114  0.077 -0.082 -0.086  0.124       
lbinom1     -0.521 -0.371 -0.036 -0.032 -0.079 -0.082 -0.043  0.142 -0.038 -0.065
optimizer (nloptwrap) convergence code: 0 (OK)
boundary (singular) fit: see help('isSingular')
```



#### 1.5 Significance values

```R
> Anova(model_shannon)

Analysis of Deviance Table (Type II Wald chisquare tests)

Response: shannon_entropy
             Chisq    Df    Pr(>Chisq)   
std_age     9.7564    1      0.001787 **
std_bci     0.0000    1      0.997383   
rank        2.0552    3      0.561028   
sex         1.4113    1      0.234836   
year        0.1191    1      0.729986   
habitat     1.7600    2      0.414783   
lbinom      0.5766    1      0.447636   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



#### 1.6 Marginal and conditional R-squared

```R
 r.squaredGLMM(model_shannon)
           R2m         R2c
[1,]   0.0758452    0.1757042
```



#### 1.7 Plot model effects

```R
plot(allEffects(model_shannon))
```



![shannon-model-effects](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/shannon-model-effects.svg)



#### 1.8 Repeatability Analysis

```R
> rpts <- rpt(shannon_entropy ~ std_age + std_bci + rank + sex  + year + habitat + (1|nest) + (1|ring_number), 
            grname = c("nest", "ring_number"), data = metadata, 
            datatype = "Gaussian",adjusted = TRUE,
            nboot = 1000, npermut = 1000)

> summary(rpts)



Repeatability estimation using the lmm method

Call = rpt(formula = shannon_entropy ~ std_age + std_bci_two + rank + sex + year + habitat + (1 | nest) + (1 | ring_number), grname = c("nest", "ring_number"), data = metadata, datatype = "Gaussian", nboot = 1000, npermut = 1000, adjusted = TRUE)

Data: 228 observations
----------------------------------------

nest (54 groups)

Repeatability estimation overview: 
      R      SE       2.5%    97.5%   P_permut   LRT_P
  0.109    0.0618      0      0.224    0.017     0.028

Bootstrapping and Permutation test: 
           N   Mean     Median     2.5%   97.5%
boot     1000 0.0940   9.14e-02      0   0.2237
permut   1000 0.0156   1.12e-12      0   0.0932

Likelihood ratio test: 
logLik full model = -239.6265
logLik red. model = -241.4516
D  = 3.65, df = 1, P = 0.028

----------------------------------------


ring_number (117 groups)

Repeatability estimation overview: 
      R     SE        2.5%   97.5%     P_permut    LRT_P
      0    0.0633      0     0.219        1        0.5

Bootstrapping and Permutation test: 
            N     Mean    Median      2.5%    97.5%
boot     1000   0.0430   1.13e-08      0      0.219
permut   1000   0.0347   1.69e-11      0      0.175

Likelihood ratio test: 
logLik full model = -239.6265
logLik red. model = -239.6265
D  = 2.27e-13, df = 1, P = 0.5

----------------------------------------
```



```R
plot(rpts, grname="nest", type="boot", cex.main=0.8, col = "#ECEFF4")
plot(rpts, grname="ring_number", type="boot", cex.main=0.8, col = "#ECEFF4")
```

![](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/shannon-repeatability-nest.svg)

![shannon-repeatability-ID](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/shannon-repeatability-ID.svg)

### Faith phylogenetic diversity

#### 2. Log transform Faith 

(model residuals not normal distributed)

```R
model_faith <- lmer(faith_pd ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)

> check_normality(model_faith)
Warning: Non-normality of residuals detected (p = 0.003).

#Log transform faith
metadata$log_faith <- log10(metadata$faith_pd)
```



#### 2.1 Model Faith PD

```R
model_faith <- lmer(log_faith ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)
```



#### 2.2 Check Normality

```R
> check_normality(model_faith)
OK: residuals appear as normally distributed (p = 0.081).
```



#### 2.3 Model Diagnostics

```R
check_model(model_faith)
```



![model-diagnostics-faith](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/faith-model-diagnostics.svg)



#### 2.4 Model Summary

```R
> summary(model_faith)

Linear mixed model fit by REML ['lmerMod']
Formula: log_faith ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1 | nest/ring_number)
   Data: metadata

REML criterion at convergence: -225.5

Scaled residuals: 
    Min      1Q     Median      3Q     Max 
-3.1497   -0.5645  0.1041     0.6145  2.5580 

Random effects:
 Groups           Name        Variance    Std.Dev.
 ring_number:nest (Intercept) 0.0000000   0.00000 
 nest             (Intercept) 0.0006666   0.02582 
 Residual                     0.0164990   0.12845 
Number of obs: 226, groups:  ring_number:nest, 117; nest, 54

Fixed effects:
               Estimate   Std. Error    t value
(Intercept)   0.9407717    0.0199731    47.102
std_age      -0.0285771    0.0102204    -2.796
std_bci      -0.0067600    0.0097165    -0.696
rank2        -0.0189343    0.0192559    -0.983
rank3        -0.0154191    0.0279228    -0.552
rank4         0.2347500    0.0988584     2.375
sexM         -0.0360834    0.0186700    -1.933
year2021     -0.0134265    0.0354902    -0.378
habitatsouth -0.0029790    0.0315251    -0.094
habitatteuto  0.0008674    0.0526736     0.016
lbinom1       0.0146748    0.0200224     0.733

Correlation of Fixed Effects:
            (Intr) std_ag std_b_ rank2  rank3  rank4  sexM   yr2021 hbttst hbtttt
std_age      0.164                                                               
std_bci     -0.079  0.018                                                        
rank2       -0.428  0.188  0.090                                                 
rank3       -0.280  0.200  0.210  0.374                                          
rank4       -0.027  0.097  0.246  0.127  0.155                                   
sexM        -0.394 -0.052 -0.097 -0.037 -0.137 -0.129                            
year2021    -0.166 -0.115  0.105 -0.031  0.035  0.045 -0.148                     
habitatsoth -0.124 -0.011  0.172  0.049  0.187  0.104 -0.220  0.136              
habitatteut -0.086 -0.142  0.251  0.076  0.125  0.087 -0.082 -0.071  0.130       
lbinom1     -0.525 -0.394 -0.043 -0.034 -0.070 -0.095 -0.049  0.152 -0.041 -0.063
optimizer (nloptwrap) convergence code: 0 (OK)
boundary (singular) fit: see help('isSingular')
```



#### 2.5 Significance values

```R
> Anova(model_faith)
Analysis of Deviance Table (Type II Wald chisquare tests)

Response: log_faith
             Chisq   Df   Pr(>Chisq)   
std_age     7.8181   1     0.005173 **
std_bci     0.4840   1     0.486605   
rank        7.5620   3     0.055986 . 
sex         3.7353   1     0.053274 . 
year        0.1431   1     0.705197   
habitat     0.0098   2     0.995125   
lbinom      0.5372   1     0.463609   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



#### 2.6 Marginal and Conditional R-squared

```R
> r.squaredGLMM(model_faith)
            R2m        R2c
[1,]   0.09268855   0.1279241
```



#### 2.7 Plot model effects

```
plot(allEffects(model_faith))
```

![faith-model-effects](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/faith-model-effects.svg)



#### 2.8 Repeatability Analysis

```R
> rpts2 <- rpt(log_faith ~ std_age + std_bci + rank + sex  + year + habitat + (1|nest) + (1|ring_number), 
            grname = c("nest", "ring_number"), data = metadata, 
            datatype = "Gaussian",adjusted = TRUE,
            nboot = 1000, npermut = 1000)

> summary(rpts2)

Repeatability estimation using the lmm method

Call = rpt(formula = log_faith ~ std_age + std_bci + rank + sex + year + habitat + (1 | nest) + (1 | ring_number), grname = c("nest", "ring_number"), data = metadata, datatype = "Gaussian", nboot = 1000, npermut = 1000, adjusted = TRUE)

Data: 228 observations
----------------------------------------

nest (54 groups)

Repeatability estimation overview: 
      R     SE      2.5%   97.5%    P_permut   LRT_P
 0.0385   0.046      0    0.158       0.17     0.222

Bootstrapping and Permutation test: 
            N     Mean     Median      2.5%      97.5%
boot     1000   0.0404    2.57e-02      0       0.1575
permut   1000   0.0167    2.04e-11      0       0.0987

Likelihood ratio test: 
logLik full model = 117.2883
logLik red. model = 116.9959
D  = 0.585, df = 1, P = 0.222

----------------------------------------


ring_number (117 groups)

Repeatability estimation overview: 
      R     SE       2.5%    97.5%    P_permut  LRT_P
      0   0.0551      0      0.188        1      0.5

Bootstrapping and Permutation test: 
            N   Mean   Median   2.5%  97.5%
boot     1000 0.0350 3.76e-11      0  0.188
permut   1000 0.0329 1.42e-12      0  0.192

Likelihood ratio test: 
logLik full model = 117.2883
logLik red. model = 117.2883
D  = 2.27e-13, df = 1, P = 0.5

----------------------------------------
```



```R
plot(rpts2, grname="nest", type="boot", cex.main=0.8, col = "#ECEFF4")
plot(rpts2, grname="ring_number", type="boot", cex.main=0.8, col = "#ECEFF4")
```





![](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/16S/model-results/faith-repeatability-nest.svg)

#### ![faith-repeatability-ID](/home/localadmin/microbiome-analysis/alpha-diversity/16S/model-results/faith-repeatability-ID.svg)

------





# B) 28S rRNA alpha diversity statistical analysis



### 1. Shannon diversity Index

```R
#Load libraries
library(tidyverse)
library(lme4)
library(MuMIn)
library (performance)
library(datawizard)
library(car)
library(effects)
library(ggpubr)
library(multcomp)
library(jtools)

#Load dataset
metadata <- readRDS("metadata-rarefied.rds")

# Calculate age and body condition index
metadata$bci_two<-resid(glm(weight~ log10(wing) + sex, gaussian, metadata, na.action="na.exclude")) #calculate body condition
metadata$std_bci <- scale(metadata$bci_two) # scale bci values

metadata$age_days <- buteo_age(df = metadata, wing = "wing", sex = TRUE, unit = c("cm"), .plot = F, decimals = 2,.show_model = T)$fit
metadata$std_age <- scale(metadata$age_days) # scale age values

saveRDS(metadata,"metadata-rarefied.rds")
```



#### 1.1. Transform Shannon 

```R
> model_shannon <- lmer(shannon_entropy ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)

> check_normality(model_shannon)
Warning: Non-normality of residuals detected (p < .001).

# Reflect and log transform shannon
metadata$log_shannon <- log10(max(metadata$shannon_entropy+1) - metadata$shannon_entropy)

# Data reflection changes direction of relationships
metadata$log_shannon <- -metadata$log_shannon # change directions of relatioships again
```



#### 1.2. Model Faith PD

```R
model_shannon <- lmer(log_shannon ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)
```



#### 1.3. Check Normality

```R
> check_normality(model_shannon_final_prev)
Warning: Non-normality of residuals detected (p = 0.046). # residuals still not normal distributed
```



#### 1.4. Model Diagnostics

```R
check_model(model_shannon)# normality of residuals identified by visual inspection
```

![shannon-model-diagnostics](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/shannon-model-diagnostics.svg)



#### 1.5. Model Summary

```R
> summary(model_shannon)
Linear mixed model fit by REML ['lmerMod']
Formula: log_shannon ~ std_age + std_bci + sex + rank + habitat + year + lbinom + (1 | nest/ring_number)
   Data: metadata

REML criterion at convergence: -28.9

Scaled residuals: 
    Min      1Q      Median      3Q        Max 
-1.7910   -0.7103   -0.0168    0.6698    2.4352 

Random effects:
 Groups           Name        Variance  Std.Dev. 
 ring_number:nest (Intercept) 2.032e-09 4.508e-05
 nest             (Intercept) 2.236e-03 4.729e-02
 Residual                     3.754e-02 1.938e-01
Number of obs: 177, groups:  ring_number:nest, 108; nest, 54

Fixed effects:
              Estimate    Std. Error    t value
(Intercept)  -0.520591     0.034071     -15.280
std_age      -0.035594     0.017628     -2.019
std_bci      -0.004052     0.016615     -0.244
sexM          0.016542     0.032398      0.511
rank2        -0.004258     0.033979     -0.125
rank3         0.029871     0.045397      0.658
rank4         0.210610     0.152755      1.379
habitatsouth  0.068654     0.053602      1.281
habitatteuto  0.112459     0.108764      1.034
year2021      0.003454     0.061376      0.056
lbinom1       0.019645     0.034963      0.562

Correlation of Fixed Effects:
            (Intr) std_ag std_b_ sexM   rank2  rank3  rank4  hbttst hbtttt yr2021
std_age      0.186                                                               
std_bci     -0.023  0.002                                                        
sexM        -0.408 -0.103 -0.112                                                 
rank2       -0.426  0.143  0.086 -0.062                                          
rank3       -0.265  0.231  0.206 -0.101  0.372                                   
rank4       -0.014  0.114  0.251 -0.136  0.123  0.178                            
habitatsoth -0.111  0.047  0.138 -0.164  0.026  0.183  0.106                     
habitatteut -0.123 -0.151  0.184  0.036  0.123  0.113  0.066  0.078              
year2021    -0.109 -0.060  0.095 -0.155 -0.078  0.009  0.049  0.125 -0.091       
lbinom1     -0.504 -0.398 -0.024 -0.099  0.032 -0.120 -0.105 -0.079 -0.070  0.107
```



#### 1.6. Significance values

```R
> Anova(model_shannon)
Analysis of Deviance Table (Type II Wald chisquare tests)

Response: log_shannon
             Chisq  Df  Pr(>Chisq)  
std_age     4.0769  1    0.04347 *
std_bci     0.0595  1    0.80735  
sex         0.2607  1    0.60965  
rank        2.3064  3    0.51129  
habitat     2.5173  2    0.28404  
year        0.0032  1    0.95512  
lbinom      0.3157  1    0.57421  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



#### 1.7. Marginal and conditional R-squared

```R
> r.squaredGLMM(model_shannon)
            R2m          R2c
[1,]    0.05739074    0.1103858
```



#### 1.8. Plot model effects

```R
plot(allEffects(model_shannon))
```

![shannon-model-effects](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/shannon-model-effects.svg)



#### 1.9. Repeatability Analysis

```R
> rpts <- rpt(shannon_entropy ~ std_age + std_bci + rank + sex  + year + habitat + (1|nest) + (1|ring_number), 
            grname = c("nest", "ring_number"), data = metadata, 
            datatype = "Gaussian",adjusted = TRUE,
            nboot = 1000, npermut = 1000)

> summary(rpts)

Repeatability estimation using the lmm method

Call = rpt(formula = log_shannon ~ std_age + std_bci + rank + sex + year + habitat + (1 | nest) + (1 | ring_number), grname = c("nest", "ring_number"), data = metadata, datatype = "Gaussian", nboot = 1000, npermut = 1000, adjusted = TRUE)

Data: 179 observations
----------------------------------------

nest (54 groups)

Repeatability estimation overview: 
      R     SE      2.5%    97.5%   P_permut   LRT_P
 0.0701   0.0643     0      0.217    0.108     0.18

Bootstrapping and Permutation test: 
            N     Mean     Median   2.5%   97.5%
boot     1000   0.0634   5.07e-02    0     0.217
permut   1000   0.0210   6.95e-12    0     0.123

Likelihood ratio test: 
logLik full model = 16.23022
logLik red. model = 15.81233
D  = 0.836, df = 1, P = 0.18

----------------------------------------


ring_number (109 groups)

Repeatability estimation overview: 
      R       SE      2.5%   97.5%   P_permut   LRT_P
      0     0.0745     0     0.243       1        1

Bootstrapping and Permutation test: 
            N   Mean      Median     2.5%  97.5%
boot     1000   0.0486   3.10e-11     0    0.243
permut   1000   0.0396   1.42e-12     0    0.233

Likelihood ratio test: 
logLik full model = 16.23022
logLik red. model = 16.23022
D  = -1.51e-12, df = 1, P = 1

----------------------------------------



```



```R
plot(rpts, grname="nest", type="boot", cex.main=0.8, col = "#ECEFF4")
plot(rpts, grname="ring_number", type="boot", cex.main=0.8, col = "#ECEFF4")
```

![shannon-repeatability-nest](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/shannon-repeatability-nest.svg)

![shannon-repeatability-ID](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/shannon-repeatability-ID.svg)



### 2. Faith phylogenetic diversity



#### 2.1. Log transform Faith 

(model residuals not normal distributed)

```R
model_faith <- lmer(faith_pd ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)

> check_normality(model_faith_final_lbinom)
Warning: Non-normality of residuals detected (p < .001).

#Log transform faith
metadata$log_faith <- log10(metadata$faith_pd)
```



#### 2.2. Model Faith PD

```R
model_faith <- lmer(log_faith ~ std_age + std_bci + rank + sex + year + habitat + lbinom + (1|nest/ring_number), data = metadata)
```



#### 2.3. Check Normality

```R
> check_normality(model_faith)
OK: residuals appear as normally distributed (p = 0.661).
```



#### 2.4. Model Diagnostics

```R
check_model(model_faith)
```

![faith-model-diagnostics](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/faith-model-diagnostics.svg)



#### 2.5. Model Summary

```R
> summary(model_faith)

Linear mixed model fit by REML ['lmerMod']
Formula: log_faith ~ std_age + std_bci + rank + sex + year + habitat +      lbinom + (1 | nest/ring_number)
   Data: metadata

REML criterion at convergence: -10

Scaled residuals: 
     Min       1Q     Median     3Q      Max 
-2.62275   -0.65997  0.01234  0.63677  2.32240 

Random effects:
 Groups           Name        Variance Std.Dev.
 ring_number:nest (Intercept) 0.00000  0.0000  
 nest             (Intercept) 0.00000  0.0000  
 Residual                     0.04422  0.2103  
Number of obs: 177, groups:  ring_number:nest, 108; nest, 54

Fixed effects:
              Estimate  Std. Error    t value
(Intercept)   0.403635   0.035173     11.476
std_age      -0.020558   0.018479     -1.113
std_bci      -0.019744   0.017132     -1.152
rank2         0.007209   0.036431      0.198
rank3        -0.014333   0.048106     -0.298
rank4         0.315343   0.158692      1.987
sexM          0.008349   0.033891      0.246
year2021      0.004901   0.060742      0.081
habitatsouth  0.094678   0.053631      1.765
habitatteuto  0.219712   0.113001      1.944
lbinom1       0.043000   0.036659      1.173

Correlation of Fixed Effects:
            (Intr) std_ag std_b_ rank2  rank3  rank4  sexM   yr2021 hbttst hbtttt
std_age      0.205                                                               
std_bci     -0.013 -0.006                                                        
rank2       -0.443  0.137  0.072                                                 
rank3       -0.293  0.224  0.204  0.372                                          
rank4       -0.014  0.101  0.266  0.125  0.165                                   
sexM        -0.414 -0.117 -0.111 -0.065 -0.093 -0.128                            
year2021    -0.098 -0.065  0.096 -0.083  0.007  0.049 -0.153                     
habitatsoth -0.095  0.048  0.147  0.025  0.196  0.112 -0.169  0.123              
habitatteut -0.125 -0.151  0.184  0.127  0.119  0.073  0.042 -0.085  0.078       
lbinom1     -0.504 -0.418 -0.036  0.030 -0.116 -0.118 -0.099  0.111 -0.087 -0.070
optimizer (nloptwrap) convergence code: 0 (OK)
boundary (singular) fit: see help('isSingular')
```



#### 2.6. Significance values

```R
> Anova(model_faith_final)
Analysis of Deviance Table (Type II Wald chisquare tests)

Response: log_faith
             Chisq Df  Pr(>Chisq)  
std_age     1.2377  1    0.26591  
std_bci     1.3281  1    0.24914  
rank        4.3856  3    0.22272  
sex         0.0607  1    0.80540  
year        0.0065  1    0.93569  
habitat     6.3988  2    0.04079 *
lbinom      1.3759  1    0.24081  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```



##### 2.6.1. Multiple comparison test for "Habitat"

```R
> library(multcomp)

> multicomp <- glht(model_faith, linfct = mcp(habitat="Tukey")) #multicomparisson for linear models

> confint(glht(model_faith_final_lbinom, mcp(habitat="Tukey")))

	 Simultaneous Confidence Intervals

Multiple Comparisons of Means: Tukey Contrasts


Fit: lmer(formula = log_faith ~ std_age + std_bci + rank + sex + 
    year + habitat + lbinom + (1 | nest/ring_number), data = metadata)

Quantile = 2.3042
95% family-wise confidence level
 

Linear Hypotheses:
                   Estimate       lwr       upr     
south - north == 0  0.09468    -0.02890   0.21826
teuto - north == 0  0.21971    -0.04067   0.48009
teuto - south == 0  0.12503    -0.15431   0.40438

> summary(multicomp, test = adjusted("holm"))

 Simultaneous Tests for General Linear Hypotheses

Multiple Comparisons of Means: Tukey Contrasts


Fit: lmer(formula = log_faith ~ std_age + std_bci_two + rank + sex + 
    year + habitat + lbinom + (1 | nest/ring_number), data = metadata)

Linear Hypotheses:
                   Estimate   Std. Error  z value  Pr(>|z|)
south - north == 0  0.09468     0.05363    1.765    0.116
teuto - north == 0  0.21971     0.11300    1.944    0.116
teuto - south == 0  0.12503     0.12123    1.031    0.302
(Adjusted p values reported -- BH method)
```



#### 2.7. Marginal and Conditional R-squared

```R
> r.squaredGLMM(model_faith)
           R2m          R2c
[1,]    0.09305071   0.09305071
```



#### 2.8. Plot model effects

```
plot(allEffects(model_faith))
```

![faith-model-effects](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/faith-model-effects.svg)



#### 2.9. Repeatability Analysis

```R
> rpts2 <- rpt(log_faith ~ std_age + std_bci + rank + sex  + year + habitat + (1|nest) + (1|ring_number), 
            grname = c("nest", "ring_number"), data = metadata, 
            datatype = "Gaussian",adjusted = TRUE,
            nboot = 1000, npermut = 1000)

> summary(rpts2)

Repeatability estimation using the lmm method

Call = rpt(formula = log_faith ~ std_age + std_bci + rank + sex + year + habitat + (1 | nest) + (1 | ring_number), grname = c("nest", "ring_number"), data = metadata, datatype = "Gaussian", nboot = 1000, npermut = 1000, adjusted = TRUE)

Data: 179 observations
----------------------------------------

nest (54 groups)

Repeatability estimation overview: 
      R     SE      2.5%  97.5%   P_permut   LRT_P
      0   0.0366     0    0.128       1        1

Bootstrapping and Permutation test: 
            N    Mean    Median    2.5%   97.5%
boot     1000  0.0204   1.22e-14    0     0.128
permut   1000  0.0214   1.12e-13    0     0.145

Likelihood ratio test: 
logLik full model = 6.67855
logLik red. model = 6.67855
D  = 0, df = 1, P = 1

----------------------------------------


ring_number (109 groups)

Repeatability estimation overview: 
      R     SE     2.5%  97.5%   P_permut   LRT_P
      0   0.0653    0    0.237       1        1

Bootstrapping and Permutation test: 
            N   Mean   Median   2.5%  97.5%
boot     1000 0.0376 2.40e-13      0  0.237
permut   1000 0.0321 3.66e-17      0  0.223

Likelihood ratio test: 
logLik full model = 6.67855
logLik red. model = 6.67855
D  = 0, df = 1, P = 1

----------------------------------------

```



```R
plot(rpts2, grname="nest", type="boot", cex.main=0.8, col = "#ECEFF4")
plot(rpts2, grname="ring_number", type="boot", cex.main=0.8, col = "#ECEFF4")
```



![faith-repeatability-nest](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/faith-repeatability-nest.svg)

![faith-repeatability-ID](/home/localadmin/MEGA/PhD/Buzzards/microbiome-analysis/alpha-diversity/28S/concat_230_185/model-results/faith-repeatability-ID.svg)