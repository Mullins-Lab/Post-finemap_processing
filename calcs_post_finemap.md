#### Make sure to the paths accordingly; paths below are indicative
Input for step 1) are all the .gz files from running Polyfun-based fine-mapping.

Input for step 2) and 3) are the outputted merged dataframes from step 1).

### Step 1.
### Merge all individual fine-mapping jobs per datasets into one file and add a LOCUS column

```
files <- list.files(path="/Users/koromm03/Downloads/polyfun_finemap_UKB_finemap/", pattern="daner_bip_pgc3_finemap_finemap.*.gz", full.names=TRUE, recursive=FALSE)
files
df <- data.frame()

library(data.table)
library(tidyverse)
for (file in files){
  f= fread(file, header=TRUE) 
  locus <- strsplit(file, split='[..]')[[1]][2] 
  print(locus) 
  f$LOCUS= locus
  df= rbind(df,f, fill=TRUE)
}

write_csv(df, file='/Users/koromm03/Downloads/polyfun_finemap_UKB_finemap/daner_bip_pgc3_polyfun_finemap.merged.csv')
```

### Step 2.
### Count the size of CSs per each finemapped locus
The input is a merged dataframe with all finemapping results from step 1.

```
library(ggplot2)
library(tidyverse)
library(vroom)
library(dplyr)

dat <- read_csv("/Users/koromm03/Desktop/dfs_ashvin_loop/polyfun_finemap_HRC.csv") %>%
  janitor::clean_names()

`%nin%` = negate(`%in%`)

test <- dat %>% 
  filter(credible_set != 0) %>%
  group_by(locus) %>%
  count(credible_set)


test %>% 
  pivot_wider(names_from = credible_set, values_from = n) %>%
  write_csv(., "~/Downloads/test.csv")
 ``` 
  
### Step 3.
### Create filtered results (subsetting SNPs per PIP and if within credible sets)

```
library(tidyverse)
library(dplyr)

finemap_files <- list.files(path='/Users/koromm03/Desktop/polyfun_processing/merged_dfs/', pattern = "*.tsv", full.names=TRUE, recursive=FALSE)
for(i in 1:length(finemap_files))
{
dat <- read_tsv(finemap_files[i]) %>%
  janitor::clean_names()

`%nin%` = negate(`%in%`)

test <- dat %>% 
  filter(credible_set != 0) %>%
  group_by(locus) %>%
  filter(pip>0.95)


test %>% 
  #pivot_wider(names_from = pip, values_from = pip) %>%
  write_tsv(., file = file.path(paste(finemap_files[i],"_filt_process_095.tsv")))
}
```
  
 
  
  
  
  
  
