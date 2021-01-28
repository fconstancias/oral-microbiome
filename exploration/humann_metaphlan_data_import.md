---
title: "Oral Microbiome - humann  metaphlan  import 2021-01-28"
author: "Florentin CONSTANCIAS"
date: "2021-01-28"
output: 
  html_document: 
    toc: yes
    toc_depth: 4
    keep_md: yes
---






```r
require(tidyverse); packageVersion("tidyverse")
```

```
## Loading required package: tidyverse
```

```
## ── Attaching packages ────────────────────────────────── tidyverse 1.3.0.9000 ──
```

```
## ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
## ✓ tibble  3.0.5     ✓ dplyr   1.0.3
## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
## ✓ readr   1.4.0     ✓ forcats 0.5.0
```

```
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
```

```
## [1] '1.3.0.9000'
```

```r
require(phyloseq); packageVersion("phyloseq")
```

```
## Loading required package: phyloseq
```

```
## [1] '1.34.0'
```




```r
sapply(c(biobakery_function, alpha, beta, taxa, norm), source, chdir = TRUE)
```

Inputs:


```r
metpahlan <- "~/Documents/GitHub/oral-microbiome/data/processed/humann/merged_metaphlan3_humann3.tsv"

pathway_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/pathabundance_cpm_joined_tables.tsv"
pathway_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/pathabundance_cpm_joined_tables.tsv"

kegg_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv"
kegg_RNA <- "~/Documents/GitHub/oral-microbiome/data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_ko_renamed_kegg-orthology.tsv"

l4c_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_level4ec_renamed_ec.tsv"
l4c_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_level4ec_renamed_ec.tsv"

go_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_go_renamed_go.tsv"
go_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_go_renamed_go.tsv"

infogo_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_go_renamed_infogo1000.tsv"
infogo_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_go_renamed_infogo1000.tsv"


metac_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_rxn_renamed_metacyc-rxn.tsv"
metac_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_rxn_renamed_metacyc-rxn.tsv"

pfam_DNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/DNA/genefamilies_cpm_joined_tables_uniref90_pfam_renamed_pfam.tsv"
pfam_RNA <- "~/Documents/GitHub/oral-microbiome//data/processed/humann/RNA/genefamilies_cpm_joined_tables_uniref90_pfam_renamed_pfam.tsv"

metadata_path <- "~/Documents/GitHub/oral-microbiome/data/metadata_all_DNA_RNA.xlsx"
```

Data loading:


```r
metadata_path %>% readxl::read_xlsx(col_types = c(rep("text", 7),rep("numeric", 11),rep("text", 5), rep("numeric", 13))) -> metadata
```

Since we are going to import a lot of humman tables let's write a function for that:


```r
import_human_2phyloseq <- function(DNA_file_path,
                                   RNA_file_path,
                                   type = '# Gene Family',
                                   n_rows = Inf,
                                   str_rm_a = "_cat_Abundance-CPM_DNA",
                                   str_rplace_a = "",
                                   str_rm_b = "_cat_Abundance-CPM_RNA",
                                   str_replace_b = "_RNA",
                                   sample_column = "Sample_ID"){
  
  humann_DNA_RNA_2phyloseq(DNA_humann_2df = DNA_file_path %>%
                             humann_2df(type = type, n_rows = n_rows) %>%
                             clean_humann_df(),
                           RNA_humann_2df = RNA_file_path %>%
                             humann_2df(type = type, n_rows = n_rows) %>%
                             clean_humann_df()) %>%
    clean_phyloseq_sample_names(str_rm = str_rm_a,
                                str_replace = str_rplace_a) %>%
    clean_phyloseq_sample_names(str_rm = str_rm_b,
                                str_replace = str_replace_b) %>%
    physeq_add_metadata(metadata,
                        sample_column = sample_column) -> physeq
  
  return(physeq)
}

physeq_count_normalize <- function(physeq,
                                   sample_idx = "Sample_ID",
                                   value_idx = "ReadNorm",
                                   return_both = FALSE){
  
  HTSSIP::OTU_qPCR_trans(physeq, 
                         physeq %>% 
                           sample_data() %>% 
                           data.frame(), 
                         sample_idx = sample_idx, value_idx = value_idx) -> physeq_count
  
  if(return_both == TRUE){
    out <- list("physeq" = physeq,
                "physeq_count" = physeq_count)
    return(out)
  }else{
    return(physeq_count)
  }
}

save_object_name <- function(object,
                             path = paste0(here::here(), "/data/processed/"),
                             return_input = FALSE){
  saveRDS(object,
          paste0(path, "/", deparse(substitute(object)), ".RDS"))
  if(return_input == TRUE){
    return(object)
  }
}
```


```r
metpahlan %>%
  metaphlan_2phyloseq(id = 'clade_name', skip_col = 1, metadata = "none") %>%
  clean_phyloseq_sample_names(str_rm = "_DNA_cat_metaphlan_bugs_list",
                              str_replace = "") %>%
  physeq_add_metadata(metadata %>% filter(Type == "DNA"),
                      sample_column = "Sample") %>%
  physeq_count_normalize(return_both = TRUE)  -> metaphlan_ps

metaphlan_ps
```

```
## $physeq
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 283 taxa and 59 samples ]:
## sample_data() Sample Data:        [ 59 samples by 36 sample variables ]:
## tax_table()   Taxonomy Table:     [ 283 taxa by 7 taxonomic ranks ]:
## phy_tree()    Phylogenetic Tree:  [ 283 tips and 282 internal nodes ]:
## taxa are rows
## 
## $physeq_count
## phyloseq-class experiment-level object
## otu_table()   OTU Table:          [ 283 taxa and 59 samples ]:
## sample_data() Sample Data:        [ 59 samples by 36 sample variables ]:
## tax_table()   Taxonomy Table:     [ 283 taxa by 7 taxonomic ranks ]:
## phy_tree()    Phylogenetic Tree:  [ 283 tips and 282 internal nodes ]:
## taxa are rows
```

```r
save_object_name(metaphlan_ps)
```


```r
import_human_2phyloseq(pathway_DNA,
                       pathway_RNA,
                       type = "# Pathway",
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE) -> humann_pw_ps

humann_pw_ps

save_object_name(humann_pw_ps)
```



```r
import_human_2phyloseq(kegg_DNA,
                       kegg_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_kegg_ps

humann_kegg_ps

save_object_name(humann_kegg_ps)
```


```r
import_human_2phyloseq(l4c_DNA,
                       l4c_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_l4c_ps

humann_l4c_ps

save_object_name(humann_l4c_ps)
```



```r
import_human_2phyloseq(go_DNA,
                       go_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_go_ps

humann_go_ps

save_object_name(humann_go_ps)
```


```r
import_human_2phyloseq(infogo_DNA,
                       infogo_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_infogo_ps

humann_infogo_ps

save_object_name(humann_infogo_ps)
```


```r
import_human_2phyloseq(metac_DNA,
                       metac_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_metac_ps

humann_metac_ps

save_object_name(humann_metac_ps)
```


```r
import_human_2phyloseq(pfam_DNA,
                       pfam_RNA,
                       type = '# Gene Family',
                       n_rows = Inf,
                       str_rm_a = "_cat_Abundance-CPM_DNA",
                       str_rplace_a = "",
                       str_rm_b = "_cat_Abundance-CPM_RNA",
                       str_replace_b = "_RNA",
                       sample_column = "Sample_ID")  %>%
  subset_samples(Experiment == "Exp2") %>%
  physeq_count_normalize(return_both = TRUE)  -> humann_pfam_ps

humann_pfam_ps

save_object_name(humann_pfam_ps)
```

Would be great to apply this directly on the two lists
<https://stackoverflow.com/questions/19002378/applying-a-function-to-two-lists>

```r
# list("kegg_DNA" = kegg_DNA, 
#               "l4c_DNA" = l4c_DNA) -> list_DNA
# list("kegg_RNA" = kegg_RNA, 
#               "l4c_RNA" = l4c_RNA) -> list_RNA
# 
#   mapply(function(DNA_file_path,RNA_file_path) {
#   sapply(length(list_RNA), FUN = import_human_2phyloseq(DNA_file_path = list_DNA[[row]],RNA_file_path=list_RNA[[row]], type = '# Gene Family', n_rows = 10))
#   }, DNA_file_path=list_DNA, RNA_file_path=list_RNA) -> out
```

save RDS objects

```r
paste0(here::here(),
       "/data/processed/",
       "humann_metaphlan_import",
       "_",
       format(Sys.time(), "%Y%b%d")
       ,".RData") %>% save.image()
```



```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] speedyseq_0.5.3.9001 reshape2_1.4.4       scales_1.1.1        
##  [4] phyloseq_1.34.0      forcats_0.5.0        stringr_1.4.0       
##  [7] dplyr_1.0.3          purrr_0.3.4          readr_1.4.0         
## [10] tidyr_1.1.2          tibble_3.0.5         ggplot2_3.3.3       
## [13] tidyverse_1.3.0.9000
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-151        fs_1.5.0            HTSSIP_1.4.1       
##  [4] lubridate_1.7.9.2   progress_1.2.2      httr_1.4.2         
##  [7] rprojroot_2.0.2     tools_4.0.3         backports_1.2.1    
## [10] vegan_2.5-7         R6_2.5.0            DT_0.17            
## [13] mgcv_1.8-33         DBI_1.1.1           BiocGenerics_0.36.0
## [16] colorspace_2.0-0    permute_0.9-5       rhdf5filters_1.2.0 
## [19] ade4_1.7-16         withr_2.4.1         tidyselect_1.1.0   
## [22] prettyunits_1.1.1   compiler_4.0.3      cli_2.2.0          
## [25] rvest_0.3.6         Biobase_2.50.0      xml2_1.3.2         
## [28] digest_0.6.27       rmarkdown_2.6       XVector_0.30.0     
## [31] pkgconfig_2.0.3     htmltools_0.5.1.1   dbplyr_2.0.0       
## [34] htmlwidgets_1.5.3   rlang_0.4.10        readxl_1.3.1       
## [37] rstudioapi_0.13     generics_0.1.0      jsonlite_1.7.2     
## [40] crosstalk_1.1.1     magrittr_2.0.1      biomformat_1.18.0  
## [43] Matrix_1.3-2        Rcpp_1.0.6          munsell_0.5.0      
## [46] S4Vectors_0.28.1    Rhdf5lib_1.12.0     fansi_0.4.2        
## [49] ape_5.4-1           lifecycle_0.2.0     stringi_1.5.3      
## [52] yaml_2.2.1          MASS_7.3-53         zlibbioc_1.36.0    
## [55] rhdf5_2.34.0        plyr_1.8.6          grid_4.0.3         
## [58] parallel_4.0.3      crayon_1.3.4        lattice_0.20-41    
## [61] splines_4.0.3       Biostrings_2.58.0   haven_2.3.1        
## [64] multtest_2.46.0     hms_1.0.0           knitr_1.30         
## [67] pillar_1.4.7        igraph_1.2.6        codetools_0.2-18   
## [70] stats4_4.0.3        reprex_1.0.0        glue_1.4.2         
## [73] evaluate_0.14       data.table_1.13.6   modelr_0.1.8       
## [76] vctrs_0.3.6         foreach_1.5.1       cellranger_1.1.0   
## [79] gtable_0.3.0        assertthat_0.2.1    xfun_0.20          
## [82] broom_0.7.3         survival_3.2-7      iterators_1.0.13   
## [85] IRanges_2.24.1      cluster_2.1.0       ellipsis_0.3.1     
## [88] here_1.0.1
```
