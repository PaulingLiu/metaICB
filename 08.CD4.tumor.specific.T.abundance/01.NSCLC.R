
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)

#--- load data

cd4.sce <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/05.CD4.sct.new.rds.gz")
cd4.sce <- cd4.sce[,Idents(cd4.sce) != 'XCL1']
DimPlot(cd4.sce[,cd4.sce$cellid %in% tda$cellid[tda$cluster == "Th1"]], reduction = "umap", label = T) + scale_color_manual(values = comb.d3[c(1,2,3,5,6,7,10,12,13,15)])
FeaturePlot(cd4.sce[,cd4.sce$patient == "P36"], reduction = "umap", features = "CXCL13")

#--- match clones

all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")

tda <- tibble(
  cellid = cd4.sce$cellid,
  patient = cd4.sce$patient,
  sample = paste(cd4.sce$patient, cd4.sce$treatment, cd4.sce$num, sep = "."),
  cluster = as.character(Idents(cd4.sce))
) %>%
  dplyr::mutate(patient = stringr::str_remove(patient, "\\.0")) %>%
  dplyr::left_join(all.tcr[,c('cellid',"clone.id")], by = "cellid") %>%
  dplyr::filter(patient %in% c("P1","P10","P13","P19","P30","P33","P35","P36","P37","P38"))

th1.clones <- tda %>%
  dplyr::filter(cluster == "CD4_C7-Th1-like") %>%
  dplyr::pull(clone.id) %>%
  unique()

response = c('PR',"PD","PR","pre","PR","pre","PR","PD","pre","PR","pre","PR","pre","PR","pre","PR","pre","PD","PD","PD")

tda <- tda %>%
  dplyr::mutate(cluster = ifelse(clone.id %in% th1.clones & !is.na(clone.id), "CD4_C7-Th1-like", cluster)) %>%
  dplyr::mutate(cluster = ifelse(cluster == "CD4_C7-Th1-like", "Th1", "other"))

tda %>%
  dplyr::mutate(sample = stringr::str_remove(sample, "\\.0")) %>%
  dplyr::count(sample, cluster) %>%
  tidyr::spread(key = "cluster", value = "n") %>%
  dplyr::mutate(Th1 = ifelse(is.na(Th1), 0, Th1)) %>%
  dplyr::mutate(prop = Th1/(Th1+other)) %>%
  dplyr::mutate(response = response) %>%
  ggplot(aes(response, prop)) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(width = 0.1) +
  theme_classic()
