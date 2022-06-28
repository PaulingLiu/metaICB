
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/08.lung.pd1/02.all.rds.gz")
all.tcr  <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/08.lung.pd1/03.all.tcr.rds.gz")

sce.filt <- sce[,colnames(sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(sce.filt))

sce.filt <- NormalizeData(sce.filt, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sce.filt)
sce.filt <- ScaleData(sce.filt, features = all.genes, do.center = F, do.scale = F)

#---------- Tex clones ----------#

all.tcr.pro <- all.tcr %>%
  tidyr::separate(cellid, c("patient","tumor","num"), sep = "[_\\.]", remove = F) %>%
  #dplyr::mutate(clone.id = paste(patient, tumor,num,`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ".")) %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ".")) %>%
  as.tibble()

use.clones <- unique(all.tcr.pro$clone.id)

use.clones.cells <- all.tcr.pro %>%
  dplyr::mutate(sample = paste(patient, tumor, num, sep = ".")) %>%
  dplyr::select(sample, patient, clone.id, cellid) %>%
  dplyr::group_by(sample, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD8A","CD4")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- sce.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(sce.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"CD8A"],
  CXCL13 = pr.matr[,"CXCL13"],
  CD4 = pr.matr[,"CD4"],
  size = use.clones.cells$num,
  #patient = use.clones.cells$patient,
  patient = use.clones.cells$sample,
  clone.id = use.clones.cells$clone.id,
  #clones = use.clones.cells$clones
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/02.NSCLC.2.rds.gz")

#------------- cutoff determination -------------#

tibble(CD8A = pda$CD4) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_density(lwd = 1, color = "#EE6363") +
  theme_classic() +
  labs(
    x = "CD8A expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

tibble(a = pda$CXCL13) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 0.4, linetype = "dashed") +
  geom_density(lwd = 1, color = "#008B8B") +
  theme_classic() +
  labs(
    x = "CXCL13 expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

#------------- Tex cell frequency -------------#

Tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id) %>%
  unique()

tex.clones.pro <- pda %>%
  dplyr::filter(CD4 > 0.25) %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  #dplyr::mutate(clone = ifelse(clone.id %in% Tex.clones, "Tex", "Other"))
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5, "Tex", "Other"))


rep.data <- tibble(
  patient = c("MD01-004","MD01-005","MD01-010","MD01-019","MD01-024","MD043-003","MD043-006",
              "MD043-008","MD043-011","NY016-007","NY016-014","NY016-015","NY016-021","NY016-022","NY016-025"),
  response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR',
               'non-MPR','non-MPR','non-MPR','MPR','MPR')
)

tex.clones.pro %>%
  tidyr::separate(patient, c("patient","tumor","num"), sep = "\\.") %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(rep.data, by = "patient") -> pda.pro

pda.pro %>%
  dplyr::mutate(residual = c(40,0,5,95,100,5,50,10,75,60,95,80,100,5,0)) %>%
  ggplot(aes(response, freq)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response), width = 0.15, size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "CD4"
  ) +
  scale_color_manual(values = c("#CD2626","#009ACD"))



wilcox.test(pda.pro$freq[pda.pro$response == "MPR"], pda.pro$freq[pda.pro$response != "MPR"])

#-------------- Classification -----------------#

Tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 0.25 & CXCL13 > 0.2, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id) %>%
  unique()

tex.clones.pro <- pda %>%
  dplyr::filter(CD8A > 0.25) %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  #dplyr::mutate(clone = ifelse(clone.id %in% Tex.clones, "Tex", "Other"))
  dplyr::mutate(clone = ifelse(CD8A > 0.25 & CXCL13 > 0.2, "Tex", "Other"))


rep.data <- tibble(
  patient = c("MD01-004","MD01-005","MD01-010","MD01-019","MD01-024","MD043-003","MD043-006",
              "MD043-008","MD043-011","NY016-007","NY016-014","NY016-015","NY016-021","NY016-022","NY016-025"),
  response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR',
               'non-MPR','non-MPR','non-MPR','MPR','MPR')
)

tex.clones.pro %>%
  tidyr::separate(patient, c("patient","tumor","num"), sep = "\\.") %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(rep.data, by = "patient") -> pda.pro.cd8

#--- classification

pda.pro.cd8 %>%
  dplyr::mutate(cd4 = pda.pro$freq) %>%
  dplyr::filter(response != "PR") %>%
  ggplot(aes(cd4, freq)) +
  geom_point(aes(color = response), size = 3) +
  #scale_x_sqrt() +
  #scale_y_sqrt() +
  theme_bw() +
  scale_color_manual(values = c("#B452CD", "#008B8B")) +
  geom_hline(yintercept = 0.55, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed")  +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "",
    y = ""
  )
