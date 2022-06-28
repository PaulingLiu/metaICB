
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

all.tcr <- all.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::filter(cellid %in% cd4.sce$cellid) %>%
  dplyr::select(patient, clone.id, cellid) %>%
  dplyr::group_by(patient, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD4","CD8A")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- cd4.sce@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(cd4.sce@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD4 = pr.matr[,"CD4"],
  CXCL13 = pr.matr[,"CXCL13"],
  CD8A = pr.matr[,"CD8A"],
  size = use.clones.cells$num,
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/01.NSCLC.1.rds.gz")
#------------- cutoff determination -------------#

tibble(CD4 = pda$CD8A) %>%
  ggplot(aes(CD4)) +
  geom_vline(xintercept = 0.6, linetype = "dashed") +
  geom_density(lwd = 1, color = "#EE6363") +
  theme_classic() +
  labs(
    x = "CD4 expression",
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

#--- boxplot

tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.4, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- pda %>%
  #dplyr::filter(CD4 > 0.25) %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(clone.id %in% tex.clones, "Tex", "Other"))

response = c('PR',"PD","PR","pre","PR","pre","PR","PD","pre","PR","pre","PR","pre","PR","pre","PR","pre","PD","PD","PD")
used.color <- c("#FA8072", "#6CA6CD", "grey50")

tex.clones.pro %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  #.[c(1:6,9:11,17:18,32:35,37:41),] %>%
  #dplyr::mutate(response = response) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

pro.pda %>%
  ggplot(aes(factor(response, levels = c("pre","PR","PD")), freq)) +
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
  scale_color_manual(values = c("#009ACD", "#CD2626", "#EE799F"))


#--- classification

pro.pda.cd8 %>%
  dplyr::mutate(cd4 = pro.pda$freq[pro.pda$patient %in% pro.pda.cd8$patient]) %>%
  dplyr::filter(stringr::str_detect(patient, "ut")) %>%
  #dplyr::filter(response != "PR") %>%
  ggplot(aes(cd4,freq)) +
  geom_point(size = 3, alpha = 0.8, shape = 21, stroke = 0.8) +
  geom_text_repel(aes(label = patient), box.padding = 0.5, max.overlaps = Inf) +
  #geom_label(aes(label = patient), nudge_x = 0.03, parse = T) +
  scale_x_sqrt() +
  scale_y_sqrt() +
  theme_bw() +
  scale_color_manual(values = c("black", "#B452CD")) +
  scale_fill_manual(values = c("#008B8B", "#B452CD")) +
  geom_hline(yintercept = 0.06, linetype = "dashed") +
  geom_vline(xintercept = 0.01, linetype = "dashed")  +
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
