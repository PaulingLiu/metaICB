
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

tcell.sce <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/02.tcell.sce.rds.gz")
all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")

tcell.sce <- tcell.sce[,tcell.sce$patient %in% c("P1","P10","P13","P19","P30","P33","P35","P36","P37","P38")]
tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))

#---------- Tex clones ----------#

all.tcr <- all.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::select(patient, clone.id, cellid) %>%
  dplyr::group_by(patient, clone.id) %>%
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
    meta.expr[i,] <- tcell.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(tcell.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
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
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/02.NSCLC2.new.0523.rds.gz")
pda %>% readr::write_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/02.NSCLC2.new.0523.rds.gz", compress = "gz")
#------------- cutoff determination -------------#

tibble(CD8A = pda$CD8A) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.6, linetype = "dashed") +
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
  geom_vline(xintercept = 1, linetype = "dashed") +
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

tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 0.6 & CXCL13 > 1, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- pda %>%
  #dplyr::filter(CD8A > 0.6) %>%
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
  #dplyr::mutate(response = response) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda.cd8

pro.pda.cd8 %>%
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
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#009ACD", "#CD2626", "#EE799F"))


wilcox.test(pro.pda$freq[pro.pda$response == "pre"], pro.pda$freq[pro.pda$response == "PD"])

#-------------- barplot -----------------#

pro.pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(freq), n = length(freq), sd = sd(freq)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("pre","PR","PD")), mean)) +
  geom_col(aes(fill = response, color = response),  alpha = 0.75, lwd = 0.7) +
  geom_jitter(aes(factor(response, levels = c("pre","PR","PD")), freq, color = response), size = 2.5, data = pro.pda, width = 0.1) +
  geom_errorbar(aes(ymin = mean, ymax = ymax), width = 0.25, lwd = 0.8) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion in T cells"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_fill_manual(values = c("#008B8B", "#4F94CD","#F08080")) +
  scale_color_manual(values = c("#008B8B", "#4F94CD","#F08080"))
