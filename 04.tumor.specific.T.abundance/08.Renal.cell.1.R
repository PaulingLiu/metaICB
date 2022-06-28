
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

tcr.data <- readr::read_tsv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/ccRCC_TCRs.txt")
sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/ccRCC_6pat_Seurat")
sce <- UpdateSeuratObject(sce)
table(sce$region, sce$type)
sce <- sce[,sce$type == "Tumor"]
table(sce$Sample)

sce <- sce[,sce$Sample %in% c("t2","t3","t4")]


#---------- T cell ------------#

sce.sub <- sce[,colnames(sce) %in% tcr.data$cell]
sce.sub <- NormalizeData(sce.sub)
sce.sub <- ScaleData(sce.sub, do.center = F, do.scale = F)

sce.sub <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/tcell.sce.rds.gz")

tcr.data <- tcr.data %>% dplyr::filter(cell %in% colnames(sce.sub))

#---------- Tex clones ----------#

all.tcr <- tibble(
  cellID = colnames(sce.sub),
  Patient = as.character(sce.sub$Sample),
  sample = sce.sub$Sample2
) %>%
  dplyr::inner_join(tcr.data, by = c("cellID" = "cell")) %>%
  dplyr::rename(clone.id = TCR_clone) 

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::group_by(Patient, sample, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cdg <- readr::read_rds("/raid1/pauling/projects/01_data/03_gene/coding_gene.rds.gz")
cdg %>%
  dplyr::filter(gene_name %in% c("CXCL13","CD8A"))

cal.genes <- c("ENSG00000156234","ENSG00000153563")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellID
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- sce.sub@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(sce.sub@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"ENSG00000153563"],
  CXCL13 = pr.matr[,"ENSG00000156234"],
  size = use.clones.cells$num,
  patient = use.clones.cells$Patient,
  sample = use.clones.cells$sample,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/08.renal.1.rds.gz")

#------------- cutoff determination -------------#

tibble(CD8A = tex.clones.pro$CD8A) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
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

tibble(a = tex.clones.pro$CXCL13) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
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

tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 1 & CXCL13 > 0.8, "Tex", "Other"))

tex.clones.pro %>%
  dplyr::group_by(patient, sample, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Tex = ifelse(is.na(Tex), 0 , Tex)) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

mda <- tibble(
  patient = c("t2","t3","t4"),
  response = c("non.rp","rp","rp")
) %>%
  dplyr::distinct(patient, response)

pro.pda <- pro.pda %>%
  dplyr::filter(Other+Tex > 49) %>%
  dplyr::left_join(mda, by = "patient") %>%
  dplyr::mutate(resp = response)


pro.pda %>%
  ggplot(aes(factor(response, levels = c("rp","non.rp")), freq)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response, shape = patient), size = 2, stroke = 1.2,width = 0.15) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#009ACD","#CD2626")) +
  scale_shape_manual(values = c(1,2,3))

uniq.resp <- unique(pro.pda$resp)
wilcox.test(pro.pda$freq[pro.pda$resp == uniq.resp[1]], pro.pda$freq[pro.pda$resp == uniq.resp[2]])

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
