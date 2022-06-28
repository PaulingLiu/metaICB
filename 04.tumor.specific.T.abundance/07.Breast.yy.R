
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

all.sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/PDL1_PTX_immune.final.sct.rds")
meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/Stable2-new-CCR3.csv")
tda <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/TCRfile-new-CCR3.csv")

meta <- meta %>%
  dplyr::filter(`Cell barcode` %in% colnames(tcell.sce)) %>%
  dplyr::filter(Treatment == "anti-PDL1+Chemo") %>%
  dplyr::filter(`Major celltype` %in% c("ILC cell", "T cell")) %>%
  dplyr::filter(`Cell barcode` %in% tda$cellID)

tda <- tda %>% dplyr::filter(cellID %in% meta$`Cell barcode`)
all.tcr <- tda[,c("cellID","clone.id")] %>%
  dplyr::inner_join(meta[,c("Cell barcode","Patient","Sample","Origin","Group","Efficacy")], by = c("cellID" = "Cell barcode")) %>%
  dplyr::filter(Origin != "blood") %>%
  dplyr::filter(Group != "Progression")

#--------- tcell object ----------#

all.sce <- all.sce[,colnames(all.sce) %in% tda$cellID]
all.sce@assays$SCT <- NA

count.data <- all.sce@assays$RNA@counts
tcell.sce <- CreateSeuratObject(counts = count.data)
tcell.sce <- NormalizeData(tcell.sce)
tcell.sce <- ScaleData(tcell.sce, do.center = F, do.scale = F)



#---------- Tex clones ----------#

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::select(Patient, Group, clone.id, cellID) %>%
  dplyr::group_by(Patient, Group, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD8A")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellID
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- tcell.sce@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(tcell.sce@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD8A = pr.matr[,"CD8A"],
  CXCL13 = pr.matr[,"CXCL13"],
  size = use.clones.cells$num,
  patient = use.clones.cells$Patient,
  group = use.clones.cells$Group,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/07.breast.yy.rds.gz")

#------------- cutoff determination -------------#

tibble(CD8A = tcell.sce@assays$RNA@scale.data["CD8A",]) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
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


tibble(a = tcell.sce@assays$RNA@scale.data["CXCL13",]) %>%
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

tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 0.5 & CXCL13 > 0.4, "Tex", "Other"))

tex.clones.pro %>%
  dplyr::group_by(patient, group, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Tex = ifelse(is.na(Tex), 0 , Tex)) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

mda <- tibble(
  patient = meta$Patient,
  response = meta$Efficacy
) %>%
  dplyr::distinct(patient, response)

pro.pda <- pro.pda %>%
  dplyr::filter(Other+Tex > 49) %>%
  dplyr::left_join(mda, by = "patient") %>%
  dplyr::mutate(resp = paste0(response,".",group))


pro.pda %>%
  ggplot(aes(factor(resp, levels = c("PR.Pre-treatment","PR.Post-treatment","SD.Pre-treatment","SD.Post-treatment")), freq)) +
  geom_boxplot(aes(color = resp), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = resp), size = 2, width = 0.15) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD"))


uniq.resp <- unique(pro.pda$resp)
wilcox.test(pro.pda$freq[pro.pda$resp == uniq.resp[4]], pro.pda$freq[pro.pda$resp == uniq.resp[1]])

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
