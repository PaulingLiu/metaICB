
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

tcell.sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1864-counts_tcell_cohort1.rds")
tcell.sce <- CreateSeuratObject(counts = tcell.sce)
tcell.sce <- NormalizeData(tcell.sce)
tcell.sce <- ScaleData(tcell.sce, do.center = F, do.scale = F)

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1863-counts_cells_cohort1-T_cell.metadata.csv", col_names = F)
meta <- meta %>% 
  dplyr::select(X1, X4, X5, X6, X11) %>%
  dplyr::rename(
    cellid = X1,
    patient = X4,
    treatment = X5,
    response = X6,
    batch = X11
  ) %>%
  dplyr::filter(response %in% c("E","NE"))

tcell.sce <- tcell.sce[, colnames(tcell.sce) %in% meta$cellid]

#---------- TCR data -----------#

all.tcr1 <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1879-BIOKEY_barcodes_vdj_combined_cohort1.xls")
all.tcr2 <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1880-BIOKEY_barcodes_vdj_combined_cohort2.xls")
all.tcr <- all.tcr1 %>% dplyr::bind_rows(all.tcr2)

all.tcr <- all.tcr %>% 
  dplyr::rename(cellid = barcode, clone.id = clonotype, sample = SAMPLE_ID) %>%
  dplyr::mutate(cellid = paste0(cellid, "-1")) %>%
  dplyr::filter(cellid %in% colnames(tcell.sce))

#---------- Tex clones ----------#

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::select(sample, clone.id, cellid) %>%
  dplyr::group_by(sample, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD4", "CD8A")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
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
  CD4 = pr.matr[,"CD4"],
  CXCL13 = pr.matr[,"CXCL13"],
  size = use.clones.cells$num,
  sample = use.clones.cells$sample,
  clone.id = use.clones.cells$clone.id
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/05.BC.batch1.rds.gz")
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
  geom_vline(xintercept = 0.5, linetype = "dashed") +
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
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5, "Tex", "Other"))

#CD8 0.5

tex.clones.pro %>%
  dplyr::filter(CD4 > 0.25) %>%
  dplyr::group_by(sample, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Tex = ifelse(is.na(Tex), 0 , Tex)) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

mda <- tibble(
  patient = meta$patient,
  response = meta$response
) %>%
  dplyr::distinct(patient, response)

pro.pda <- pro.pda %>%
  tidyr::separate(sample, c("bio","num","tr")) %>%
  dplyr::mutate(patient = paste0(bio,"_",num)) %>%
  #dplyr::filter(Other+Tex > 49) %>%
  dplyr::left_join(mda, by = "patient") %>%
  dplyr::mutate(resp = paste0(response,".",tr))


pro.pda %>%
  #dplyr::filter(batch == 2) %>%
  ggplot(aes(factor(resp, levels = c("E.Pre","E.On","NE.Pre","NE.On")), freq)) +
  geom_boxplot(aes(color = resp), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = resp), size = 2, width = 0.15) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "CD4"
  ) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD"))


uniq.resp <- unique(pro.pda$resp)
t.test(pro.pda$freq[pro.pda$resp == uniq.resp[1]], pro.pda$freq[pro.pda$resp == uniq.resp[3]])

#-------------- Classification

tex.clones.pro <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 0.5 & CXCL13 > 0.5, "Tex", "Other"))

#CD8 0.5

tex.clones.pro %>%
  dplyr::filter(CD8A > 0.5) %>%
  dplyr::group_by(sample, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Tex = ifelse(is.na(Tex), 0 , Tex)) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda.cd8


pro.pda %>%
  dplyr::mutate(cd4 = pro.pda.cd8$freq) %>%
  dplyr::filter(resp %in% c("E.Pre","NE.Pre","NE.On")) %>%
  ggplot(aes(freq, cd4)) +
  geom_point(aes(fill = resp, color = resp), size = 3, alpha = 0.8, shape = 21, stroke = 0.8) +
  theme_bw() +
  scale_color_manual(values = c("#B452CD","black","#008B8B")) +
  scale_fill_manual(values = c("#B452CD","#008B8B","#008B8B")) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_vline(xintercept = 0.23, linetype = "dashed") +
  scale_y_sqrt() +
  scale_x_sqrt() +
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
