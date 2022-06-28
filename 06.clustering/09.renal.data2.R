
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

#---------- load data -----------#

sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/03.RCC.data.2/01.tcell.rds.gz")
tibble(
  CXCL13 = sce@assays$RNA@scale.data["CXCL13",],
  patient = sce$patient,
  clusters = Idents(sce),
  cells = colnames(sce)
) %>%
  dplyr::mutate(type = ifelse(clusters %in% c(0,9), "CD8", "CD4")) -> pda

Idents(sce) <- pda$type
FeaturePlot(sce[,sce$patient %in% c("p55")], features = c("GZMK","CXCL13","HAVCR2","LAG3"))

#---------- Tex clones ----------#

pda <- pda %>%
  dplyr::mutate(clone = ifelse(type == "CD8" & CXCL13 > 0.8, "Tex", "Other"))

sce.tex <- sce[,pda$cells[pda$clone == "Tex"]]
sce.pro <- sce.tex[,sce.tex$patient %in% c("p55","p915","p913","p906")]


 c("p55","p915","p913","p906")
c("rp","rp","non.rp","non.rp")


#----------- batch correction

seurat.process <- function(.x){
  .x <- NormalizeData(.x)
  .x <- ScaleData(.x, do.scale = F, do.center = F)
  .x <- FindVariableFeatures(.x, nfeatures = 2000)
  .x <- RunPCA(.x, features = VariableFeatures(.x))
  .x <- RunUMAP(.x, dims = 1:15)
  return(.x)
}

cell.count <- table(sce.pro$patient)
sce.pro.sub <- sce.pro[,sce.pro$patient %in% names(cell.count[cell.count>29])]
sce.pro.sub <- seurat.process(sce.pro.sub)
sce.pro.sub <- bbknn.batch(sce.pro.sub)

DimPlot(sce.pro.sub, label = T)
FeaturePlot(sce.pro.sub[,sce.pro.sub$patient == "p915"], features = c("LAG3"))

#---------- Annotation

clusters <- Idents(sce.pro.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(0,1,2,3,4,6,7)] <- "TD"
cluster.name[clusters %in% c(5)] <- "Prolif."

sce.pro.sub$clustesrs <- Idents(sce.pro.sub)
Idents(sce.pro.sub) <- cluster.name

sce.pro.sub %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/03.RCC.data.2/01.tumor.reactive.cells.batch1.rds.gz", compress = "gz")

DimPlot(sce.pro.sub)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/01.tumor.reactive.cells.batch1.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro))
) %>%
  dplyr::arrange(cluster) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 1.5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#F08080","grey90","#008B8B","#3CB371"))

#----------- prop analysis

tmp <- table(paste(sce.pro$response, sce.pro$treatment, sce.pro$patient, sep = "."), Idents(sce.pro))
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 19,]

as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(sample = stringr::str_remove(sample,"-treatment")) %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  ggplot(aes(factor(response, levels = c("PR.Pre","PR.Post","SD.Pre","SD.Post")), `0.GZMK+Texp`+`0.IL7R+Texp`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response),  width = 0) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of precursor cells in all\nCXCL13+ tumor-reactive CD8 T cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD")) +
  ylim(0,1)

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","tumor"), sep = "\\.") %>%
  dplyr::select(response, treatment,`0.GZMK+Texp`)

#a$prop <- rowSums(a[,-1])
t.test(a$`0.GZMK+Texp`[1:8], a$`0.GZMK+Texp`[17:19])
