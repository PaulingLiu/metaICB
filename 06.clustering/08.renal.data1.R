
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

#---------- load data -----------#

tcr.data <- readr::read_tsv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/ccRCC_TCRs.txt")
tcell.sce <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/tcell.sce.rds.gz")
tcr.data <- tcr.data %>% dplyr::filter(cell %in% colnames(tcell.sce))

tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% tcr.data$cell]
all.tcr <- tcr.data %>% dplyr::filter(cell %in% colnames(tcell.filt)) %>% dplyr::rename(clone.id = TCR_clone)

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/08.renal.1.rds.gz")
pda <- pda %>% 
  dplyr::filter(clone == "Tex") %>% 
  dplyr::arrange(desc(size))

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% pda$clone.id)
sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cell]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cell
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]
sce.pro$trb <- all.tcr.pro$clone.id

#---------- samples ---------#

rep.data <- tibble(patient = sce.pro$Sample) %>%
  dplyr::mutate(response = ifelse(patient %in% c("t3","t4"), "PR", "PD"))


sce.pro$response <- rep.data$response
sce.pro$patient <- as.character(sce.pro$Sample2)


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
sce.pro.sub <- sce.pro[,sce.pro$patient %in% names(cell.count[cell.count>15])]
sce.pro.sub <- seurat.process(sce.pro.sub)
sce.pro.sub <- bbknn.batch(sce.pro.sub)

DimPlot(sce.pro.sub, label = T)

cdg <- readr::read_rds("/raid1/pauling/projects/01_data/03_gene/coding_gene.rds.gz")
cdg <- cdg %>% dplyr::select(gene_name, gene_id) %>% dplyr::distinct(gene_name, .keep_all = T)
cda <- as.data.frame(cdg)
rownames(cdg) <- cdg$gene_name

getgene <- function(.x){
  cdg[.x,] %>% dplyr::pull(gene_id)
}

FeaturePlot(sce.pro.sub[,sce.pro.sub$Sample == "t3"], features = getgene("ENTPD1"), pt.size = 0.5)

#---------- Annotation

clusters <- Idents(sce.pro.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(0,1,2,3,4,5,6)] <- "TD"
cluster.name[clusters %in% c(7)] <- "Prolif."

sce.pro.sub$clustesrs <- Idents(sce.pro.sub)
Idents(sce.pro.sub) <- cluster.name

sce.pro.sub %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/01.tumor.reactive.cells.batch1.rds.gz", compress = "gz")

DimPlot(sce.pro.sub)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/01.tumor.reactive.cells.batch1.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  patient = as.character(sce.pro$Sample)
) %>%
  dplyr::filter(patient == "t3") %>%
  dplyr::arrange(cluster) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 1.5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#F08080","grey90","#008B8B","#3CB371"))

#----------- prop analysis

tmp <- table(paste(sce.pro$response, sce.pro$patient, sep = "."), Idents(sce.pro))
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
