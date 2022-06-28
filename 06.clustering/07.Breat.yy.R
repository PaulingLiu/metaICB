
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

#---------- load data -----------#

tcell.sce <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/tcell.rds.gz")
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

tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellID]
all.tcr <- all.tcr %>% dplyr::filter(cellID %in% colnames(tcell.filt))

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/07.breast.yy.rds.gz")
pda <- pda %>% 
  dplyr::filter(clone == "Tex") %>% 
  dplyr::arrange(desc(size))

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% pda$clone.id)
sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cellID]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellID
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]
sce.pro$trb <- all.tcr.pro$clone.id

#---------- samples ---------#

rep.data <- tibble(cellID = colnames(sce.pro)) %>%
  dplyr::left_join(all.tcr, by = "cellID")

sce.pro$response <- rep.data$Efficacy
sce.pro$patient <- rep.data$Patient
sce.pro$treatment <- rep.data$Group

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
sce.pro.sub <- bbknn.batch(sce.pro, r = 1.5)

DimPlot(sce.pro.sub[,sce.pro.sub$patient == "P019" & sce.pro.sub$treatment == "Post-treatment"], label = T)
FeaturePlot(sce.pro[,sce.pro$patient == "P019" & sce.pro$treatment == "Post-treatment"], features = c("GZMK"))

#---------- Annotation

clusters <- Idents(sce.pro)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(2,0,4)] <- "0.GZMK+Texp"
cluster.name[clusters %in% c(1,6,8)] <- "0.IL7R+Texp"
cluster.name[clusters %in% c(5,7,3)] <- "TD"
#cluster.name[Idents(sce.pro.sub) == 15] <- "Prolif."

sce.pro$clustesrs <- Idents(sce.pro)
Idents(sce.pro) <- cluster.name

sce.pro %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/01.tumor.reactive.cells.batch1.rds.gz", compress = "gz")

DimPlot(sce.pro.sub)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- clustering new

markers <- FindAllMarkers(sce.pro, only.pos = T)
sce.pro <- FindVariableFeatures(sce.pro, nfeatures = 2000)
sce.pro <- RunPCA(sce.pro, features = unique(markers$gene))
#sce.pro <- RunUMAP(sce.pro, dims = 1:15)

pca <- sce.pro@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca[,1:(20-1)], obs=sce.pro$patient)
sc$tl$pca(adata)
adata$obsm$X_pca = pca[,1:(20-1)]
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 1.4)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce.sub <- loadumap(sce.pro, umap, bb.res)
DimPlot(sce.sub, label = T)
#---------- Annotation

clusters <- Idents(sce.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(1,0,6,10,11,4)] <- "0.GZMK+Texp"
cluster.name[clusters %in% c(2,8)] <- "0.IL7R+Texp"
cluster.name[clusters %in% c(3,5,9,7,12)] <- "TD"
cluster.name[clusters %in% c(13)] <- "Prolif."

sce.sub$clustesrs <- Idents(sce.sub)
Idents(sce.sub) <- cluster.name

sce.sub %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/02.tumor.reactive.cells.batch1.rds.gz", compress = "gz")

DimPlot(sce.sub[,sce.sub$patient == "P012" & sce.sub$treatment == "Pre-treatment"], label = T)


#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/02.tumor.reactive.cells.batch1.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  response = sce.pro$response,
  treatment = sce.pro$treatment
) %>%
  dplyr::arrange(desc(cluster)) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 2, alpha = 0.7) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = comb.d3[c(10,13,7,12)])

FeaturePlot(sce.pro, features = "MKI67", cols = c("grey92", "#EE0000", "#CD2626"), pt.size = 1)

#----------- prop analysis
sce.sub <- breast3.sce
tmp <- table(paste(sce.sub$response, sce.sub$treatment, sce.sub$patient, sep = "."), Idents(sce.sub))
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 19,]

p7 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(sample = stringr::str_remove(sample,"-treatment")) %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  ggplot(aes(factor(response, levels = c("PR.Pre","PR.Post","SD.Pre","SD.Post")), `0.GZMK+Texp`+`0.IL7R+Texp`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response),  width = 0, size = 2) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of precursor cells in all\nCXCL13+ tumor-reactive CD8 T cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.y = element_text(size = 0, colour = "black")) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD")) +
  ylim(0,1)

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","tumor"), sep = "\\.") %>%
  dplyr::select(response, treatment,`0.GZMK+Texp`,`0.IL7R+Texp`)

a$prop <- rowSums(a[,-c(1,2)])
wilcox.test(a$prop[1:2], a$prop[6:7])
