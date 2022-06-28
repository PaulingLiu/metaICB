
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

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
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ".")) %>%
  as.tibble()

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/02.NSCLC2.rds.gz")
pda <- pda %>% dplyr::filter(clone == "Tex") %>% dplyr::filter(size > 5)

all.tcr.pro <- all.tcr.pro %>% dplyr::filter(clone.id %in% pda$clone.id)

sce.pro <- sce.filt[,colnames(sce.filt) %in% all.tcr.pro$cellid]

#---------- samples ---------#

rep.data <- tibble(
  patient = c("MD01-004","MD01-005","MD01-010","MD01-019","MD01-024","MD043-003","MD043-006",
              "MD043-008","MD043-011","NY016-007","NY016-014","NY016-015","NY016-021","NY016-022","NY016-025"),
  response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR',
               'non-MPR','non-MPR','non-MPR','MPR','MPR')
)

mda <- tibble(cells = colnames(sce.pro)) %>%
  tidyr::separate(cells, c("sample","barcode"), sep = "\\.", remove = F) %>%
  tidyr::separate(sample, c("patient","tumor","num"), sep = "\\_", remove = F) %>%
  dplyr::mutate(tumor = paste(patient, tumor, sep = "_")) %>%
  dplyr::left_join(rep.data, by = "patient")
  
sce.pro$sample <- mda$sample
sce.pro$tumor <- mda$tumor
sce.pro$response <- mda$response

#----------- batch correction

sce.pro <- FindVariableFeatures(sce.pro, nfeatures = 2000)
sce.pro <- RunPCA(sce.pro, features = unique(markers$gene))
sce.pro <- RunUMAP(sce.pro, dims = 1:15)
sce.sub <- sce.pro[,sce.pro$sample %in% unique(sce.pro$sample)[table(sce.pro$sample) > 100]]

pca <- sce.sub@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca[,1:(20-1)], obs=sce.sub$sample)
sc$tl$pca(adata)
adata$obsm$X_pca = pca[,1:(20-1)]
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 1.5)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce.sub <- loadumap(sce.sub, umap, bb.res)

#sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/02.data/01.lung.nature/04.cxcl13.sce.rds.gz")

#---------- Annotation

cluster.name[clusters %in% c(1)] <- "GZMK+Texp"
cluster.name[clusters %in% c(0,9,14)] <- "IL7R+Texp"
cluster.name[clusters %in% c(6,13)] <- "Prolif."
cluster.name[clusters %in% c(12,4,11,8)] <- "HSP+ TD"
cluster.name[cluster.name == "cells"] <- "TD"

sce.sub$clustesrs <- Idents(sce.sub)
Idents(sce.sub) <- cluster.name
sce.sub <- readr::read_rds("projects/07_pan.cancer.treatment/02.data/01.lung.nature/04.sample.batch.cxcl13.sce.rds.gz")

DimPlot(sce.sub, label = T)
FeaturePlot(sce.sub, features = "GZMK", pt.size = 0.2, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()
VlnPlot(sce.pro.new, features = "HAVCR2", pt.size = 0)

#---------- UMAP plot

tibble(
  umap1 = lung2.sce@reductions$umap@cell.embeddings[,1],
  umap2 = lung2.sce@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(lung2.sce)),
  response = lung2.sce$response
) %>%
  dplyr::arrange(desc(cluster)) %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 0.8, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  #scale_color_manual(values = comb.d3[c(10,13,7,12,12)]) +
  scale_color_manual(values = comb.d3[c(10,12,13,7,12)])

tibble(
  umap1 = sce.sub@reductions$umap@cell.embeddings[,1],
  umap2 = sce.sub@reductions$umap@cell.embeddings[,2],
  cluster = Idents(sce.sub),
  response = sce.sub$response
) %>%
  dplyr::arrange(cluster) %>%
  dplyr::filter(response != "MPR") %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 0.6) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#3CB371", "grey90","grey90",comb.d3[10], "grey90"))

#----------- Scibet classification
sce.sub <- lung2.sce
matr <- sce.sub@assays$RNA@counts
matr <- as.matrix(matr)
matr <- t(matr)
matr <- as.data.frame(matr)
clusters <- Idents(sce.sub)
clusters <- as.character(clusters)
matr$label <- clusters

#------- query set -------#

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/02.data/01.lung.nature/04.cxcl13.sce.rds.gz")
q.expr <- sce.pro@assays$RNA@counts
q.expr <- t(as.matrix(q.expr))
q.expr <- as.data.frame(q.expr)
prd <- scibet::SciBet(matr, q.expr, k = 2000)

#----------- prop analysis

tmp <- table(paste(sce.pro$response, sce.pro$tumor, sep = "."), prd)
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 49,]
as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  ggplot(aes(response, `IL7R+Texp`+`GZMK+Texp`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response), width = 0.1, size = 2) +
  theme_classic() +
  labs(
    x = "",
    y = "Proportion of precursor cells in all\nCXCL13+ tumor-reactive CD8 T cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black")) +
  scale_color_manual(values = c("#CD2626", "#009ACD")) +
  ylim(0,1)

as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::arrange(tumor) %>%
  dplyr::mutate(ratio = c(40,0,5,95,100,5,10,75,75,60,95,80,100,5,0)) %>%
  ggplot(aes(ratio, `IL7R+Texp`+`GZMK+Texp`)) +
  #geom_vline(xintercept = 10, linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se =F, color = "grey10", lwd = 0.6) +
  theme_classic() +
  labs(
    x = "% residual tumor",
    y = "Proportion of precursor cells in all\nCXCL13+ tumor-reactive CD8 T cells"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 10, colour = "black")) +
  scale_color_manual(values = c("#CD2626", "#009ACD")) +
  ylim(0,1)

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::select(response, `IL7R+Texp`,`GZMK+Texp`)

a$prop <- rowSums(a[,-1])
t.test(a$prop[1:6], a$prop[7:15])

cda <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::arrange(tumor) %>%
  dplyr::mutate(ratio = c(40,0,5,95,100,5,10,75,75,60,95,80,100,5,0))

cor.test(cda$`IL7R+Texp`+cda$`GZMK+Texp`, cda$ratio)
