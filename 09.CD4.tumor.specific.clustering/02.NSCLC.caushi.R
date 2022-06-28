
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
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ".")) %>%
  as.tibble()

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/02.NSCLC.2.rds.gz")
Tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone.id.new = paste(patient,clone.id, sep = ".")) %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5 & CD8A == 0, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 3) %>%
  dplyr::pull(clone.id.new) %>%
  unique()

all.tcr.pro <- all.tcr.pro %>% 
  dplyr::mutate(clone.id.new = paste(patient, tumor,num,clone.id, sep = ".")) %>%
  dplyr::filter(clone.id.new %in% Tex.clones)

tda <- data.frame(
  cellid = all.tcr.pro$cellid,
  clone.id = all.tcr.pro$clone.id.new,
  row.names = all.tcr.pro$cellid
)

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
sce.pro <- RunPCA(sce.pro, features = VariableFeatures(sce.pro))
sce.pro <- RunUMAP(sce.pro, dims = 1:15)
sce.sub <- sce.pro[,sce.pro$sample %in% unique(sce.pro$sample)[table(sce.pro$sample) > 50]]
tda <- tda[sce.sub$cellid,]
sce.sub$clone <- tda$clone.id

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
DimPlot(sce.sub, label = T)
markers <- FindMarkers(sce.sub, ident.1 = 3, only.pos = T)
FeaturePlot(sce.sub, features = c("TCF7","CXCL13","BCL6","PDCD1","IL7R","HAVCR2"), ncol = 3)


a <- sce.sub[,Idents(sce.sub) == 1]
rank.clones <- table(a$clone) %>%
  as.data.frame() %>%
  as.tibble() %>%
  dplyr::arrange(desc(Freq))


sce.sub.new <- sce.sub[,!(Idents(sce.sub) %in% c(1,9))]
FeaturePlot(sce.sub.new, features = "HAVCR2")

sce.sub.new %>% readr::write_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/10.CD4.tumor.reactive.T.clustering/02.NSCLC.caushi.rds.gz", compress = "gz")
