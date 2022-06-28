
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(reticulate)

#---------- data process -----------#

in.path <- "projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/01.geo/01.all.data"
out.path <- "/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/01.geo/02.samples"
files <- list.files(path = in.path)

fda <- tibble(files = files) %>%
  tidyr::separate(files, c("GEO","patient","antigen","file"), sep = "\\_", remove = F) %>%
  dplyr::mutate(sample = paste(GEO,patient, antigen,sep = "."))

for (i in 1:nrow(fda)) {
  if(i%%4 == 1){
    dir.create(file.path(out.path, fda$sample[i]))
  }
  
  if(i%%4 == 0){
    file.copy(from = file.path(in.path, fda$files[i]), to = file.path(out.path, fda$sample[i-1], fda$files[i]))
  }else{
    file.copy(from = file.path(in.path, fda$files[i]), to = file.path(out.path, fda$sample[i], fda$files[i]))
  }
}

#----------- TCR data -------------#

tcr.path <- "/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/01.geo/01.all.data"
tcr.files <- files[stringr::str_detect(files, "contig")]
tcr <- list()


for (i in 2:length(tcr.files)) {
  print(i)
  tda <- readr::read_csv(file.path(tcr.path, tcr.files[i]))
  tda.pro <- process.config(tda)
  tda.pro <- process.tcr(tda.pro)
  tcr[[i]] <- tda.pro %>% 
    dplyr::filter(!is.na(`CDR3(Alpha1)`)) %>%
    dplyr::filter(!is.na(`CDR3(Beta1)`)) %>%
    dplyr::mutate(CellName = unlist(CellName)) %>%
    dplyr::mutate(CellName = paste0("Sample_",i,"_",CellName)) %>%
    dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":"))
}

all.tcr.data <- Reduce(rbind, tcr)
all.tcr.data <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/01.geo/all.tcr.data.rds.gz")

#---------- Load data -----------#

dirs <- list.dirs(out.path)
dirs <- dirs[-1]

load10x <- function(path, gene.count = 500, umi.low = 500, umi.high = 25000){ 
  
  gene_path <- "/home/pauling/projects/01_data/03_gene/coding_gene.rds.gz"
  cd_gene <- readr::read_rds(gene_path)
  
  matr <- Read10X(path)
  if(length(matr) == 2){
    matr <- matr[[1]]
  }
  ig.name <- rownames(matr)
  ig.name <- ig.name[stringr::str_detect(ig.name,"^IGH*")]
  over_gene <- c(intersect(rownames(matr), cd_gene$gene_name))
  matr <- matr[over_gene,]
  total_umi <- colSums(as.matrix(matr))
  gene_count <- colSums(as.matrix(matr) > 0, na.rm = T)
  matr <- matr[, c(total_umi < umi.high & total_umi > umi.low & gene_count > gene.count)]
  
  return(matr)
}

matr <- list()
for (i in 1:length(dirs)) {
  matr[[i]] <- load10x(path = dirs[i])
  print("done")
}

meta <- list()
for (i in 1:length(matr)) {
  colnames(matr[[i]]) <- paste("Sample",i,colnames(matr[[i]]), sep = "_")
  meta[[i]] <- tibble(cellid = colnames(matr[[i]])) %>%
    dplyr::mutate(sample = paste0("Sample.",i))
}

all.matr <- Reduce(Matrix::cBind, matr)

meta <- Reduce(rbind, meta) %>% as.data.frame()
meta <- as.data.frame(meta)
rownames(meta) <- meta$cellid

sce <- CreateSeuratObject(counts = all.matr, meta.data = meta)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes, do.center = F, do.scale = F)
sce <- FindVariableFeatures(sce, nfeatures = 1500)
sce <- RunPCA(sce, features = VariableFeatures(sce))
sce <- RunUMAP(sce, dims = 1:15)
sce <- FindNeighbors(sce, reduction = "pca", dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5)

sce <- sce[,Idents(sce) != 6]

#----------- batch information

bda <- tibble(
  batch = c("P1.tumor","P1.tumor","P2.tumor","P2.LNmet","P2.LNmet","P3.tumor","P3.tumor","P3.LNmet",
            "P4.tumor","P4.LNmet","P4.LNmet","P5.tumor","P5.LNmet"),
  sample = paste0("Sample.",1:13)
)

bda <- tibble(sample = sce$sample) %>% dplyr::left_join(bda, by = "sample")
sce$batch <- bda$batch

#----------- batch correction

pca <- sce@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca, obs=sce$batch)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 3)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce <- loadumap(sce, umap, bb.res)
DimPlot(sce.pro, label = T)

#----------- annotation

sce <- sce[,Idents(sce) != "22"]
clusters <- Idents(sce)
cluster.name <- rep("cells",ncol(sce))

cluster.name[clusters %in% c(7,17,2)] <- "stem-like"
cluster.name[clusters %in% c(10,8)] <- "Transitory"
cluster.name[clusters %in% c(23)] <- "Prolif."
cluster.name[cluster.name == "cells"] <- "TD"

Idents(sce)

sce.pro <- sce
Idents(sce.pro) <- cluster.name

sce.pro@reductions$umap@cell.embeddings %>%
  as.tibble() %>%
  dplyr::mutate(label = Idents(sce.pro)) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = label), size=0.8, alpha = 1) +
  theme_classic() +
  scale_color_manual(values = c("#009ACD","#FF7F50","#CD2626", "#8968CD"))

#----------- Feature plot

FeaturePlot(sce.pro, features = "CXCL13", pt.size = 0.3, cols = c("grey92", "#EE0000", "#CD2626")) + NoAxes()



#----- TCR + scRNA-seq

sce.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/sce.rds.gz")
overlap.cells <- intersect(colnames(sce.pro), all.tcr.data$CellName)

sce.pro <- sce.pro[,colnames(sce.pro) %in% overlap.cells]
tda <- data.frame(
  row.names = all.tcr.data$CellName,
  cellname = all.tcr.data$CellName,
  clone.id = all.tcr.data$clone.id
)
tda <- tda[colnames(sce.pro),]
sce.pro$clone.id <- tda$clone.id

sce.pro %>% readr::write_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/02.HNSCC.nature/sce.tcr.rds.gz", compress = "gz")
