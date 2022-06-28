
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#--- load data

tcell.sce <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/02.tcell.sce.rds.gz")
all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")

all.tcr <- all.tcr %>%
  dplyr::mutate(clone.id = paste(`Identifier(Alpha1)`,`Identifier(Beta1)`, sep = ":")) %>%
  as.tibble()

#--- tumor-reactive CD4 T cell clones

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/01.NSCLC.1.rds.gz")
tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.4 & CD8A == 0, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% tex.clones)

sce.pro <- tcell.sce[,tcell.sce$cellid %in% all.tcr.pro$cellid]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellid
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]

sce.pro$trb <- all.tcr.pro$clone.id

sce.pro <- sce.pro[,sce.pro$patient %in% c("P1","P10","P13","P19","P30","P33","P35","P36","P37","P38")]
FeaturePlot(sce.pro, features = "PDCD1")

#--- batch correction

sce.pro <- FindVariableFeatures(sce.pro, nfeatures = 1500)
sce.pro <- RunPCA(sce.pro, features = VariableFeatures(sce.pro))
sce.pro <- RunUMAP(sce.pro, dims = 1:15)

sce.sub <- sce.pro[,sce.pro$patient %in% unique(sce.pro$patient)[table(sce.pro$patient) > 19]]

use_python("/home/pauling/anaconda3/bin/python")
pca <- sce.sub@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca[,1:(20-1)], obs=sce.sub$patient)
sc$tl$pca(adata)
adata$obsm$X_pca = pca[,1:(20-1)]
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 0.8)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce.sub <- loadumap(sce.sub, umap, bb.res)
DimPlot(sce.sub, label = T)

FeaturePlot(sce.sub[,sce.sub$patient == "P13" & sce.sub$num == 0], features = c("IL7R"))
markers <- FindAllMarkers(sce.sub, only.pos = T)

sce.sub <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/10.CD4.tumor.reactive.T.clustering/01.NSCLC1.rds.gz")

#--- gene heatmap

cd4.pro <- tcell.sce[,tcell.sce$cellid %in% cd4.sce$cellid]
mda <- tibble(
  label = as.character(Idents(cd4.sce)),
  cellid = cd4.sce$cellid
)

mda <- mda %>%
  dplyr::filter(label != "CD4_C9-Prolif.") %>%
  dplyr::mutate(label = ifelse(label %in% c("CD4_C1-Naive","CD4_C7-Th1-like","CD4_C8-Treg"), label, "Other")) %>%
  dplyr::filter(cellid %in% cd4.pro$cellid)

genes <- c("CXCL13","PDCD1","BHLHE40","TOX","GNG4","IFNG","IL21","CD82","TPI1","HAVCR2","TIGIT","ENTPD1","TNFRSF4","BATF")
clusters <- unique(mda$label)
mean.matr <- Matrix::Matrix(data = 0, nrow = length(clusters), ncol = length(genes))

for (i in 1:length(clusters)) {
  use.cells <- mda %>%
    dplyr::filter(label == clusters[i]) %>%
    dplyr::pull(cellid)
  
  mean.matr[i,] <- Matrix::rowMeans(cd4.pro@assays$RNA@scale.data[genes,use.cells])
}

rownames(mean.matr) <- clusters
colnames(mean.matr) <- genes

as.matrix(mean.matr) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "label") %>%
  dplyr::mutate_if(is.double, funs((. - mean(.))/sd(.))) %>%
  tidyr::gather(key = "gene", value = "expr", -label) %>%
  ggplot(aes(factor(label, levels = c("CD4_C7-Th1-like","CD4_C8-Treg","Other","CD4_C1-Naive")), factor(gene, levels = rev(genes)), fill = expr)) +
  geom_tile() +
  theme(axis.title = element_text(size = 0)) +
  theme(axis.text = element_text(size = 10)) +
  theme(legend.title = element_text(size = 0)) +
  theme(legend.text = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  theme(axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color="black",angle = 45, hjust = 1)) +
  scale_fill_distiller(palette = "RdBu")
