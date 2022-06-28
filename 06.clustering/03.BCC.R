
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

#---------- load data -----------#

bcc.sce <- readr::read_rds("/raid1/pauling/projects/01_data/12_NM_PD1/02.bcc.tcell.sce.rds.gz")
mda <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/02.meta.txt.gz")

all.tcr <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/GSE123813_bcc_tcr.txt.gz")
all.tcr <- all.tcr %>%
  dplyr::filter(stringr::str_detect(cdr3s_aa, "TRB")) %>%
  dplyr::mutate(
    sep = purrr::map(
      .x = cdr3s_aa,
      .f = function(.x){
        stringr::str_split(.x, pattern = ";", n = Inf, simplify = F)[[1]]
      }
    )
  ) %>%
  dplyr::mutate(
    TRB = purrr::map_chr(
      .x = sep,
      .f = function(.x){
        for (i in 1:length(.x)) {
          if(stringr::str_detect(.x[i], "TRB")){
            return(.x[i])
            break
          }
        }
      }
    )
  ) %>%
  dplyr::mutate(
    n = purrr::map_dbl(
      .x = sep,
      .f = function(.x){
        length(.x)
      }
    )
  )

tcell.filt <- bcc.sce[,colnames(bcc.sce) %in% all.tcr$cdr3s_nt]
all.tcr <- all.tcr %>% dplyr::filter(cdr3s_nt %in% colnames(tcell.filt))
all.tcr <- all.tcr %>% dplyr::rename(clone.id = TRB, cellid = cdr3s_nt)
all.tcr <- all.tcr %>% dplyr::inner_join(mda[,c(1,2,3)], by = c("cellid" = "cell.id"))

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/03.BCC.rds.gz")
pda <- pda %>% dplyr::filter(clone == "Tex") %>% dplyr::filter(size > 0) %>% dplyr::arrange(desc(size))

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% pda$clone.id)

sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cellid]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellid
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]

sce.pro$trb <- all.tcr.pro$clone.id

#---------- samples ---------#

rep.data <- tibble(patient = sce.pro$patient) %>%
  dplyr::mutate(response = ifelse(patient %in% c("su001","su002","su003","su004","su009","su012"), "PR", "PD"))

sce.pro$response <- rep.data$response

#----------- batch correction

sce.pro <- FindVariableFeatures(sce.pro, nfeatures = 1500)
sce.pro <- RunPCA(sce.pro, features = VariableFeatures(sce.pro))
sce.pro <- RunUMAP(sce.pro, dims = 1:15)

sce.sub <- sce.pro[,sce.pro$patient %in% unique(sce.pro$patient)[table(sce.pro$patient) > 19]]

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
sc$tl$leiden(adata, resolution = 1.5)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce.sub <- loadumap(sce.sub, umap, bb.res)
DimPlot(sce.sub, label = T)

FeaturePlot(sce.sub, features = c("GZMK","HAVCR2","LAG3","IL7R"))

#---------- Annotation

clusters <- Idents(sce.sub)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(5,6)] <- "0.GZMK+Texp"
cluster.name[clusters %in% c(0)] <- "0.IL7R+Texp"
cluster.name[clusters %in% c(8)] <- "Prolif."
cluster.name[cluster.name == "cells"] <- "TD"

sce.sub$clustesrs <- Idents(sce.sub)
Idents(sce.sub) <- cluster.name

sce.sub %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/05.bcc/01.tumor.reactive.cells.rds.gz", compress = "gz")


DimPlot(sce.sub)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro1 <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/05.bcc/01.tumor.reactive.cells.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  response = sce.pro$response,
  treatment = sce.pro$treatment
) %>%
  dplyr::arrange(desc(cluster)) %>%
  dplyr::filter(response == "PD" & treatment == "post") %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 3, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#3CB371",comb.d3[10],"grey88","grey88")) +
  xlim(1.5,10.2) +
  ylim(1.9,9.5)

#----------- prop analysis

sce.pro1 <- bcc.sce
tmp1 <- table(paste(sce.pro1$response, sce.pro1$treatment, sce.pro1$patient, sep = "."), Idents(sce.pro1))
tmp1 <- as.data.frame.array(tmp1)
tmp1 <- tmp1[rowSums(tmp1) > 19,]
as.data.frame.matrix(tmp1/rowSums(tmp1)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  ggplot(aes(factor(response, levels = c("PR.pre","PR.post","PD.pre","PD.post")), `0.GZMK+Texp`+`0.IL7R+Texp`)) +
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
    axis.title = element_text(size = 13, colour = "black")) +
  scale_color_manual(values = c("#009ACD", "#4F94CD", "#CD2626", "#EE799F")) +
  ylim(0,1)

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::select(response, `IL7R+Texp`,`GZMK+Texp`)

a$prop <- rowSums(a[,-1])
t.test(a$prop[1:6], a$prop[7:15])
