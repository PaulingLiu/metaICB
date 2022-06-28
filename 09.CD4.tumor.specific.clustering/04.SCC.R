
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

bcc.sce <- readr::read_rds("/raid1/pauling/projects/01_data/12_NM_PD1/04.scc.tcell.rds.gz")
bcc.sce <- NormalizeData(bcc.sce)
bcc.sce <- ScaleData(bcc.sce, features = rownames(bcc.sce), do.center = F, do.scale = F)
mda <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/GSE123813_scc_metadata.txt.gz")

all.tcr <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/bcc_scc_tcr_vdj.txt.gz")
all.tcr <- all.tcr %>%
  dplyr::filter(stringr::str_detect(cdr3s_nt, "TRB")) %>%
  dplyr::mutate(
    sep = purrr::map(
      .x = cdr3s_nt,
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

all.tcr <- all.tcr %>% 
  dplyr::select(clonotype.id, TRB, n) %>%
  dplyr::rename(clone.id = TRB, cellid = clonotype.id)

tcell.filt <- bcc.sce[,colnames(bcc.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))
all.tcr <- all.tcr %>% dplyr::inner_join(mda[,c(1,2,3)], by = c("cellid" = "cell.id"))

#---------- Tex clones ----------#

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/04.SCC.rds.gz")
tex.clones <- pda %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5 & CD8A == 0, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% tex.clones)

sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cellid]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellid
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]

sce.pro$trb <- all.tcr.pro$clone.id

#---------- samples ---------#

rep.data <- tibble(patient = sce.pro$patient) %>%
  dplyr::mutate(response = ifelse(patient %in% c("su010","su011"), "PR", "PD"))

sce.pro$response <- rep.data$response

#----------- batch correction

sce.pro <- FindVariableFeatures(sce.pro, nfeatures = 2000)
sce.pro <- RunPCA(sce.pro, features = VariableFeatures(sce.pro))
sce.pro <- RunUMAP(sce.pro, dims = 1:15)

pca <- sce.pro@reductions$pca@cell.embeddings
anndata = import("anndata",convert=FALSE)
sc = import("scanpy",convert=FALSE)
np = import("numpy",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
adata = anndata$AnnData(X=pca, obs=sce.pro$patient)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata, batch_key=0)
sc$tl$umap(adata)
sc$tl$leiden(adata, resolution = 1)
umap = py_to_r(adata$obsm['X_umap'])

tibble::tibble(
  leiden.res.1.bbknn = py_to_r(np$asarray(adata$obs$leiden)), 
  UMAP1 = umap[,1], 
  UMAP2 = umap[,2]) -> bb.res

sce.pro <- loadumap(sce.pro, umap, bb.res)
DimPlot(sce.pro, label = T)
FeaturePlot(sce.pro, features = c("CXCL13","FOXP3","HAVCR2","IL7R","MKI67","TCF7"), ncol = 3)
FeaturePlot(sce.pro, features = c("TCF7"))

sce.sub <- sce.pro[,Idents(sce.pro) != 6]
sce.sub %>% readr::write_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/10.CD4.tumor.reactive.T.clustering/04.SCC.rds.gz", compress = "gz")
#---------- Annotation

clusters <- Idents(sce.pro)
cluster.name <- rep("cells", length(clusters))
cluster.name[clusters %in% c(0)] <- "0.GZMK+Texp"
cluster.name[clusters %in% c(3,6)] <- "0.IL7R+Texp"
cluster.name[clusters %in% c(9)] <- "Prolif."
cluster.name[cluster.name == "cells"] <- "TD"

sce.pro$clustesrs <- Idents(sce.pro)
Idents(sce.pro) <- cluster.name

sce.pro %>% readr::write_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/06.scc/01.tumor.reactive.cells.rds.gz", compress = "gz")

DimPlot(sce.pro)
FeaturePlot(sce.pro[,sce.pro$treatment == "pre" & sce.pro$response == "PR"], features = "HAVCR2", pt.size = 3, cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#---------- UMAP plot

sce.pro <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/06.scc/01.tumor.reactive.cells.rds.gz")

tibble(
  umap1 = sce.pro@reductions$umap@cell.embeddings[,1],
  umap2 = sce.pro@reductions$umap@cell.embeddings[,2],
  cluster = as.character(Idents(sce.pro)),
  response = sce.pro$response,
  treatment = sce.pro$treatment
) %>%
  dplyr::arrange(desc(cluster)) %>%
  dplyr::filter(response == "PD" & treatment == "pre") %>%
  ggplot(aes(umap1, umap2)) +
  geom_point(aes(color = cluster), size = 3, alpha = 0.8) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#3CB371",comb.d3[10],"grey88","grey88")) #+
xlim(-1.8,8.5) +
  ylim(-2,9)

#----------- prop analysis
sce.pro <- scc.sce
tmp <- table(paste(sce.pro$response, sce.pro$treatment, sce.pro$patient, sep = "."), Idents(sce.pro))
tmp <- as.data.frame.array(tmp)
tmp <- tmp[rowSums(tmp) > 19,]

tmp <- rbind(tmp, tmp1)

p3 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(tumor = c(rep("SCC",7),rep("BCC",10))) %>%
  ggplot(aes(factor(response, levels = c("PR.pre","PR.post","PD.pre","PD.post")), `0.IL7R+Texp`+`0.GZMK+Texp`)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_jitter(aes(color = response, shape = tumor), width = 0.1, size = 1.8, stroke = 1) +
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
  scale_color_manual(values = c("#009ACD", "#4F94CD", "#CD2626", "#EE799F")) +
  ylim(0,1) +
  scale_shape_manual(values = c(2,1))

#----- statistical test

a <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::select(response, `IL7R+Texp`,`GZMK+Texp`)

a$prop <- rowSums(a[,-1])
t.test(a$prop[1:6], a$prop[7:15])
