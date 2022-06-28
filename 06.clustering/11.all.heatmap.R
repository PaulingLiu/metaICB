
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)

#--- load data

lung1.sce <- readr::read_rds("projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/07.Tex.rds.gz")
lung2.sce <- readr::read_rds("projects/07_pan.cancer.treatment/02.data/01.lung.nature/04.sample.batch.cxcl13.sce.rds.gz")
bcc.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/05.bcc/01.tumor.reactive.cells.rds.gz")
scc.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/06.scc/01.tumor.reactive.cells.rds.gz")
breast1.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/01.tumor.reactive.cells.batch1.rds.gz")
breast2.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/02.tumor.reactive.cells.batch2.rds.gz")
breast3.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/02.tumor.reactive.cells.batch1.rds.gz")
renal1.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/02.RCC.data.1/01.tumor.reactive.cells.batch1.rds.gz")
renal2.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/03.RCC.data.2/01.tumor.reactive.cells.batch1.rds.gz")

#--- lung clustering

lung1.sce <- lung1.sce[,lung1.sce$patient %in% c("P1","P10","P13","P19","P30","P33","P35")]

#--- renal data genes

cdg <- readr::read_rds("/raid1/pauling/projects/01_data/03_gene/coding_gene.rds.gz")
cdg <- cdg %>% dplyr::select(gene_name, gene_id) %>% dplyr::distinct(gene_name, .keep_all = T)
cda <- as.data.frame(cdg)
rownames(cdg) <- cdg$gene_name

getgene <- function(.x){
  cdg[.x,] %>% dplyr::pull(gene_id)
}

#--- heatmap matrix

get.expr <- function(object, genes, data){
  object <- ScaleData(object, do.center = T, do.scale = T, features = genes, scale.max = 0.92)
  tmp.expr <- as.tibble(t(object@assays$RNA@scale.data[genes,]))
  tmp.expr$cluster = as.character(Idents(object))
  tmp.expr <- tmp.expr %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_if(is.double, ~ mean(.x, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster = paste(cluster, data, sep = "."))
}

genes <- c("UBE2C","RRM2","TYMS","CDK1","STMN1","MKI67","TUBB","DUSP4","PRF1","ITGAE","CTSW","GZMB","HAVCR2","ENTPD1",
           "LAYN","TIGIT","LAG3","CTLA4","PDCD1","TOX",
           "LMNA","ANXA1","CST7","ENC1",
           "CD44","DUSP2","GPR183","KLRG1","GZMK","CD28","TCF7","CCR7","SELL","IL7R")

lung2.clusters <- as.character(Idents(lung2.sce))
lung2.clusters[lung2.clusters == "HSP+ TD"] <- "TD"
lung2.sce.pro <- lung2.sce
Idents(lung2.sce.pro) <- lung2.clusters

expr.lung1 <- get.expr(lung1.sce, genes, "lung1")
expr.lung2 <- get.expr(lung2.sce.pro, genes, "lung2")
expr.bcc <- get.expr(bcc.sce, genes, "bcc")
expr.scc <- get.expr(scc.sce, genes, "scc")
expr.breast1 <- get.expr(breast1.sce, genes, "breast1")
expr.breast2 <- get.expr(breast2.sce, genes, "breast2")
expr.breast3 <- get.expr(breast3.sce, genes, "breast3")
expr.rc1 <- get.expr(renal1.sce, getgene(genes), "rc1")
expr.rc2 <- get.expr(renal2.sce, genes, "rc2")

colnames(expr.rc1) <- c("cluster", genes)

expr <- rbind(
  expr.lung1, expr.lung2, expr.bcc, expr.scc,
  expr.breast1,expr.breast2,expr.breast3,
  expr.rc1, expr.rc2
)

expr <- expr %>% dplyr::mutate_if(is.double,funs((. - mean(.))/sd(.)))

#--- heatmap plot

clusters <- expr$cluster
clusters <- c(clusters[1:5],NA,clusters[6:31])
cluster.order <- c(clusters[c(c(1,7,11,15,25),c(2,5,10,14,18,21,24),c(4,9,13,17,20,23,26+1,28+1,30+1),c(3,8,12,16,19,22,26,27+1,29+1))])

expr %>%
  tidyr::gather(key = "gene", value = "expr", -cluster) %>%
  dplyr::mutate(expr = ifelse(expr < -2, -2, expr)) %>%
  dplyr::mutate(expr = ifelse(expr > 2, 2, expr)) %>%
  ggplot(aes(factor(cluster, levels = cluster.order), factor(gene, levels = genes))) +
  geom_tile(aes(fill = expr)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black", angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = "", y = "")


##################################################################################################

#-- Prop of Texp cells

p1 <- a %>%
  ggplot(aes(factor(group, levels = c("pre","responsive","non.res")), prop)) +
  geom_boxplot(aes(color = group),outlier.size = -1) +
  geom_jitter(aes(color = group), width = 0.1, size = 2, height = 0) +
  theme_classic() +
  labs(x = "", y = "Proportion of Texp cells in all Tex cells") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.y = element_text(size = 0, colour = "black")) +
  scale_color_manual(values = c("#009ACD", "#EE799F","#CD2626","#4F94CD")) +
  ylim(0,1)

ggarrange(p1,p2,p3,p5,p6,p7, nrow = 1, widths = c(1.5,1.1,1.85,1.85,1.1,1.85))


##################################################################################################

#-- barplot

prop1 <- a %>% 
  dplyr::filter(group == "responsive") %>% 
  dplyr::rename(response = group) %>% 
  dplyr::mutate(group = "01.lung") %>%
  dplyr::select(-patient)

prop2 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","tumor"), sep = "\\.") %>%
  dplyr::mutate(prop = `GZMK+Texp` + `IL7R+Texp`) %>%
  dplyr::filter(response == "MPR") %>%
  dplyr::select(response, prop) %>%
  dplyr::mutate(group = "02.lung")

prop3 <- as.data.frame.matrix(tmp1/rowSums(tmp1)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(prop = `0.GZMK+Texp`+`0.IL7R+Texp`) %>%
  dplyr::select(response, prop) %>%
  dplyr::filter(response == "PR.post") %>%
  dplyr::mutate(group = "03.BCC")

prop4 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(prop = `0.GZMK+Texp`+`0.IL7R+Texp`) %>%
  dplyr::select(response, prop) %>%
  dplyr::filter(response == "PR.post") %>%
  dplyr::mutate(group = "04.SCC")

prop5 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(prop = `0.GZMK+Texp`) %>%
  dplyr::select(response, prop) %>%
  dplyr::filter(response == "E.On") %>%
  dplyr::mutate(group = "05.breast")

prop6 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(prop = `0.GZMK+Texp`) %>%
  dplyr::select(response, prop) %>%
  dplyr::filter(response == "E.On") %>%
  dplyr::mutate(group = "06.breast")

prop7 <- as.data.frame.matrix(tmp/rowSums(tmp)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(sample = stringr::str_remove(sample,"-treatment")) %>%
  tidyr::separate(sample, c("response","treatment","patient"), sep = "\\.") %>%
  dplyr::mutate(response = paste(response, treatment, sep = ".")) %>%
  dplyr::mutate(prop = `0.GZMK+Texp`+`0.IL7R+Texp`) %>%
  dplyr::select(response, prop) %>%
  dplyr::filter(response == "PR.Post") %>%
  dplyr::mutate(group = "07.breast")

prop8 <- tibble(
  prop = c(0,0,0),
  response = c("R","R","R"),
  group = "08.renal"
)

prop9 <- tibble(
  prop = c(0,0),
  response = c("R","R"),
  group = "09.renal"
)

pda <- rbind(prop1, prop2, prop3, prop4, prop5, prop6, prop7, prop8, prop9)
pda <- readr::read_rds("projects/07_pan.cancer.treatment/02.data/04.all.Texp.post.prop.rds.gz")
rank.level <- c('07.breast',"01.lung","02.lung","06.breast","04.SCC","03.BCC","05.breast","08.renal","09.renal")
pda %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(mean = mean(prop), n = length(prop), sd = sd(prop)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(group, levels = rank.level), mean)) +
  geom_col(aes(fill = group, color = group), alpha = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_jitter(aes(group, prop, colour = group), size = 2, data = pda, width = 0.1) +
  geom_errorbar(aes(ymin = mean, ymax = ymax), width = 0.4) +
  theme_classic() +
  scale_fill_manual(values = comb.d3[c(1:7,10,12)]) +
  scale_color_manual(values = comb.d3[c(1:7,10,12)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 35, hjust = 1)
  )


prop.aov <- aov(prop~group, data=pda)
summary(prop.aov)



#--- TOX

lung2.clusters <- as.character(Idents(lung2.sce))
lung2.clusters[lung2.clusters %in% c("HSP+ TD", "Prolif.")] <- "TD"

a <- tibble(
  expr = lung2.sce@assays$RNA@scale.data["TOX",],
  label = lung2.clusters,
  patient = lung2.sce$patient
) %>%
  dplyr::filter(patient != "MD01-005") %>%
  dplyr::filter(expr > 0.4) %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(n = n())

b <- tibble(
  expr = lung2.sce@assays$RNA@scale.data["TOX",],
  label = lung2.clusters,
  patient = lung2.sce$patient
) %>%
  dplyr::filter(patient != "MD01-005") %>%
  #dplyr::sample_n(700, replace = T) %>%
  dplyr::ungroup() %>%
  #dplyr::filter(expr > 0.4) %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(n = n(), mea = mean(expr))

a$n/b$n

tibble(
  Frac = c(0.28,0.39,0.54),
  Mean = c(0.5,0.62,0.84),
  group = c("TOXlow","TOXint","Terminal")
) %>%
  ggplot(aes(factor(group, levels = c("TOXlow","TOXint","Terminal")), Mean)) +
  geom_point(aes(fill = Frac), size = 5, shape = 21, color = "black") +
  theme_classic() +
  labs(
    x = "",
    y = "log2(CPM + 1)"
  ) +
  theme(
    axis.title = element_text(size = 13, color = "black"), 
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(color = "black"), 
    axis.text.x = element_text(color = "black")) +
  scale_fill_distiller(palette = "GnBu", direction = 1)

tibble(
  Frac = c(0.2288752,0.2197,0.2778869),
  Mean = c(0.335,0.344,0.391),
  group = c("TOXint","TOXlow","Terminal")
) %>%
  ggplot(aes(factor(group, levels = c("TOXlow","TOXint","Terminal")), Frac)) +
  geom_point(aes(fill = Mean), size = 5, shape = 21, color = "black") +
  theme_classic() +
  labs(
    x = "",
    y = "Fraction of TOX"
  ) +
  theme(
    axis.title = element_text(size = 13, color = "black"), 
    axis.text = element_text(size = 12, color = "black"),
    axis.text.y = element_text(color = "black"), 
    axis.text.x = element_text(color = "black")) +
  scale_fill_distiller(palette = "GnBu", direction = 1) +
  ylim(0.215,0.28)

#--- TOX

FeaturePlot(lung1.sce, features = "HAVCR2", cols = c("grey92", "#EE0000", "#CD2626"), reduction = "umap", pt.size = 0.5) + NoAxes()

library(Nebulosa)
plot_density(breast3.sce, c("ZNF683"), joint = TRUE, combine = FALSE, adjust = 0.9, method = "wkde")
