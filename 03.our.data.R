
#----- library -----#

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#----- load data ------#

cd8.sce <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/06.CD8.rds.gz")
all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")

DimPlot(cd8.sce)
table(cd8.sce$patient)

#----- responsive tumors -----#

cd8.rp <- cd8.sce[,cd8.sce$patient %in% c("P1","P10","P13","P19","P30","P33","P35")]
cd8.rp <- cd8.rp[,cd8.rp$num != 2]
DimPlot(cd8.rp, label = T)

levels(cd8.rp)
new.cluster <- c("pre","pre","pre","Tex","pre","pre","Prolif.")
names(new.cluster) <- levels(cd8.rp)

cd8.rp <- RenameIdents(cd8.rp, new.cluster)
DimPlot(cd8.rp, cols = c("#E0E0E0", "#EE6363", "#009ACD")) + NoLegend() + NoAxes()

#----- TCR clones -----#

tex.sce <- cd8.rp[,Idents(cd8.rp) %in% c("Tex")]

clone.cells <- function(.x, .y){
  .y %>%
    dplyr::filter(cellid %in% colnames(.x)) %>%
    dplyr::pull(clone.id) -> tmp.clones
  
  .y %>%
    dplyr::filter(clone.id %in% tmp.clones) %>%
    dplyr::pull(cellid) -> tmp.cells
  
  return(tmp.cells)
}

cd8.sub <- cd8.rp[,intersect(colnames(cd8.rp), all.tcr$cellid)]
tex.cells <- clone.cells(tex.sce, all.tcr)

Idents(cd8.sub) <- cd8.sub$major.cluster
DimPlot(cd8.sub, cols = c("#6CA6CD", "#CD5C5C", "#008B8B"))

tex.sce <- cd8.sub[,colnames(cd8.sub) %in% tex.cells]
cd8.sub <- cd8.clone
clusters <- ifelse(cd8.sub$clone.id %in% pre.data$clone, "Tex", "Not qualified")
clusters <- ifelse(cd8.sub$clone.id %in% pre.data$clone[pre.data$group == "nonTex"], "nonTex", clusters)
Idents(cd8.sub) <- clusters

cd8.sub@reductions$umap@cell.embeddings %>%
  as.tibble() %>%
  dplyr::mutate(label = Idents(cd8.sub), num = cd8.sub$num) %>%
  dplyr::arrange(label) %>%
  dplyr::filter(num != 0) %>%
  ggplot(aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = label), size = 0.1) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#D9D9D9", "#9F79EE", "#3CB371"))
c("#3CB371", "#7B68EE")
c("#D02090", "#EE5C42", "#9F79EE", "#7B68EE")
#----- pre-treatment tumor clones -----#

all.tcr <- all.tcr[all.tcr$cellid %in% colnames(cd8.rp),]
cd8.clone <- cd8.rp[,all.tcr$cellid]
cd8.clone$clone.id <- all.tcr$clone.id
cd8.clone$clone.pro <- paste(cd8.clone$patient, cd8.clone$num, cd8.clone$clone.id, sep = ".")
cd8.clone.pre <- cd8.clone[,cd8.clone$num == 0]

use.clones <- tibble(clones = cd8.clone.pre$clone.id) %>%
  dplyr::count(clones) %>%
  dplyr::filter(n >= 10)

clone.data <- tibble(
  clone = cd8.clone.pre$clone.id,
  cluster = Idents(cd8.clone.pre)
) %>%
  dplyr::filter(clone %in% use.clones$clones) %>%
  dplyr::count(clone, cluster) %>%
  tidyr::spread(key = "cluster", value = "n") %>%
  dplyr::mutate_if(is.integer, funs(ifelse(is.na(.), 0, .)))

clone.data[,2:4] <- clone.data[,2:4]/rowSums(clone.data[,2:4])
clone.data <- clone.data %>% 
  dplyr::mutate(group = ifelse(Tex > 0.4, "Tex", "other")) %>%
  dplyr::mutate(group = ifelse(Tex == 0, "nonTex", group)) %>%
  dplyr::select(clone, group)

#-----  Classification -----#

cd8.clone.pre.scale <- ScaleData(cd8.clone.pre, do.scale = F, do.center = F)
cd8.clone.pre.scale <- cd8.clone.pre.scale[,cd8.clone.pre.scale$clone.id %in% use.clones$clones]

pre.data <- tibble(
  cxcl13 = cd8.clone.pre.scale@assays$RNA@scale.data["CXCL13",],
  clone = cd8.clone.pre.scale$clone.id,
  #type = cd8.clone.scale$type,
  tex = cd8.clone.pre.scale@assays$RNA@scale.data["ENTPD1",] + cd8.clone.pre.scale@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.pre.scale@assays$RNA@scale.data["ITGAE",] + cd8.clone.pre.scale@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::inner_join(clone.data, by = "clone") %>%
  dplyr::filter(group != "other") %>%
  dplyr::filter(group == "Tex" | (group == "nonTex" & m < 0.2))

pre.data %>%
  #dplyr::filter(clone %in% c(clone1$trb, clone2$trb)) %>%
  ggplot(aes(m, m2)) +
  geom_point(aes(color = group), size = 2.5) +
  geom_vline(xintercept = 0.2) +
  geom_hline(yintercept = 2) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "CXCL13 expression",
    y = "Exhaustion score"
  ) +
  scale_color_manual(values = c("#00868B", "#CD2626")) +
  ylim(0,5.8)

pre.data %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(n = sum(n))

#----- Differential expression of post-treatment Tex clones vs. non-Tex clones-----#

cd8.clone.post <- cd8.clone[,cd8.clone$num != 0]
cd8.clone.post <- cd8.clone.post[,cd8.clone.post$clone.id %in% pre.data$clone]

cd8.clone.post$type = ifelse(cd8.clone.post$clone.id %in% clone.data$clone[clone.data$group == "Tex"], "Tex", "nonTex")
Idents(cd8.clone.post) <- cd8.clone.post$type
#markesr <- FindAllMarkers(cd8.clone.post, only.pos = T)

tibble(
  expr1 = Matrix::rowMeans(cd8.clone.post@assays$RNA@scale.data[,cd8.clone.post$type == "Tex"]),
  expr2 = Matrix::rowMeans(cd8.clone.post@assays$RNA@scale.data[,cd8.clone.post$type != "Tex"]),
  gene = rownames(cd8.clone.post)
) -> pda

fit <- loess(expr1 ~ expr2, data = pda, span = 0.8)
prd <- predict(fit, pda$expr2)

pda %>%
  dplyr::mutate(prd = prd) %>%
  dplyr::mutate(name = ifelse(gene %in% c("CTLA4","PDCD1","CXCL13","ENTPD1","LAG3","TIGIT","BTAF","HAVCR2","TNFRSF9","GZMB","CCL5"), gene, NA)) %>%
  dplyr::mutate(sig = ifelse(gene %in% markesr$gene[markesr$p_val_adj < 1e-3 & markesr$avg_logFC > 0.5],"sig","not")) %>%
  dplyr::arrange(sig) %>%
  ggplot(aes(expr2, expr1)) +
  geom_point(aes(color = sig)) +
  geom_line(aes(expr2, prd), lwd = 0.9, color = "black", linetype = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("#DEDEDE", "#CD2626")) +
  labs(
    x = "Average expresion in Tex clones",
    y = "Average expresion in non-Tex clones"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  geom_text_repel(aes(label = name),
                  size = 4.3,
                  box.padding = 0.7,
                  segment.color = "black",
                  show.legend = FALSE)


#-----  Classification -----#

tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  type = cd8.clone.post$type,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post@assays$RNA@scale.data["ITGAE",] + cd8.clone.post@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  ggplot(aes(m, m2)) +
  geom_point(aes(color = type), size= 3) +
  geom_vline(xintercept = 0.15) +
  geom_hline(yintercept = 1) +
  #scale_y_log10(lim = c(1,3), breaks = c(1,3)) +
  #scale_x_log10(lim = c(0.9,5), breaks = c(1,5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "CXCL13 expression",
    y = "Exhaustion score"
  ) +
  scale_color_manual(values = c("#00868B", "#CD2626"))

tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  type = cd8.clone.post$type,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post@assays$RNA@scale.data["ITGAE",] + cd8.clone.post@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  ggplot(aes(1+m, 1+m2)) +
  geom_point(aes(shape = type), size= 2.5, stroke = 1) +
  geom_vline(xintercept = 1.16) +
  geom_hline(yintercept = 3.15) +
  scale_y_log10(lim = c(1,7), breaks = c(1,5)) +
  scale_x_log10(lim = c(0.9,5), breaks = c(1,5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "CXCL13 expression",
    y = "Exhaustion score"
  ) +
  scale_shape_manual(values = c(3,1))

tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  type = cd8.clone.post$type,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post@assays$RNA@scale.data["ITGAE",] + cd8.clone.post@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(n = sum(n))

cal.accuracy <- function(.x){
  tmp.1 <- tibble(
    cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
    trb = cd8.clone.post$clone.id,
    type = cd8.clone.post$type,
    tex = cd8.clone.post@assays$RNA@scale.data[.x,]
  ) %>%
    dplyr::group_by(trb, type) %>%
    dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
    dplyr::filter(n > 5)
  
  cutoff <- sort(tmp.1$m2[tmp.1$type == "nonTex"], decreasing = T)
  tpr <- c()
  for (i in 1:6) {
    n1 <- tmp.1 %>%
      dplyr::filter(type == "Tex") %>%
      dplyr::filter(m2 > cutoff[i]) %>%
      dplyr::ungroup() %>%
      nrow(.)
    
    tpr[i] <- n1/nrow(tmp.1[tmp.1$type == "Tex",])
  }
  tibble(tpr = tpr,
         fpr = 0:5/51) %>%
    dplyr::mutate(gene = .x)
}

tibble(
  gene = c("CXCL13","CTLA4","PDCD1","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5")
) %>%
  dplyr::mutate(res = purrr::map(gene, cal.accuracy)) %>%
  dplyr::pull(res) %>%
  reduce(rbind) -> pda

pda$gene <- factor(pda$gene, levels = c("CXCL13","CTLA4","PDCD1","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5"))

pda %>%
  ggplot(aes(fpr, tpr)) +
  geom_line(aes(color = gene), lwd = 0.6) +
  geom_point(aes(color = gene)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlab("FPR") +
  ylab("TPR") +
  scale_color_manual(values = comb.d3[1:12]) +
  xlim(0,0.1)

#----- AUC score for each gene -----#

rda <- tibble(
  type = cd8.clone.post$type,
  trb = cd8.clone.post$clone.id,
  CXCL13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],CTLA4 = cd8.clone.post@assays$RNA@scale.data["CTLA4",],
  PDCD1 = cd8.clone.post@assays$RNA@scale.data["PDCD1",],ENTPD1 = cd8.clone.post@assays$RNA@scale.data["ENTPD1",],
  LAG3 = cd8.clone.post@assays$RNA@scale.data["LAG3",],LAYN = cd8.clone.post@assays$RNA@scale.data["LAYN",],
  TIGIT = cd8.clone.post@assays$RNA@scale.data["TIGIT",],BAG3 = cd8.clone.post@assays$RNA@scale.data["BAG3",],
  BATF = cd8.clone.post@assays$RNA@scale.data["BATF",],HAVCR2 = cd8.clone.post@assays$RNA@scale.data["HAVCR2",],
  TNFRSF9 = cd8.clone.post@assays$RNA@scale.data["TNFRSF9",],GZMB = cd8.clone.post@assays$RNA@scale.data["GZMB",],
  CCL5 = cd8.clone.post@assays$RNA@scale.data["CCL5",], 
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["CTLA4",] + 
    cd8.clone.post@assays$RNA@scale.data["LAYN",] + cd8.clone.post@assays$RNA@scale.data["HAVCR2",]
)

rda <- rda %>%
  dplyr::group_by(type, trb) %>%
  dplyr::summarise(
    CXCL13 = mean(CXCL13),CTLA4 = mean(CTLA4),PDCD1 = mean(PDCD1),
    ENTPD1 = mean(ENTPD1),LAG3 = mean(LAG3),LAYN = mean(LAYN),
    TIGIT = mean(TIGIT),BAG3 = mean(BAG3),BATF = mean(BATF),
    HAVCR2 = mean(HAVCR2),TNFRSF9 = mean(TNFRSF9),GZMB = mean(GZMB),
    CCL5 = mean(CCL5),n=n(), tex = mean(tex),
  ) %>%
  dplyr::filter(n > 5) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = ifelse(type == "Tex", 1, 0))

rda$response <- as.factor(rda$response)


roc.list <- roc(response ~ CXCL13+CTLA4 + PDCD1 + ENTPD1 + 
                  LAG3 + LAYN + TIGIT + BATF + HAVCR2 + 
                  TNFRSF9 + GZMB + CCL5, data = rda, legacy.axes = TRUE)

pROC::ggroc(roc.list, lwd = 0.65) +
  theme_classic() +
  scale_color_manual(values = comb.d3[1:12]) +
  theme(
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  xlab("FPR") 

#----- Feature plot -----#

FeaturePlot(cd8.clone.post[,cd8.clone.post$type == "Tex"], features = "CXCL13", cols = c("grey92", "#EE0000", "#CD2626")) + NoLegend() + NoAxes()

#----- compare with Nature paper -----#

oda <- tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  type = cd8.clone.post$type,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["CTLA4",] + 
    cd8.clone.post@assays$RNA@scale.data["LAYN",] + cd8.clone.post@assays$RNA@scale.data["HAVCR2",]
) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5)


eda %>%
  ggplot(aes(factor(type, levels = c("Tex","MANA","nonTex","Viral")), m)) +
  geom_hline(yintercept = 0.155) +
  geom_jitter(aes(color = type, fill = type), width = 0.15, height = 0, size = 2, stroke = 1.2, shape = 1) +
  #scale_y_sqrt(breaks = c(0,1,2)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#CD2626","#6CA6CD","#EE6363", "#00868B")) +
  scale_fill_manual(values = c("#CD2626","#6CA6CD","#EE6363", "#00868B")) +
  labs(
    x = "",
    y = "CXCL13 expression"
  )


#--- Differential expression analysis

clone.rank <- tibble(clone = cd8.clone[,cd8.clone$patient == "P19"]$clone.id) %>%
  dplyr::count(clone) %>%
  dplyr::arrange(desc(n))

#clone1.sce <- cd8.clone[,cd8.clone$clone.id %in% use.clones$clones]
#clone1.sce <- clone1.sce[,clone1.sce$patient == "P35"]
clone1.sce <- cd8.clone[,cd8.clone$clone.id == clone.rank$clone[1]]
clone1.sce <- SCTransform(clone1.sce, return.only.var.genes = F)


Idents(clone1.sce) <- clone1.sce$treatment

markers <- FindAllMarkers(clone1.sce, only.pos = T, logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)

exhaustion.gene <- c("DUSP4","RGCC","NR4A2","ITGAE","RGS1","CTLA4","PDCD1","ENTPD1","LAG3","TIGIT","BTAF","HAVCR2","GZMB","LAYN")
mem.gene <- c("GZMK","IL7R","IL6ST","GZMA","CCND3","FKBP5","TXNIP")

markers %>%
  dplyr::rename(FC = avg_logFC) %>%
  dplyr::mutate(name = ifelse(gene %in% c(exhaustion.gene, mem.gene), gene, NA)) %>%
  dplyr::mutate(FC = ifelse(cluster == "ut", -FC,FC)) %>%
  dplyr::mutate(sig = ifelse(p_val_adj < 0.01 & abs(FC) > 0.15, "yes", "no")) %>%
  dplyr::filter(FC < 2) %>%
  ggplot(aes(FC, -log10(p_val_adj))) +
  geom_point(aes(color = sig), size = 0.9) +
  geom_text_repel(aes(label = name),
                  parse = T,
                  size = 3.5,
                  box.padding = 0.52,
                  segment.color = "black",
                  show.legend = FALSE) +
  theme_classic() +
  xlim(-1.6,1.6) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#D6D6D6","#CD2626"))
c("#EE2C2C", "#CD2626")
#--- individual clones

use.clonees <- table(p19.sce$clone.id, p19.sce$treatment) %>%
  as.data.frame.array() %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::arrange(desc(ut)) %>%
  dplyr::filter(tr > 20 & ut > 20)

pda <- list()
for (i in 1:nrow(use.clonees)) {
  clone.sce <- p19.sce[,p19.sce$clone.id %in% use.clonees$cloneid[i]]
  pda[[i]] <- t(clone.sce@assays$RNA@scale.data[exhaustion.gene,]) %>%
    as.tibble() %>%
    dplyr::mutate(group = clone.sce$treatment) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_if(is.double, mean)
}

Reduce(rbind, pda) %>%
  dplyr::arrange(desc(group)) %>%
  dplyr::mutate(clone = c(1:nrow(use.clonees),1:nrow(use.clonees))) %>%
  dplyr::filter(clone %in% c(1,2,3,5:7,12:14,16:19,22,24:26,29)) %>%
  dplyr::mutate_if(is.double, funs((. - mean(.))/sd(.))) %>%
  tidyr::gather(key = gene, value = expr, -group, -clone) %>%
  dplyr::mutate(group = paste0(group, ".",clone)) -> pda.pro

pda.pro %>%
  ggplot(aes(factor(group, levels = unique(pda.pro$group)), gene)) +
  geom_tile(aes(fill = expr), color = 'white') +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic()

#--- individual patients

patients <- unique(cd8.clone$patient)
patients <- patients[-4]

sct.list <- list()
for (i in 1:length(patients)) {
  sct.list[[i]] <- cd8.clone[,cd8.clone$patient == patients[i]]
  sct.list[[i]] <- SCTransform(sct.list[[i]], return.only.var.genes = F)
  Idents(sct.list[[i]]) <- sct.list[[i]]$treatment
}


p19.sce <- cd8.clone[,cd8.clone$patient == "P19"]
p19.sce <- SCTransform(p19.sce, return.only.var.genes = F)
Idents(p19.sce) <- p19.sce$treatment

VlnPlot(sct.list[[1]], features = "HAVCR2")
exhaustion.gene <- c("DUSP4","NR4A2","ITGAE","RGS1","CTLA4","ENTPD1","LAG3","TIGIT","HAVCR2","GZMB","LAYN")
t(tmp.sce@assays$RNA@scale.data[exhaustion.gene,]) %>%
  as.tibble() %>%
  dplyr::mutate(group = tmp.sce$treatment) %>%
  tidyr::gather(key = "gene", value = "expr", -group) %>%
  dplyr::mutate(expr = expr + rnorm(nrow(.))/1000) %>%
  ggplot(aes(gene, expr)) +
  #geom_point(aes(color = group), size = 0.01, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9)) +
  #geom_quasirandom(size = 0.1) +
  #geom_boxplot(outlier.size = -1, fill = NA) +
  geom_violin(aes(fill = factor(group, levels = c("ut","tr"))), scale = "width", trim = T, draw_quantiles = 0.5, alpha = 0.8, quantile.lwd = 1) +
  theme_classic() +
  scale_fill_manual(values = c("#3CB371", "#FF7F50")) +
  scale_color_manual(values = c("grey70", "grey70"))

tmp.sce <- sct.list[[6]]

comb.sce <- list(sct.list[[2]],p19.sce,sct.list[[4]],sct.list[[5]],sct.list[[6]])

pda <- list()
for (i in 1:length(comb.sce)) {
  pda[[i]] <- t(comb.sce[[i]]@assays$RNA@scale.data[exhaustion.gene,]) %>%
    as.tibble() %>%
    dplyr::mutate(group = comb.sce[[i]]$treatment) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise_if(is.double, mean)
}

Reduce(rbind, pda) %>%
  dplyr::mutate(patient = c(rep("P10",2),rep("P19",2),rep("P30",2),rep("P33",2),rep("P35",2))) %>%
  dplyr::mutate_if(is.double, funs((. - mean(.))/sd(.))) %>%
  tidyr::gather(key = gene, value = expr, -group, -patient) %>%
  dplyr::arrange(desc(group), patient) %>%
  dplyr::mutate(group = paste0(group, ".",patient)) %>%
  ggplot(aes(group,gene)) +
  geom_tile(aes(fill = expr)) +
  scale_fill_distiller(palette = "RdBu") +
  theme_classic()




