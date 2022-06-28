
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#--- load data

tcell.sce <- readr::read_rds("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/tcell.rds.gz")
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

all.tcr <- all.tcr %>%
  dplyr::filter(Patient %in% c("P012","P019"))

all.tcr <- as.data.frame(all.tcr)
rownames(all.tcr) <- all.tcr$cellID

length(intersect(colnames(tcell.sce), all.tcr$cellID))
tcell.sce <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellID]
all.tcr <- all.tcr[colnames(tcell.sce),]
tcell.sce$clone.id <- all.tcr$clone.id
tcell.sce$treatment <- all.tcr$Group

#--- CD8 T cells

cd8.clones <- tex.clones.pro %>%
  dplyr::filter(group == "Pre-treatment") %>%
  dplyr::filter(CD8A > 0.5)

#--- pre-treatment samples

bcc.sce.pre <- tcell.sce[,tcell.sce$treatment == "Pre-treatment" & tcell.sce$clone.id %in% cd8.clones$clone.id]

pre.data.sta <- tibble(
  cxcl13 = bcc.sce.pre@assays$RNA@scale.data["CXCL13",],
  clone = bcc.sce.pre$clone.id,
  #type = cd8.clone.scale$type,
  tex = bcc.sce.pre@assays$RNA@scale.data["ENTPD1",] + bcc.sce.pre@assays$RNA@scale.data["LAYN",] + 
    bcc.sce.pre@assays$RNA@scale.data["ITGAE",] + bcc.sce.pre@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::mutate(label = ifelse(m>0.05 & m2>0.8, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.05 & m2<0.8, "nonTex", label)) %>%
  dplyr::mutate(label = ifelse(m<0.05 & label == "Other", "Other.1", label)) 

pre.data.sta %>%
  ggplot(aes(1+m, m2)) +
  geom_point(aes(color = label), size = 2) +
  geom_vline(xintercept = 1+0.05) +
  geom_hline(yintercept = 0.8) +
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
  scale_color_manual(values = c("#00868B", "grey85","#CD2626")) +
  scale_x_log10()

pre.data.sta %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(sum = sum(n), n = n())

pre.data <- tibble(
  cxcl13 = bcc.sce.pre@assays$RNA@scale.data["CXCL13",],
  clone = bcc.sce.pre$clone.id,
  #type = cd8.clone.scale$type,
  tex = bcc.sce.pre@assays$RNA@scale.data["ENTPD1",] + bcc.sce.pre@assays$RNA@scale.data["LAYN",] + 
    bcc.sce.pre@assays$RNA@scale.data["ITGAE",] + bcc.sce.pre@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::mutate(label = ifelse(m>0.05 & m2>0.8, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.05 & m2<0.8, "nonTex", label)) %>%
  dplyr::filter(label != "Other")

#--- post-treatment samples

cd8.clone.post <- tcell.sce[,tcell.sce$treatment == "Post-treatment" & tcell.sce$clone.id %in% pre.data$clone]

bcc.post.res <- tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post@assays$RNA@scale.data["ITGAE",] + cd8.clone.post@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::mutate(type = ifelse(trb %in% pre.data$clone[pre.data$label == "Tex"], "Tex", "nonTex")) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 1)

cd8.clone.post.filt <- cd8.clone.post[,cd8.clone.post$clone.id %in% bcc.post.res$trb]

bcc.post.res %>%
  ggplot(aes(m, m2)) +
  geom_point(aes(color = type), size= 3) +
  geom_vline(xintercept = 0.5) +
  geom_hline(yintercept = 2) +
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

#--- classification

cd8.clone.post.filt$type <- ifelse(cd8.clone.post.filt$clone.id %in% pre.data$clone[pre.data$label == "Tex"], "Tex", "nonTex")

rda.bcc <- tibble(
  type = cd8.clone.post.filt$type,
  trb = cd8.clone.post.filt$clone.id,
  CXCL13 = cd8.clone.post.filt@assays$RNA@scale.data["CXCL13",],CTLA4 = cd8.clone.post.filt@assays$RNA@scale.data["CTLA4",],
  PDCD1 = cd8.clone.post.filt@assays$RNA@scale.data["PDCD1",],ENTPD1 = cd8.clone.post.filt@assays$RNA@scale.data["ENTPD1",],
  LAG3 = cd8.clone.post.filt@assays$RNA@scale.data["LAG3",],LAYN = cd8.clone.post.filt@assays$RNA@scale.data["LAYN",],
  TIGIT = cd8.clone.post.filt@assays$RNA@scale.data["TIGIT",],BAG3 = cd8.clone.post.filt@assays$RNA@scale.data["BAG3",],
  BATF = cd8.clone.post.filt@assays$RNA@scale.data["BATF",],HAVCR2 = cd8.clone.post.filt@assays$RNA@scale.data["HAVCR2",],
  TNFRSF9 = cd8.clone.post.filt@assays$RNA@scale.data["TNFRSF9",],GZMB = cd8.clone.post.filt@assays$RNA@scale.data["GZMB",],
  CCL5 = cd8.clone.post.filt@assays$RNA@scale.data["CCL5",], 
  tex = cd8.clone.post.filt@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post.filt@assays$RNA@scale.data["CTLA4",] + 
    cd8.clone.post.filt@assays$RNA@scale.data["LAYN",] + cd8.clone.post.filt@assays$RNA@scale.data["HAVCR2",]
)


rda <- rda.bcc %>%
  dplyr::group_by(type, trb) %>%
  dplyr::summarise(
    CXCL13 = mean(CXCL13),CTLA4 = mean(CTLA4),PDCD1 = mean(PDCD1),
    ENTPD1 = mean(ENTPD1),LAG3 = mean(LAG3),LAYN = mean(LAYN),
    TIGIT = mean(TIGIT),BAG3 = mean(BAG3),BATF = mean(BATF),
    HAVCR2 = mean(HAVCR2),TNFRSF9 = mean(TNFRSF9),GZMB = mean(GZMB),
    CCL5 = mean(CCL5),n=n(), tex = mean(tex),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(response = ifelse(type == "Tex", 1, 0))

rda <- rda.bcc %>%
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
