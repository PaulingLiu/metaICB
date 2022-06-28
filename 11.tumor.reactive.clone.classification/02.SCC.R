
#----- library -----#

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#--- load data

scc.sce <- readr::read_rds("/raid1/pauling/projects/01_data/12_NM_PD1/04.scc.tcell.rds.gz")
scc.sce <- NormalizeData(scc.sce)
scc.sce <- ScaleData(scc.sce, features = rownames(scc.sce), do.center = F, do.scale = F)
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

tcell.filt <- scc.sce[,colnames(scc.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))
all.tcr <- all.tcr %>% dplyr::inner_join(mda[,c(1,2,3)], by = c("cellid" = "cell.id"))

#--- CD8 T cells

cd8.clones <- tex.clones.pro %>%
  dplyr::filter(treatment == "pre") %>%
  dplyr::filter(CD8A > 0.5) %>%
  dplyr::filter(patient %in% c("su010","su011"))

#--- pre-treatment samples

rownames(all.tcr) <- all.tcr$cellid
all.tcr <- all.tcr[tcell.filt$cell.id,]
tcell.filt$clone.id <- all.tcr$clone.id

scc.sce.pre <- tcell.filt[,tcell.filt$treatment == "pre" & tcell.filt$clone.id %in% cd8.clones$clone.id]

pre.data.sta.SCC <- tibble(
  cxcl13 = scc.sce.pre@assays$RNA@scale.data["CXCL13",],
  clone = scc.sce.pre$clone.id,
  #type = cd8.clone.scale$type,
  tex = scc.sce.pre@assays$RNA@scale.data["ENTPD1",] + scc.sce.pre@assays$RNA@scale.data["LAYN",] + 
    scc.sce.pre@assays$RNA@scale.data["ITGAE",] + scc.sce.pre@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::mutate(label = ifelse(m>0.1 & m2>2, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.1 & m2<2, "nonTex", label)) %>%
  dplyr::mutate(label = ifelse(m<0.1 & label == "Other", "Other.1", label)) 

pre.data.sta.SCC %>%
  ggplot(aes(1+m, m2)) +
  geom_point(aes(color = label)) +
  geom_vline(xintercept = 1+0.1) +
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
  scale_color_manual(values = c("#00868B", "grey85","grey85","#CD2626")) +
  scale_x_log10()


pre.data <- tibble(
  cxcl13 = scc.sce.pre@assays$RNA@scale.data["CXCL13",],
  clone = scc.sce.pre$clone.id,
  #type = cd8.clone.scale$type,
  tex = scc.sce.pre@assays$RNA@scale.data["ENTPD1",] + scc.sce.pre@assays$RNA@scale.data["LAYN",] + 
    scc.sce.pre@assays$RNA@scale.data["ITGAE",] + scc.sce.pre@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::mutate(label = ifelse(m>0.1 & m2>2, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.1 & m2<2, "nonTex", label)) %>%
  dplyr::filter(label != "Other")

#--- post-treatment samples

cd8.clone.post.SCC <- tcell.filt[,tcell.filt$treatment == "post" & tcell.filt$clone.id %in% pre.data$clone]

bcc.post.res.scc <- tibble(
  cxcl13 = cd8.clone.post.SCC@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post.SCC$clone.id,
  tex = cd8.clone.post.SCC@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post.SCC@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post.SCC@assays$RNA@scale.data["ITGAE",] + cd8.clone.post.SCC@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::mutate(type = ifelse(trb %in% pre.data$clone[pre.data$label == "Tex"], "Tex", "nonTex")) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 2)

cd8.clone.post.SCC.filt <- cd8.clone.post.SCC[,cd8.clone.post.SCC$clone.id %in% bcc.post.res.scc$trb]

bcc.post.res.scc %>%
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

#--- combination

pre.data.sta %>%
  dplyr::bind_rows(pre.data.sta.SCC) %>%
  ggplot(aes(1+m, m2)) +
  geom_point(aes(color = label), size = 2) +
  geom_vline(xintercept = 1+0.1) +
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
  scale_color_manual(values = c("#00868B", "grey85","grey85","#CD2626")) +
  scale_x_log10()

pre.data.sta %>%
  dplyr::bind_rows(pre.data.sta.SCC) %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(sum = sum(n), n = n())


#----- AUC score for each gene -----#


cd8.clone.post.filt$type <- ifelse(cd8.clone.post.filt$clone.id %in% pre.data$clone[pre.data$label == "Tex"], "Tex", "nonTex")
cd8.clone.post.SCC.filt$type <- ifelse(cd8.clone.post.SCC.filt$clone.id %in% bcc.post.res.scc$trb[bcc.post.res.scc$type == "Tex"], "Tex", "nonTex")

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

rda.scc <- tibble(
  type = cd8.clone.post.SCC.filt$type,
  trb = cd8.clone.post.SCC.filt$clone.id,
  CXCL13 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["CXCL13",],CTLA4 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["CTLA4",],
  PDCD1 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["PDCD1",],ENTPD1 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["ENTPD1",],
  LAG3 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["LAG3",],LAYN = cd8.clone.post.SCC.filt@assays$RNA@scale.data["LAYN",],
  TIGIT = cd8.clone.post.SCC.filt@assays$RNA@scale.data["TIGIT",],BAG3 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["BAG3",],
  BATF = cd8.clone.post.SCC.filt@assays$RNA@scale.data["BATF",],HAVCR2 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["HAVCR2",],
  TNFRSF9 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["TNFRSF9",],GZMB = cd8.clone.post.SCC.filt@assays$RNA@scale.data["GZMB",],
  CCL5 = cd8.clone.post.SCC.filt@assays$RNA@scale.data["CCL5",], 
  tex = cd8.clone.post.SCC.filt@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post.SCC.filt@assays$RNA@scale.data["CTLA4",] + 
    cd8.clone.post.SCC.filt@assays$RNA@scale.data["LAYN",] + cd8.clone.post.SCC.filt@assays$RNA@scale.data["HAVCR2",]
)

rda <- rbind(rda.scc, rda.bcc) %>%
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

rda <- rbind(rda.scc, rda.bcc) %>%
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
