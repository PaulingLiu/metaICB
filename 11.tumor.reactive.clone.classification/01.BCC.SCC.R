
#----- library -----#

library(Seurat)
library(tidyverse)
library(reticulate)
library(ComplexHeatmap)
library(pROC)
library(ggROC)

#--- load data

bcc.sce <- readr::read_rds("/raid1/pauling/projects/01_data/12_NM_PD1/02.bcc.tcell.sce.rds.gz")
mda <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/02.meta.txt.gz")

all.tcr <- readr::read_tsv("/raid1/pauling/projects/01_data/12_NM_PD1/GSE123813_bcc_tcr.txt.gz")
all.tcr <- all.tcr %>%
  dplyr::filter(stringr::str_detect(cdr3s_aa, "TRB")) %>%
  dplyr::mutate(
    sep = purrr::map(
      .x = cdr3s_aa,
      .f = function(.x){
        .x <- stringr::str_replace_all(.x, "\\\t", ";")
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

#--- CD8 T cells

cd8.clones <- tex.clones.pro %>%
  dplyr::filter(treatment == "pre") %>%
  dplyr::filter(CD8A > 0.5) %>%
  dplyr::filter(patient %in% c("su001","su002","su003","su004","su009","su012"))

#--- pre-treatment samples

rownames(all.tcr) <- all.tcr$cellid
all.tcr <- all.tcr[tcell.filt$cell.id,]
tcell.filt$clone.id <- all.tcr$clone.id

bcc.sce.pre <- tcell.filt[,tcell.filt$treatment == "pre" & tcell.filt$clone.id %in% cd8.clones$clone.id]

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
  dplyr::mutate(label = ifelse(m>0.1 & m2>2, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.1 & m2<2, "nonTex", label)) %>%
  dplyr::mutate(label = ifelse(m<0.1 & label == "Other", "Other.1", label)) 

rbind(pre.data.sta, pre.data.sta.SCC) %>%
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
  cxcl13 = bcc.sce.pre@assays$RNA@scale.data["CXCL13",],
  clone = bcc.sce.pre$clone.id,
  #type = cd8.clone.scale$type,
  tex = bcc.sce.pre@assays$RNA@scale.data["ENTPD1",] + bcc.sce.pre@assays$RNA@scale.data["LAYN",] + 
    bcc.sce.pre@assays$RNA@scale.data["ITGAE",] + bcc.sce.pre@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::group_by(clone) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 5) %>%
  dplyr::mutate(label = ifelse(m>0.1 & m2>2, "Tex", "Other")) %>%
  dplyr::mutate(label = ifelse(m<0.1 & m2<2, "nonTex", label)) %>%
  dplyr::filter(label != "Other")

#--- post-treatment samples

cd8.clones.post <- tex.clones.pro %>%
  dplyr::filter(treatment == "post") %>%
  dplyr::filter(CD8A > 0.5) %>%
  dplyr::filter(patient %in% c("su001","su002","su003","su004","su009","su012"))

cd8.clone.post <- tcell.filt[,tcell.filt$treatment == "post" & tcell.filt$clone.id %in% pre.data$clone]

bcc.post.res <- tibble(
  cxcl13 = cd8.clone.post@assays$RNA@scale.data["CXCL13",],
  trb = cd8.clone.post$clone.id,
  tex = cd8.clone.post@assays$RNA@scale.data["ENTPD1",] + cd8.clone.post@assays$RNA@scale.data["LAYN",] + 
    cd8.clone.post@assays$RNA@scale.data["ITGAE",] + cd8.clone.post@assays$RNA@scale.data["BATF",]
) %>%
  dplyr::mutate(type = ifelse(trb %in% pre.data$clone[pre.data$label == "Tex"], "Tex", "nonTex")) %>%
  dplyr::group_by(trb, type) %>%
  dplyr::summarise(m = mean(cxcl13), n = n(), m2 = mean(tex)) %>%
  dplyr::filter(n > 2)

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
