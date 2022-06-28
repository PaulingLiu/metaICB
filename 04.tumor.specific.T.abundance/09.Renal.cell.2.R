
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/03.RCC.data.2/01.tcell.rds.gz")
tibble(
  CXCL13 = sce@assays$RNA@scale.data["CXCL13",],
  patient = sce$patient,
  clusters = Idents(sce)
) %>%
  dplyr::mutate(type = ifelse(clusters %in% c(0,9), "CD8", "CD4")) -> pda

#------------- cutoff determination -------------#

tibble(CD8A = sce@assays$RNA@scale.data["CD8A",]) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_density(lwd = 1, color = "#EE6363") +
  theme_classic() +
  labs(
    x = "CD8A expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

tibble(a = sce@assays$RNA@scale.data["CXCL13",]) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 0.8, linetype = "dashed") +
  geom_density(lwd = 1, color = "#008B8B") +
  theme_classic() +
  labs(
    x = "CXCL13 expression",
    y = "Density"
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  )

#------------- Tex cell frequency -------------#

tex.clones.pro <- pda %>%
  dplyr::mutate(clone = ifelse(type == "CD8" & CXCL13 > 0.8, "Tex", "Other"))

tex.clones.pro %>%
  dplyr::group_by(patient, clone) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Tex = ifelse(is.na(Tex), 0 , Tex)) %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

pro.pda <- pro.pda %>%
  dplyr::filter(patient %in% c("p55","p913","p906","p915")) %>%
  dplyr::mutate(response = c("rp","non.rp","non.rp","rp"))

pro.pda %>%
  ggplot(aes(factor(response, levels = c("rp","non.rp")), freq)) +
  geom_boxplot(aes(color = response), outlier.size = -1, lwd = 0.6) +
  geom_point(aes(color = response), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "Frequency of CD8+ tumour-reactive\ncells in all T cells"
  ) +
  scale_color_manual(values = c("#009ACD","#CD2626"))

uniq.resp <- unique(pro.pda$response)
wilcox.test(pro.pda$freq[pro.pda$response == uniq.resp[1]], pro.pda$freq[pro.pda$response == uniq.resp[2]])
