
#--- library

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)
library(scibet)
library(ggsci)

#--- load data

lung1.sce <- readr::read_rds("projects/07_pan.cancer.treatment/")
bcc.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/05.bcc/01.tumor.reactive.cells.rds.gz")
scc.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/06.scc/01.tumor.reactive.cells.rds.gz")
breast1.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/01.tumor.reactive.cells.batch1.rds.gz")
breast2.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/02.tumor.reactive.cells.batch2.rds.gz")
breast3.sce <- readr::read_rds("projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/02.tumor.reactive.cells.batch1.rds.gz")

#--- lung cancer data

tcell.sce <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/02_clustering/revision/00.re.clustering/02.tcell.sce.rds.gz")
all.tcr <- readr::read_rds("/raid1/pauling/projects/04_lung_immune_therapy/04_TCR_analysis/05.scTCR.processed/merge.tcr.new.rds.gz")

tcell.sce <- tcell.sce[,tcell.sce$patient %in% c("P1","P10","P13","P19","P30","P33","P35")]
tcell.filt <- tcell.sce[,colnames(tcell.sce) %in% all.tcr$cellid]
all.tcr <- all.tcr %>% dplyr::filter(cellid %in% colnames(tcell.filt))
pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/01.NSCLC1.rds.gz")
pda <- pda %>% dplyr::filter(clone == "Tex")

all.tcr.pro <- all.tcr %>% dplyr::filter(clone.id %in% pda$clone.id)
sce.pro <- tcell.filt[,colnames(tcell.filt) %in% all.tcr.pro$cellid]

all.tcr.pro <- as.data.frame(all.tcr.pro)
rownames(all.tcr.pro) <- all.tcr.pro$cellid
all.tcr.pro <- all.tcr.pro[colnames(sce.pro),]
sce.pro$trb <- all.tcr.pro$clone.id
lung1.sce <- sce.pro[,sce.pro$num %in% c(0,1)]

#--- data processing

breast1.sce$trb <- stringr::str_remove(breast1.sce$trb, "_Pre")
breast1.sce$trb <- stringr::str_remove(breast1.sce$trb, "_On")

breast2.sce$trb <- stringr::str_remove(breast2.sce$trb, "_Pre")
breast2.sce$trb <- stringr::str_remove(breast2.sce$trb, "_On")

#--- calculation of clonal revival vs. clonal replacement

# pre + post > 19 (cell counts)
# post > 10 (cell counts)

replac.prop <- function(.x, cancer = ""){
  .x$treatment[.x$treatment %in% c("Pre-treatment","Pre","ut")] <- "pre"
  .x$treatment[.x$treatment %in% c("Post-treatment","On","tr")] <- "post"
  
  tmp.sta <- table(.x$patient, .x$treatment) %>%
    as.data.frame.array() %>%
    tibble::rownames_to_column(var = "patient") %>%
    dplyr::filter(pre + post > 19) %>%
    dplyr::filter(post > 0)
  
  used.patients <- tmp.sta$patient
  post.cell.counts <- tmp.sta$post
  
  .x <- .x[,.x$patient %in% used.patients]
  pda <- tibble(
    clonotype = .x$trb,
    patient = .x$patient,
    treatment = .x$treatment
  ) %>%
    dplyr::count(treatment, patient, clonotype) %>%
    tidyr::spread(key = "treatment", value = "n") %>%
    dplyr::mutate(pre = replace_na(pre, 0), post = replace_na(post, 0))
  
  res <- list()
  for (i in 1:length(used.patients)){
    tmp.pda <- pda %>%
      dplyr::filter(pre > 1 | post > 1) %>%
      dplyr::mutate(group = ifelse(pre == 0, "New", "Other")) %>%
      dplyr::filter(patient == used.patients[i])
    
    n.pre <- tmp.pda %>% dplyr::filter(pre > 1) %>% nrow()
    n.new <- tmp.pda %>% dplyr::filter(post > 1) %>% dplyr::filter(group == "New") %>% nrow()
    n.share <- tmp.pda %>% dplyr::filter(post > 1) %>% dplyr::filter(group == "Other") %>% nrow()
    res[[i]] <- tibble(n.pre = n.pre, n.new = n.new, n.share = n.share, sample = used.patients[i], prop = n.new/(n.new+n.share))
  }
  
  Reduce(rbind, res) %>%
    dplyr::mutate(cancer = cancer) %>%
    dplyr::mutate(counts = post.cell.counts)
}

return.all <- function(.x, cancer = ""){
  .x$treatment[.x$treatment %in% c("Pre-treatment","Pre","ut")] <- "pre"
  .x$treatment[.x$treatment %in% c("Post-treatment","On","tr")] <- "post"
  
  tmp.sta <- table(.x$patient, .x$treatment) %>%
    as.data.frame.array() %>%
    tibble::rownames_to_column(var = "patient") %>%
    dplyr::filter(pre + post > 19) %>%
    dplyr::filter(post > 0)
  
  used.patients <- tmp.sta$patient
  post.cell.counts <- tmp.sta$post
  
  .x <- .x[,.x$patient %in% used.patients]
  pda <- tibble(
    clonotype = .x$trb,
    patient = .x$patient,
    treatment = .x$treatment
  ) %>%
    dplyr::count(treatment, patient, clonotype) %>%
    tidyr::spread(key = "treatment", value = "n") %>%
    dplyr::mutate(pre = replace_na(pre, 0), post = replace_na(post, 0))
  
  res <- list()
  for (i in 1:length(used.patients)){
    tmp.pda <- pda %>%
      dplyr::filter(pre > 1 | post > 1) %>%
      dplyr::mutate(group = ifelse(pre == 0, "New", "Other")) %>%
      dplyr::filter(patient == used.patients[i])
    
    n.pre <- tmp.pda %>% dplyr::filter(pre > 1) %>% nrow()
    n.new <- tmp.pda %>% dplyr::filter(post > 1) %>% dplyr::filter(group == "New") %>% nrow()
    n.share <- tmp.pda %>% dplyr::filter(post > 1) %>% dplyr::filter(group == "Other") %>% nrow()
    res[[i]] <- tmp.pda
  }
  return(res)
}

lung.res <- replac.prop(lung1.sce, "01.lung")
bcc.res <- replac.prop(bcc.sce, "02.bcc")
scc.res <- replac.prop(scc.sce, "03.scc")
bc1.res <- replac.prop(breast1.sce, "04.breast1")
bc2.res <- replac.prop(breast2.sce, "05.breast2")
bc3.res <- replac.prop(breast3.sce, "06.breast3")

pda <- rbind(lung.res, bcc.res, scc.res, bc1.res, bc2.res, bc3.res) %>% dplyr::arrange(desc(cancer))
  
pda %>%
  ggplot(aes(0.8+n.pre, prop)) +
  geom_smooth(method = "lm", color = "black", se = F, linetype = "dashed") +
  geom_point(aes(color = cancer), stroke = 1, size = 2.5, alpha = 0.9) +
  ylim(0,1) +
  theme_classic() +
  scale_x_log10() +
  scale_color_manual(values = comb.d3[c(1,3:7,10,12)]) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, colour = "black")
  ) +
  labs(
    x = "Number of pre-treatment clonotypes",
    y = "Proportion of new clonotypes\nafter treatment"
  )

rank.level <- c('05.breast2',"04.breast1","02.bcc","01.lung","03.scc","06.breast3")
pda %>%
  dplyr::filter(!is.na(prop)) %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(mean = mean(prop), n = length(prop), sd = sd(prop)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(cancer, levels = rank.level), mean)) +
  geom_col(aes(fill = cancer, color = cancer), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_jitter(aes(cancer, prop, colour = cancer), size = 2, data = pda, width = 0.1) +
  geom_errorbar(aes(ymin = mean, ymax = ymax), width = 0.35, lwd = 0.6) +
  theme_classic() +
  scale_fill_manual(values = comb.d3[c(1,3:7,10,12)]) +
  scale_color_manual(values = comb.d3[c(1,3:7,10,12)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 35, hjust = 1)
  )

prop.aov <- aov(prop~cancer, data=pda)
summary(prop.aov)


#--- breast cancer

bda <- rbind(bc1.res, bc2.res, bc3.res)
bda <- bda %>% dplyr::filter(n.pre > 10)

bda %>%
  ggplot(aes(cancer, n.pre)) +
  geom_boxplot(aes(color = cancer), outlier.size = -1, lwd = 0.7) +
  geom_jitter(aes(color = cancer), width = 0.15, size = 2) +
  theme_classic() +
  scale_color_manual(values = comb.d3[c(5,7)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, colour = "black")
  ) +
  labs(
    x = "",
    y = "Baseline clonotype counts"
  ) +
  scale_y_log10(limits = c(10,300))

bda %>%
  ggplot(aes(cancer, prop)) +
  geom_boxplot(aes(color = cancer), outlier.size = -1, lwd = 0.7) +
  geom_jitter(aes(color = cancer), width = 0.15, size = 2) +
  theme_classic() +
  scale_color_manual(values = comb.d3[c(5,7)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, colour = "black")
  ) +
  labs(
    x = "",
    y = "Baseline clonotype counts"
  ) +
  ylim(0,1)

wilcox.test(bda$prop[1:6], bda$prop[7:9])

#--- lung cancer

lung.res %>%
  dplyr::mutate(time = c(42, 61, 94, 62, 57, 80, 151)) %>%
  ggplot(aes(time, prop)) +
  geom_point(size = 1.8) +
  theme_classic() +
  geom_smooth(method = "lm", se = F, lwd = 0.8, color = "grey65", linetype = "dashed") +
  scale_x_log10()  +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, color = "black")
  )

a <- lung.res %>%
  dplyr::mutate(time = c(42, 61, 94, 62, 57, 80, 151)) 

cor.test(log(1+a$n.pre), a$prop)

#--- calculation of clonal revival vs. clonal replacement

# Cell-based calculation

# pre + post > 19 (cell counts)
# post > 10 (cell counts)

replac.cell.prop <- function(.x, cancer = ""){
  .x$treatment[.x$treatment %in% c("Pre-treatment","Pre","ut")] <- "pre"
  .x$treatment[.x$treatment %in% c("Post-treatment","On","tr")] <- "post"
  
  tmp.sta <- table(.x$patient, .x$treatment) %>%
    as.data.frame.array() %>%
    tibble::rownames_to_column(var = "patient") %>%
    dplyr::filter(pre + post > 19) %>%
    dplyr::filter(post > 0)
  
  used.patients <- tmp.sta$patient
  post.cell.counts <- tmp.sta$post
  
  .x <- .x[,.x$patient %in% used.patients]
  pda <- tibble(
    clonotype = .x$trb,
    patient = .x$patient,
    treatment = .x$treatment
  ) %>%
    dplyr::count(treatment, patient, clonotype) %>%
    tidyr::spread(key = "treatment", value = "n") %>%
    dplyr::mutate(pre = replace_na(pre, 0), post = replace_na(post, 0))
  
  res <- list()
  for (i in 1:length(used.patients)){
    tmp.pda <- pda %>%
      dplyr::filter(pre > 0 | post > 0) %>%
      dplyr::mutate(group = ifelse(pre == 0, "New", "Other")) %>%
      dplyr::filter(patient == used.patients[i])
    
    n.pre <- tmp.pda %>% dplyr::filter(pre > 0) %>% dplyr::pull(pre) %>% sum()
    n.new <- tmp.pda %>% dplyr::filter(post > 0) %>% dplyr::filter(group == "New") %>% dplyr::pull(post) %>% sum()
    n.share <- tmp.pda %>% dplyr::filter(post > 0) %>% dplyr::filter(group == "Other") %>% dplyr::pull(post) %>% sum()
    res[[i]] <- tibble(n.pre = n.pre, n.new = n.new, n.share = n.share, sample = used.patients[i], prop = n.new/(n.new+n.share))
  }
  
  Reduce(rbind, res) %>%
    dplyr::mutate(cancer = cancer) %>%
    dplyr::mutate(counts = post.cell.counts)
}

lung.cell.res <- replac.cell.prop(lung1.sce, "01.lung")
bcc.cell.res <- replac.cell.prop(bcc.sce, "02.bcc")
scc.cell.res <- replac.cell.prop(scc.sce, "03.scc")
bc1.cell.res <- replac.cell.prop(breast1.sce, "04.breast1")
bc2.cell.res <- replac.cell.prop(breast2.sce, "05.breast2")
bc3.cell.res <- replac.cell.prop(breast3.sce, "06.breast3")

pda <- rbind(lung.cell.res, bcc.cell.res, scc.cell.res, bc1.cell.res, bc2.cell.res, bc3.cell.res) %>% dplyr::arrange(desc(cancer))

pda %>%
  ggplot(aes(0.8+n.pre, prop)) +
  geom_smooth(method = "lm", color = "black", se = F, linetype = "dashed") +
  geom_point(aes(color = cancer), stroke = 1, alpha = 0.9) +
  ylim(0,1) +
  theme_classic() +
  scale_x_log10() +
  scale_color_manual(values = comb.d3[c(1,3:7,10,12)]) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 12, colour = "black")
  ) +
  labs(
    x = "Number of pre-treatment clonotypes",
    y = "Proportion of new clonotypes\nafter treatment"
  )

rank.level <- c('05.breast2',"04.breast1","02.bcc","01.lung","03.scc","06.breast3")
pda %>%
  dplyr::filter(!is.na(prop)) %>%
  dplyr::group_by(cancer) %>%
  dplyr::summarise(mean = mean(prop), n = length(prop), sd = sd(prop)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(cancer, levels = rank.level), mean)) +
  geom_col(aes(fill = cancer, color = cancer), alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_jitter(aes(cancer, prop, colour = cancer), size = 2, data = pda, width = 0.1) +
  geom_errorbar(aes(ymin = mean, ymax = ymax), width = 0.35, lwd = 0.5) +
  theme_classic() +
  scale_fill_manual(values = comb.d3[c(1,3:7,10,12)]) +
  scale_color_manual(values = comb.d3[c(1,3:7,10,12)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 35, hjust = 1)
  )

prop.aov <- aov(prop~cancer, data=pda)
summary(prop.aov)
