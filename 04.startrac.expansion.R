
#---------- library -----------#

library(Seurat)
library(tidyverse)
library(reticulate)
library(RColorBrewer)

#---------- load data -----------#

tex.clones.pro1 <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/01.NSCLC1.rds.gz")
tex.clones.pro1 <- tex.clones.pro1 %>% dplyr::filter(clone == "Tex")

tibble(
  patient = unique(tex.clones.pro1$patient),
  response = c('PR',"PD","PR","pre","PR","pre","PR","PD","pre","PR","pre","PR","pre","PR","pre","PR","pre","PD","PD","PD")
) -> mda

pda <- tex.clones.pro1 %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient) %>%
  dplyr::right_join(mda, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n))

View(pda)

pda %>%
  dplyr::rename(size = n) %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("pre","PR","PD")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  #geom_crossbar(aes(ymin = mean+1, ymax = mean+1), width = 0.6, lwd = 0.5) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_fill_manual(values = c("#009ACD", "#CD2626", "#EE799F")) +
  scale_color_manual(values = c("#009ACD", "#CD2626", "#EE799F"))

wilcox.test(pda$n[pda$response == "PD"], pda$n[pda$response == "PR"])
t.test(pda$n[pda$response == "PD"], pda$n[pda$response == "PR"])

#---------- NSCLC2 -----------#

tex.clones.pro2 <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/02.NSCLC2.new.rds.gz")

rep.data <- tibble(
  patient = c("MD01-004","MD01-005","MD01-010","MD01-019","MD01-024","MD043-003","MD043-006",
              "MD043-008","MD043-011","NY016-007","NY016-014","NY016-015","NY016-021","NY016-022","NY016-025"),
  response = c('non-MPR','MPR','MPR','non-MPR','non-MPR','MPR','non-MPR','MPR','non-MPR','non-MPR',
               'non-MPR','non-MPR','non-MPR','MPR','MPR')
)

pda <- tex.clones.pro2 %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient) %>%
  tidyr::separate(patient, c('patient',"tumor","num"), sep = "\\.") %>%
  dplyr::group_by(patient) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::right_join(rep.data, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n))

pda %>%
  dplyr::rename(size = n) %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("MPR","non-MPR")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#CD2626","#009ACD"))

wilcox.test(pda$n[pda$response == "MPR"], pda$n[pda$response != "MPR"])

#---------- BCC -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/03.BCC.rds.gz")

rep.data <- tibble(patient = unique(tex.clones.pro$patient)) %>%
  dplyr::mutate(response = ifelse(patient %in% c("su001","su002","su003","su004","su009","su012"), "PR", "PD"))

pda <- tex.clones.pro %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 1) %>%
  dplyr::count(patient, treatment) %>%
  dplyr::group_by(patient, treatment) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::right_join(rep.data, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n)) %>%
  dplyr::ungroup()

pda$treatment[16] <- "pre"
pda$treatment[is.na(pda$treatment)] = "pre"

pda <- pda %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(resp = paste0(response, ".", treatment)) %>%
  #dplyr::select(patient, n, resp) %>%
  tidyr::spread(key = treatment, value = n) %>%
  dplyr::mutate_if(is.double, .funs = function(.x){ifelse(is.na(.x), 0, .x)}) %>%
  tidyr::gather(key = "treatment", value = "n", -patient, -response)

pda <- pda %>%
  dplyr::rename(size = n) %>%
  dplyr::mutate(response = paste0(response, ".", treatment))

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("PR.pre","PR.post","PD.pre","PD.post")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#009ACD", "#4F94CD", "#CD2626", "#EE799F"))

uniq.resp <- unique(pda$response)
wilcox.test(pda$size[pda$response == uniq.resp[4]], pda$size[pda$response == uniq.resp[1]])
t.test(pda$size[pda$response == uniq.resp[4]], pda$size[pda$response == uniq.resp[1]])

#---------- SCC -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/04.SCC.rds.gz")

rep.data <- tibble(patient = unique(tex.clones.pro$patient)) %>%
  dplyr::mutate(response = ifelse(patient %in% c("su010","su011"), "PR", "PD"))

pda <- tex.clones.pro %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient, treatment) %>%
  dplyr::group_by(patient, treatment) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::right_join(rep.data, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n)) %>%
  dplyr::ungroup()


pda <- pda %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(resp = paste0(response, ".", treatment)) %>%
  #dplyr::select(patient, n, resp) %>%
  tidyr::spread(key = treatment, value = n) %>%
  dplyr::mutate_if(is.double, .funs = function(.x){ifelse(is.na(.x), 0, .x)}) %>%
  tidyr::gather(key = "treatment", value = "n", -patient, -response)

pda <- pda %>%
  dplyr::rename(size = n) %>%
  dplyr::mutate(response = paste0(response, ".", treatment))

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("PR.pre","PR.post","PD.pre","PD.post")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#009ACD", "#4F94CD", "#CD2626", "#EE799F"))

uniq.resp <- unique(pda$response)
uniq.resp
wilcox.test(pda$size[pda$response == uniq.resp[1]], pda$size[pda$response == uniq.resp[4]])

#---------- breast cancer batch 1 -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/05.breast.cancer.batch1.rds.gz")

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1863-counts_cells_cohort1-T_cell.metadata.csv", col_names = F)
meta <- meta %>% 
  dplyr::select(X1, X4, X5, X6, X11) %>%
  dplyr::rename(
    cellid = X1,
    patient = X4,
    treatment = X5,
    response = X6,
    batch = X11
  ) %>%
  dplyr::filter(response %in% c("E","NE")) %>%
  dplyr::distinct(patient, response)

pda <- tex.clones.pro %>%
  tidyr::separate(sample, c("biok","num","treatment"), sep = "\\_", remove = F) %>%
  dplyr::mutate(patient = paste(biok,num, sep = "_")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient,treatment) %>%
  dplyr::group_by(patient,treatment) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::right_join(meta, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n)) %>%
  dplyr::mutate(treatment = ifelse(is.na(treatment), "Pre", treatment)) %>%
  dplyr::ungroup()

pda <- pda %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(resp = paste0(response, ".", treatment)) %>%
  #dplyr::select(patient, n, resp) %>%
  tidyr::spread(key = treatment, value = n) %>%
  dplyr::mutate_if(is.double, .funs = function(.x){ifelse(is.na(.x), 0, .x)}) %>%
  tidyr::gather(key = "treatment", value = "n", -patient, -response)

pda <- pda %>%
  dplyr::rename(size = n) %>%
  dplyr::mutate(response = paste0(response, ".", treatment))

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("E.Pre","E.On","NE.Pre","NE.On")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#CD2626", "#EE799F", "#009ACD", "#4F94CD"))

uniq.resp <- unique(pda$response)
uniq.resp
wilcox.test(pda$size[pda$response == uniq.resp[1]], pda$size[pda$response == uniq.resp[4]])

#---------- breast cancer batch 2 -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/06.breast.cancer.batch2.rds.gz")

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/04.breast.pdl1/1867-counts_cells_cohort2-T_cell.metadata.csv", col_names = F)
meta <- meta %>% 
  dplyr::select(X1, X4, X5, X6, X11) %>%
  dplyr::rename(
    cellid = X1,
    patient = X4,
    treatment = X5,
    response = X6,
    batch = X11
  ) %>%
  dplyr::filter(response %in% c("E","NE")) %>%
  dplyr::distinct(patient, response)

pda <- tex.clones.pro %>%
  tidyr::separate(sample, c("biok","num","treatment"), sep = "\\_", remove = F) %>%
  dplyr::mutate(patient = paste(biok,num, sep = "_")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient,treatment) %>%
  dplyr::group_by(patient,treatment) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::right_join(meta, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n)) %>%
  dplyr::mutate(treatment = ifelse(is.na(treatment), "Pre", treatment)) %>%
  dplyr::ungroup()

pda <- pda %>%
  dplyr::ungroup() %>%
  #dplyr::mutate(resp = paste0(response, ".", treatment)) %>%
  #dplyr::select(patient, n, resp) %>%
  tidyr::spread(key = treatment, value = n) %>%
  dplyr::mutate_if(is.double, .funs = function(.x){ifelse(is.na(.x), 0, .x)}) %>%
  tidyr::gather(key = "treatment", value = "n", -patient, -response)

pda <- pda %>%
  dplyr::rename(size = n) %>%
  dplyr::mutate(response = paste0(response, ".", treatment))

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("E.Pre","E.On","NE.Pre","NE.On")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#CD2626", "#EE799F", "#009ACD", "#4F94CD"))

uniq.resp <- unique(pda$response)
uniq.resp
wilcox.test(pda$size[pda$response == uniq.resp[2]], pda$size[pda$response == uniq.resp[3]])

#---------- breast cancer yy -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/07.breast.yy.rds.gz")

meta <- readr::read_csv("/raid1/pauling/projects/07_pan.cancer.treatment/01.scRNAseq.collection/01.breast.yy/Stable2-new-CCR3.csv")
meta <- tibble(
  patient = meta$Patient,
  response = meta$Efficacy,
  group = meta$Group
) %>%
  dplyr::distinct(patient, response) %>%
  dplyr::filter(patient %in% tex.clones.pro$patient)

pda <- tex.clones.pro %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient,group) %>%
  dplyr::group_by(patient, group) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::ungroup() %>%
  dplyr::right_join(meta, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n))

pda <- pda %>%
  dplyr::rename(size = n)

pda$group[is.na(pda$group)] <- "Pre-treatment"

pda.add1 <- tibble(patient = c("P005","P017"), group = c("Pre-treatment","Post-treatment"), size = c(0,0), response = c("SD","SD"))
pda.add1 <- tibble(patient = c("P005","P017","P002"), group = c("Pre-treatment","Post-treatment","Pre-treatment"), size = c(0,0,0), response = c("SD","SD","SD"))

pda <- rbind(pda.add1, pda) %>%
  dplyr::mutate(response = paste0(response, ".", group)) %>%
  dplyr::arrange(patient, group)

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("PR.Pre-treatment","PR.Post-treatment","SD.Pre-treatment","SD.Post-treatment")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#CD2626", "#EE799F","#009ACD", "#4F94CD"))

uniq.resp <- unique(pda$response)
uniq.resp
wilcox.test(pda$size[pda$response == uniq.resp[4]], pda$size[pda$response == uniq.resp[1]])

#---------- renal cell -----------#

tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/08.renal.1.rds.gz")

meta <- mda <- tibble(
  patient = c("t2","t3","t4"),
  response = c("non.rp","rp","rp")
)

pda <- tex.clones.pro %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::filter(size > 19) %>%
  dplyr::count(patient,sample) %>%
  dplyr::group_by(patient, sample) %>%
  dplyr::summarise(n = mean(n)) %>%
  dplyr::ungroup() %>%
  dplyr::right_join(meta, by = "patient") %>%
  dplyr::mutate(n = ifelse(is.na(n),0, n))

pda <- pda %>%
  dplyr::rename(size = n)

View(pda)

pda %>%
  dplyr::group_by(response) %>%
  dplyr::summarise(mean = mean(size), n = length(size), sd = sd(size)) %>%
  dplyr::mutate(ymax = mean + sd/sqrt(n), ymin = mean - sd/sqrt(n)) %>%
  ggplot(aes(factor(response, levels = c("rp","non.rp")), mean)) +
  geom_pointrange(aes(ymin = ymin, ymax = ymax, color = response), size = 1) +
  theme_classic() +
  labs(
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#009ACD","#CD2626"))

uniq.resp <- unique(pda$response)
uniq.resp
t.test(pda$size[pda$response == uniq.resp[1]], pda$size[pda$response == uniq.resp[2]])
