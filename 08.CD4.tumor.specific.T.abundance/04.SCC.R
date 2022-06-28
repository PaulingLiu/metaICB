
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

use.clones <- unique(all.tcr$clone.id)

use.clones.cells <- all.tcr %>%
  dplyr::select(patient, treatment, clone.id, cellid) %>%
  dplyr::group_by(patient, treatment, clone.id) %>%
  tidyr::nest() %>%
  dplyr::mutate(num = purrr::map_dbl(data, nrow)) %>%
  dplyr::filter(num >= 1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(cloneid = paste0("clone.",1:nrow(.)))

mcells <- unique(use.clones.cells$cloneid)

cal.genes <- c("CXCL13","CD4","CD8A")
meta.expr <- Matrix::Matrix(data = 0, nrow = length(mcells), ncol = length(cal.genes))

for (i in 1:length(mcells)) {
  tmp.cell <- use.clones.cells$data[[i]]$cellid
  print(i)
  if (length(tmp.cell) == 1) {
    meta.expr[i,] <- tcell.filt@assays$RNA@scale.data[cal.genes,tmp.cell]
  }else if(length(tmp.cell) > 1){
    meta.expr[i,] <- Matrix::rowMeans(tcell.filt@assays$RNA@scale.data[cal.genes,tmp.cell])
  }
}

pr.matr <- as.matrix(meta.expr)
rownames(pr.matr) <- mcells
colnames(pr.matr) <- cal.genes

tibble(
  CD4 = pr.matr[,"CD4"],
  CXCL13 = pr.matr[,"CXCL13"],
  CD8A = pr.matr[,"CD8A"],
  size = use.clones.cells$num,
  patient = use.clones.cells$patient,
  clone.id = use.clones.cells$clone.id,
  treatment = use.clones.cells$treatment
) -> pda

pda <- as.data.frame(pda)
rownames(pda) <- rownames(pr.matr)

pda <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/09.CD4.analysis/01.clone.expr/04.SCC.rds.gz")
tex.clones.pro <- readr::read_rds("/home/pauling/projects/07_pan.cancer.treatment/02.data/00.all.tmp.cxcl13.cd8.clone/04.SCC.rds.gz")

#------------- cutoff determination -------------#

tibble(CD8A = pda$CD4) %>%
  ggplot(aes(CD8A)) +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
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


tibble(a = pda$CXCL13) %>%
  ggplot(aes(a)) +
  geom_vline(xintercept = 0.3, linetype = "dashed") +
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

tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD4 > 0.25 & CXCL13 > 0.5, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- pda %>%
  dplyr::filter(CD4 > 0.25) %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(clone.id %in% tex.clones, "Tex", "Other"))

tex.clones.pro %>%
  dplyr::group_by(patient, treatment, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda

rep.data <- tibble(patient = unique(pro.pda$patient)) %>%
  dplyr::mutate(response = ifelse(patient %in% c("su010","su011"), "PR", "PD"))

pro.pda <- pro.pda %>% dplyr::left_join(rep.data, by = "patient") %>% dplyr::mutate(resp = paste(response, treatment, sep = "."))

pro.pda %>%
  dplyr::filter(paste0(patient,".",treatment) != "su013.post") %>%
  ggplot(aes(factor(resp, levels = c("PR.pre","PR.post","PD.pre","PD.post")), freq)) +
  geom_boxplot(aes(color = resp), outlier.size = -1, lwd = 0.6) +
  geom_point(aes(color = resp), size = 2) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black")) +
  labs(
    x = "",
    y = "CD4"
  ) +
  scale_color_manual(values = c("#009ACD", "#4F94CD", "#CD2626", "#EE799F"))

c("#009ACD", "#4F94CD", "#1874CD")
uniq.resp <- unique(pro.pda$resp)
t.test(pro.pda$freq[pro.pda$resp == uniq.resp[4]], pro.pda$freq[pro.pda$resp == uniq.resp[2]])

#-------------- barplot -----------------#
pda <- tex.clones.pro %>% dplyr::select(-cloneid, -clone)
tex.clones <- pda %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(CD8A > 0.5 & CXCL13 > 0.2, "Tex", "Other")) %>%
  dplyr::filter(clone == "Tex") %>%
  dplyr::pull(clone.id)

tex.clones.pro <- pda %>%
  dplyr::filter(CD8A > 0.5) %>%
  tibble::rownames_to_column(var = "cloneid") %>%
  as.tibble() %>%
  dplyr::mutate(clone = ifelse(clone.id %in% tex.clones, "Tex", "Other"))

tex.clones.pro %>%
  dplyr::group_by(patient, treatment, clone) %>%
  dplyr::summarise(n = sum(size)) %>%
  tidyr::spread(key = clone, value = n) %>%
  dplyr::mutate(freq = Tex/(Tex + Other)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq = ifelse(is.na(freq), 0 , freq)) -> pro.pda.cd8

pro.pda %>%
  dplyr::mutate(cd4 = pro.pda.cd8$freq) %>%
  dplyr::filter(resp %in% c("PR.pre","PD.pre","PD.post")) %>%
  ggplot(aes(freq, cd4)) +
  geom_point(aes(fill = resp, color = resp), size = 3, alpha = 0.8, shape = 21, stroke = 0.8) +
  theme_bw() +
  scale_color_manual(values = c("black","#008B8B","#B452CD")) +
  scale_fill_manual(values = c("#008B8B","#008B8B","#B452CD")) +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_vline(xintercept = 0.45, linetype = "dashed") +
  scale_y_sqrt() +
  scale_x_sqrt() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 13,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "",
    y = ""
  )
