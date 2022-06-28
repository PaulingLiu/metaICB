
gene <- "RPS20"
gene <- "GPI"
da1 <- tibble(
  CXCL13 = sce.nrp@assays$RNA@scale.data[gene,],
  label = Idents(sce.nrp)
) %>%
  dplyr::filter(label == "MANA") %>%
  dplyr::mutate(data = "NSCLC") %>%
  dplyr::select(-label)

da2 <- tibble(
  CXCL13 = sce.hn.pro@assays$RNA@scale.data[gene,]
) %>%
  dplyr::mutate(data = "HNSCC")

rbind(da1, da2) %>%
  ggplot(aes(CXCL13)) +
  geom_density(aes(fill = data), alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#009ACD", "#8968CD")) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  theme(
    #legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  labs(
    x = "RPS20"
  )

pda %>%
  ggplot(aes(type, m)) +
  geom_jitter(aes(color = type),width = 0.15, alpha = 0.8, size = 2) +
  theme_classic() +
  geom_hline(yintercept = 0.15, linetype = "dashed") +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 15,color="black"),
    axis.text = element_text(size = 12,color="black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(color="black"),
    axis.text.x = element_text(color="black")
  ) +
  scale_color_manual(values = c("#9F79EE", "#9F79EE", "#00868B")) +
  scale_y_sqrt(breaks = c(0, 0.2, 1, 2.5, 5)) +
  coord_flip() +
  labs(
    y = "CXCL13"
  )

pda %>%
  dplyr::filter(m != 0) %>%
  dplyr::filter(type != "Viral") %>%
  dplyr::count(type)

pda %>%
  dplyr::filter(m != 0) %>%
  dplyr::filter(type != "Viral") %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(n = sum(size))


sce.hn.filt <- sce.hn[,sce.hn$clone.id %in% filt.pda$trb]

sce.pro.filt <- sce.pro[, colnames(sce.pro) %in% colnames(sce.hn)]
DimPlot(sce.pro.filt, cols = comb.d3[c(12,13,7,10)], pt.size = 0.8)
FeaturePlot(sce.pro.filt, features = "CXCL13", pt.size = 0.3, cols = c("grey92", "#EE0000", "#CD2626")) + NoAxes()

a <- sce.pro.filt[,sce.pro.filt@assays$RNA@scale.data["CXCL13",] > 0.75]
b <- sce.pro.filt[,sce.pro.filt@assays$RNA@scale.data["CXCL13",] < 0.75]
d <- sce.pro.filt[,colnames(sce.pro.filt) %in% colnames(sce.hn.filt)]
e <- sce.pro.filt[,colnames(sce.pro.filt) %in% colnames(sce.hn)]
aa <- table(Idents(a))
bb <- table(Idents(b))
dd <- table(Idents(d))
ee <- table(Idents(e))

aa/sum(aa)
bb/sum(bb)
dd/sum(dd)
ee/sum(ee)

sce.pro <- sce
sce.pro <- AddModuleScore(sce.pro, features = markers.subsets, name = "Tcell_CD8_Seurat")

cal.sig <- function(.x){
  tibble(
    stem = .x$Tcell_CD8_Seurat1,
    terminal = .x$Tcell_CD8_Seurat2,
    cluster = as.character(Idents(.x))
  )
}

sig.score <- cal.sig(sce.pro)

sig.score %>%
  dplyr::filter(!(cluster %in% c("Prolif.","Prolf"))) %>%
  ggplot(aes(factor(cluster, levels = c("stem-like","Transitory","TD")), stem)) +
  geom_quasirandom(aes(color = cluster),size = 0.5) +
  geom_boxplot(outlier.size = -1, fill = NA) +
  theme_classic() +
  scale_color_manual(values = comb.d3[c(13,12,10)]) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 0, colour = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 35, hjust = 1)
  )

#-- CXCL13+ clones

filt.pda <- pda %>%
  dplyr::filter(type == "tumor") %>%
  dplyr::filter(m > 0.15)
