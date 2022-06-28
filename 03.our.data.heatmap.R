

#--- heatmap for clones from pre- and post-ICB responsive tumors

cd8.clone.post.filt <- ScaleData(cd8.clone.post, do.center = T, do.scale = T)

tibble(cell = colnames(cd8.clone.post.filt), type = cd8.clone.post.filt$type) %>%
  dplyr::arrange(type) %>%
  dplyr::mutate(cellid = stringr::str_sub(cell, -18, -1)) %>%
  dplyr::arrange(desc(type), desc(cellid)) -> tmp.rank

use.genes <- c("CXCL13","CTLA4","PDCD1","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3")
tmp.matr <- cd8.clone.post.filt@assays$RNA@scale.data[use.genes,tmp.rank$cell]


colnames(tmp.matr) <- NULL

tmp.matr[tmp.matr < -2] <- -2
tmp.matr[tmp.matr > 3.5] <- 3.5

ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("Tex" = "#CD853F","nonTex" = "#009ACD")),
  simple_anno_size = unit(0.4, "cm")
)

ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F, 
                        col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                        top_annotation = ha)


#--- heatmap for clones from pre- and post-ICB responsive tumors

cd8.clone.pre.scale <- ScaleData(cd8.clone.pre, do.center = T, do.scale = T)

cd8.clone.pre.scale <- cd8.clone.pre.scale[,cd8.clone.pre.scale$clone.id %in% pre.data$clone]
cd8.clone.pre.scale$type = ifelse(cd8.clone.pre.scale$clone.id %in% pre.data$clone[pre.data$group == "Tex"], "Tex", "nonTex")
Idents(cd8.clone.pre.scale) <- cd8.clone.pre.scale$type

tibble(cell = colnames(cd8.clone.pre.scale), type = cd8.clone.pre.scale$type) %>%
  dplyr::arrange(desc(type)) %>%
  dplyr::mutate(cellid = stringr::str_sub(cell, -18, -1))  -> tmp.rank

use.genes <- c("CXCL13","CTLA4","PDCD1","ENTPD1","LAG3","LAYN","TIGIT","BATF","HAVCR2","TNFRSF9","GZMB","CCL5","BAG3")
tmp.matr <- cd8.clone.pre.scale@assays$RNA@scale.data[use.genes,tmp.rank$cell]


colnames(tmp.matr) <- NULL

tmp.matr[tmp.matr < -2.5] <- -2.5
tmp.matr[tmp.matr > 3.5] <- 3.5

ha = HeatmapAnnotation(
  type = tmp.rank$type,
  col = list(
    cluster = c("Tex" = "#CD853F","nonTex" = "#009ACD")),
  simple_anno_size = unit(0.4, "cm")
)

ComplexHeatmap::Heatmap(tmp.matr, cluster_columns = F, cluster_rows = F, 
                        col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                        top_annotation = ha)
