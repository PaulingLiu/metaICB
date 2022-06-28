

matr <- readr::read_tsv("projects/07_pan.cancer.treatment/02.data/GSE156728_RC_10X.CD8.counts.txt.gz")

gene <- matr$X1
matr <- matr[,-1]
matr <- as.matrix(matr)
rownames(matr) <- gene

mda <- tibble(cells = colnames(matr)) %>%
  tidyr::separate(cells, c("cellID","patient"), sep = "\\.", remove = F)

mda <- as.data.frame(mda)
rownames(mda) <- mda$cells

rc.sce <- CreateSeuratObject(counts = matr, meta.data = mda)
rc.sce <- NormalizeData(rc.sce)
rc.sce <- ScaleData(rc.sce, do.center = F, do.scale = F)

ggplot(aes(rc.sce@assays$RNA@scale.data["CXCL13",]), data = NULL) +
  geom_density() +
  geom_vline(xintercept = 0.3)

tex.sce <- rc.sce[,rc.sce@assays$RNA@scale.data["CXCL13",] > 0.3]
tex.sce <- seurat.process(tex.sce)
tex.sce <- bbknn.batch(tex.sce)

FeaturePlot(tex.sce, features = "GZMK")
