# ===== Minimal libraries =====
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(DoubletFinder)
})

set.seed(12)

# ===== Read 6 samples & basic QC =====
zero_ABC  <- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132543/SRR33132543_out/outs/filtered_feature_bc_matrix/"), project="0h_ABC")
zero_DEF  <- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132542/SRR33132542_out/outs/filtered_feature_bc_matrix/"), project="0h_DEF")
four_ABC  <- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132541/SRR33132541_out/outs/filtered_feature_bc_matrix/"), project="4h_ABC")
four_DEF  <- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132540/SRR33132540_out/outs/filtered_feature_bc_matrix/"), project="4h_DEF")
twelve_ABC<- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132539/SRR33132539_out/outs/filtered_feature_bc_matrix/"), project="12h_ABC")
twelve_DEF<- CreateSeuratObject(Read10X("/home/qian@vitra.bio/ava/data/huang_data_mouse_whole_ovary/single_cell/SRR33132538/SRR33132538_out/outs/filtered_feature_bc_matrix/"), project="12h_DEF")

for (obj in list(zero_ABC, zero_DEF, four_ABC, four_DEF, twelve_ABC, twelve_DEF)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
}

# ===== Merge =====
alldata <- merge(zero_ABC, y = list(zero_DEF, four_ABC, four_DEF, twelve_ABC, twelve_DEF),
                 add.cell.ids = c("0h_ABC","0h_DEF","4h_ABC","4h_DEF","12h_ABC","12h_DEF"),
                 project = "alldata")

# mito & ribo %
alldata <- PercentageFeatureSet(alldata, "^mt-",   col.name = "percent_mito")

alldata <- PercentageFeatureSet(alldata, "^Rp[sl]", col.name = "percent_ribo")
VlnPlot(
  alldata,
  features = c("percent_mito", "percent_ribo",'nFeature_RNA'),
  pt.size = 0.001  # 控制散点大小，0 就不显示
)

dim(alldata)

# QC
alldata <- subset(alldata, subset = nFeature_RNA > 3000 & percent_mito < 20)
dim(alldata)

# ===== Standard workflow =====
alldata <- NormalizeData(alldata)
alldata <- FindVariableFeatures(alldata, selection.method = "vst")
alldata <- ScaleData(alldata, features = rownames(alldata))
alldata <- RunPCA(alldata, features = VariableFeatures(alldata))
alldata <- RunUMAP(alldata, dims = 1:20)
UMAPPlot(alldata)
alldata <- FindNeighbors(alldata, dims = 1:20)
alldata <- FindClusters(alldata, resolution = 0.1)
alldata <- JoinLayers(alldata)  

# ===== DoubletFinder =====
nExp <- round(ncol(alldata) * 0.05)
alldata <- doubletFinder(alldata, pN = 0.5, pK = 0.09, nExp = nExp, PCs = 1:20)
DF.name <- grep("DF.classification", colnames(alldata@meta.data), value = TRUE, fixed = FALSE)
alldata.filt <- alldata[, alldata@meta.data[, DF.name] == "Singlet"]
UMAPPlot(alldata.filt)
# ===== Reclustering after filtering =====
alldata.filt <- FindNeighbors(alldata.filt, graph.name = "test", dims = 1:20)
alldata.filt <- FindClusters(alldata.filt, resolution = 0.5, graph.name = "test")
DimPlot(alldata.filt,label = T,repel = T)

# ===== FindAllMarkers =====
alldata.filt <- JoinLayers(alldata.filt)
all_markers <- FindAllMarkers(alldata.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top20_per_cluster <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  ungroup()
DimPlot(alldata.filt,label = T,repel = T)



# ===== generate timepoints / timep =====
alldata.filt$timepoints <- dplyr::case_when(
  grepl("^0h",  alldata.filt$orig.ident) ~ "00hrs",
  grepl("^4h",  alldata.filt$orig.ident) ~ "04hrs",
  grepl("^12h", alldata.filt$orig.ident) ~ "12hrs",
  TRUE ~ NA_character_
)
alldata.filt$timep <- dplyr::case_when(
  grepl("^0h",  alldata.filt$orig.ident) ~ 0L,
  grepl("^4h",  alldata.filt$orig.ident) ~ 4L,
  grepl("^12h", alldata.filt$orig.ident) ~ 12L,
  TRUE ~ NA_integer_
)



DimPlot(alldata.filt,label = T,repel = T)
Idents(alldata.filt) <- alldata.filt$seurat_clusters
alldata.filt <- FindSubCluster(alldata.filt, "7", graph.name = "test", subcluster.name = "cl7", resolution = 0.1, algorithm = 1)
DimPlot(alldata.filt, reduction = "umap", group.by = c("cl7"), label = TRUE, label.size = 6)
unique(Idents(alldata.filt))
Idents(alldata.filt) <- "cl7"
# cluster70.markers <- FindMarkers(alldata.filt, ident.1 = "7_0", min.pct = 0.25)
# head(cluster70.markers[order(cluster70.markers$avg_log2FC,decreasing = TRUE),],20)
# cluster71.markers <- FindMarkers(alldata.filt, ident.1 = "7_1", min.pct = 0.25)
# head(cluster71.markers[order(cluster71.markers$avg_log2FC,decreasing = TRUE),],20)
Idents(alldata.filt) <- alldata.filt$cl7
levels(alldata.filt)
table(alldata.filt$cl7)
VlnPlot(alldata.filt,features = c('Ptprd','Grb14'))
VlnPlot(alldata.filt,features = c('Inhbb'))
VlnPlot(alldata.filt,features = c('Amh','Slc18a2'))
VlnPlot(alldata.filt,features = c('Sult1e1'))
VlnPlot(alldata.filt,features = c('Spp1'))
VlnPlot(alldata.filt,features = c('Cyp11a1','Gas6'))
VlnPlot(alldata.filt,features = c('Lama2','Plac8'))

new.cluster.ids <- c(
  "Meiotic Granulosa",             # 0
  "Stroma early",                    # 5
  "Large Antral Granulosa",          # 1
  "Myeloid/Macrophage",              # 11
  "small Luteal",                    # 8
  "Epithelial (OSE)",                # 14
  "Steroidogenic Theca early",        # 9
  "Endothelial",                     # 10
  "Stroma late",                      # 13
  "Preantral Granulosa",   # 7_0
  "Unknown",      # 7_1
  "Lymphatic Endothelial",           # 15
  "Steroidogenic Theca late",     # 6
  "Stroma late",                     # 4
  "COC Expansion Early",             # 2
  "COC Expansion Late",              # 12
  "large Luteal"                     # 3
)

names(new.cluster.ids) <- levels(alldata.filt)
alldata.filt <- RenameIdents(alldata.filt, new.cluster.ids)

DimPlot(alldata.filt, reduction = "umap", label = TRUE, label.size = 6)

table(Idents(alldata.filt))
alldata.filt <- subset(
  alldata.filt, 
  idents = setdiff(levels(alldata.filt), "Stroma-like (perivascular)")
)

DimPlot(alldata.filt, reduction = "umap", label = TRUE, label.size = 6)
levels(alldata.filt)


# 你当前 levels(alldata.filt) 列表里精确对应的映射（左 -> 右）
map_left_to_right <- c(
  "Meiotic Granulosa"          = "Granulosa 3",
  "Stroma early"               = "Stroma 1",
  "Large Antral Granulosa"     = "Granulosa 1",
  "Myeloid/Macrophage"         = "Myeloid",
  "small Luteal"               = "Luteal 1",
  "Epithelial (OSE)"           = "Epithelial",
  "Steroidogenic Theca early"  = "Theca 1",
  "Endothelial"                = "Endothelial 1",
  "Stroma late"                = "Stroma 2",
  "Preantral Granulosa"        = "Granulosa 2",
  "Lymphatic Endothelial"      = "Endothelial 2",
  "Steroidogenic Theca late"   = "Theca 2",
  "COC Expansion Early"        = "Cumulus 1",
  "COC Expansion Late"         = "Cumulus 2",
  "large Luteal"               = "Luteal 2"
)

# 诊断：哪些 level 没在映射里（会保持原名）
missing_in_map <- setdiff(levels(alldata.filt), names(map_left_to_right))
if (length(missing_in_map) > 0) {
  message("这些 levels 当前没有映射（将保留原名）: ",
          paste(missing_in_map, collapse = ", "))
}

# 进行稳妥的显式替换
old_ids <- Idents(alldata.filt)
new_ids <- as.character(old_ids)
for (k in names(map_left_to_right)) {
  new_ids[old_ids == k] <- map_left_to_right[[k]]
}
Idents(alldata.filt) <- factor(new_ids)
alldata.filt$celltype_new <- Idents(alldata.filt)
table(Idents(alldata.filt))
Idents(alldata.filt) <- factor(
  Idents(alldata.filt),
  levels = c("Cumulus 1","Cumulus 2","Endothelial 1","Endothelial 2","Epithelial",
             "Granulosa 1","Granulosa 2","Granulosa 3","Luteal 1","Luteal 2",
             "Myeloid","Oocyte","Stroma 1","Stroma 2","Theca 1","Theca 2",
             "Unknown 1","Unknown 2")
)

DimPlot(alldata.filt, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(alldata.filt, reduction = "umap", label = TRUE, pt.size = 0.1,repel=T,
        cols = c("#8dd96d",#cumu1
                 "#3f6b00",#cumu2
                 "#e00089",#endo1
                 "#970f84",#endo2
                 "#9e2409",#epi
                 "#f5be37",#gc1
                 "#ff7d3e",#gc2
                 "#ffa972",#gc3
                 "#f2b1ed",#luteal1
                 "#ff91d4",#luteal2
                 "#b151d8",#my
                 "#777200",#ooc
                 "#006ec9",#stroma1
                 "#355286",#stroma2
                 "#7dda8a",#theca1
                 "#1e5f28"#theca2
                 ))

# savedata
setwd("~/ava/data/huang_data_mouse_whole_ovary/single_cell")
saveRDS(alldata.filt, "10x_filt_ava_annotated.rds")


##### cellchat
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(DoubletFinder)
})

setwd("~/ava/data/huang_data_mouse_whole_ovary/single_cell")
alldata.filt <- readRDS("~/ava/data/huang_data_mouse_whole_ovary/single_cell/10x_filt.rds")
DimPlot(alldata.filt,label = T,repel = T)
table(Idents(alldata.filt))
table(Idents(alldata.filt),alldata.filt$orig.ident)
library(pheatmap)
tab <- table(Idents(alldata.filt), alldata.filt$orig.ident)
pheatmap(log1p(tab), 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         scale = "row", 
         main = "Cell type abundance heatmap (log1p)")



alldata.filt$celltype <- Idents(alldata.filt)
seurat_obj <- subset(alldata.filt, subset = orig.ident == "0h_DEF")
seurat_obj
DimPlot(seurat_obj,label = T,repel = T)

# create input for aws tangram input 
# X <- as.matrix(GetAssayData(seurat_obj, assay = "RNA", slot = "data"))
# write.csv(X, "X_matrix.csv")
# write.csv(rownames(X), "gene_names.csv", row.names = FALSE)
# write.csv(data.frame(cell_type = Idents(seurat_obj)), "cell_type.csv")



data.counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
labels <- seurat_obj$celltype
labels <- labels[colnames(seurat_obj)]
stopifnot(identical(names(labels), colnames(seurat_obj)))

labels <- factor(labels)
labels <- droplevels(labels)
if (any(is.na(labels))) {
  keep <- !is.na(labels)
  labels <- labels[keep]
  data.counts <- data.counts[, keep, drop = FALSE]
}

library(CellChat)
meta <- data.frame(group = labels, row.names = names(labels), stringsAsFactors = FALSE)
stopifnot(identical(rownames(meta), colnames(data.counts)))
cellchat <- createCellChat(object = data.counts, meta = meta, group.by = "group")
cellchat <- setIdent(cellchat, ident.use = "group")
cellchat@idents <- droplevels(factor(cellchat@idents))
cellchatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(cellchatDB, search = c("Secreted Signaling"))
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)                    
cellchat <- identifyOverExpressedGenes(cellchat)     
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 5)
df.net1 <- subsetCommunication(cellchat)
unique(df.net1$target)
df.net <- df.net1[df.net1$target %in% c('Large Antral Granulosa','Preantral Granulosa'),]
df.net <- df.net[!df.net$source %in% 'Unknown (neuronal/epithelial)',]
l <- unique(df.net$ligand)

# read in the spatial deg computed from aws merfish data analysis
spatial_deg <- read.csv('moranI_top2500_genes.csv')
spatial_deg$X[1:50]
spatial_deg$X <- tools::toTitleCase(spatial_deg$X)
l %in% spatial_deg$X
seurat_obj_filtered <- subset(seurat_obj, subset = celltype != "Unknown (neuronal/epithelial)")
DimPlot(seurat_obj_filtered,label = T,repel = T)


# read in cellchat results from HTS (yuichiro data)
cat <- c(
  "Bmp2", "Bmp4", "Bmp6", "Bmp7", "Gdf11", "Amh", "Inhibin B", "Wnt10a", 
  "Wnt10b", "Wnt2", "Wnt4", "Wnt6", "Wnt9a", "Wnt5a", "Tgfa", "Areg", 
  "Btc", "Hbegf", "Ereg", "Fgf1", "Fgf2", "Fgf18", "Fgf9", "Fgf16", 
  "Fgf21", "Igf1", "Dhh", "Ihh", "Cxcl12", "Ctf1", "Lif", "Osm", 
  "Tnf", "Tnfsf12", "Spp1", "Retn", "Nampt", "Angptl2", "Angptl4", "Angpt2", 
  "Mdk", "Edn1", "Edn2", "Nppc", "Gas6", "Pros1"
)
cat %in% spatial_deg$X

cat[cat %in% spatial_deg$X]

cat[cat %in% l]
