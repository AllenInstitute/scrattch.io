library(dplyr)
library(purrr)
library(Matrix)
devtools::load_all()
#library(scrattch.io)
options(stringsAsFactors = FALSE)

# Get dataset from tasic2016data
library(tasic2016data)

select_cells <- tasic_2016_anno %>%
  filter(primary_type != "unclassified") %>%
  filter(grepl("Sst",primary_type)) %>%
  select(sample_name) %>%
  unlist()

# raw_anno will be used by formats like loom
raw_anno <- tasic_2016_anno %>%
  filter(sample_name %in% select_cells)

# anno will be used for tome and feather
anno <- raw_anno %>%
  annotate_cat("mouse_line") %>%
  annotate_cat("cre_driver_1") %>%
  annotate_cat("cre_driver_2") %>%
  annotate_cat("cre_reporter") %>%
  annotate_cat("dissection") %>%
  annotate_cat("tdTomato") %>%
  annotate_cat("pass_qc_checks") %>%
  annotate_cat("broad_type") %>%
  annotate_cat("core_intermediate") %>%
  rename("primary_type_label" = "primary_type") %>%
  rename("secondary_type_label" = "secondary_type")

# everything will use the data matrix
data <- tasic_2016_counts[,select_cells]
datavar <- apply(data, 1, var)
keepdata <- names(datavar[order(-datavar)][1:2000])
data <- data[keepdata,]

sparse_data <- as(data, "dgCMatrix")

# desc will be used by tome and feather
desc <- data.frame(base = c("primary_type","secondary_type","core_intermediate","broad_type",
                            "mouse_line","cre_driver_1","cre_driver_2","cre_reporter",
                            "dissection","tdTomato","pass_qc_checks"),
                   name = c("Primary Cell Type","Secondary Cell Type","Core/Intermediate","Broad Cell Type",
                            "Donor Line","Donor Driver 1","Donor Driver 2","Donor Reporter Line",
                            "Dissected Layers","tdTomato expression","Passes QC Checks"),
                   type = rep("cat",11))

# Make cluster medians
cluster_ids <- unique(anno$primary_type_id)
medians <- map_dfc(cluster_ids,
               function(x) {
                 samples <- anno$sample_name[anno$primary_type_id == x]
                 apply(data[,samples], 1, median)
               })
names(medians) <- paste0("primary_type_",cluster_ids)

# hclust the clusters

cluster_dist <- dist(t(as.matrix(medians)))
cluster_hc <- hclust(cluster_dist)
cluster_hc$height <- cluster_hc$height/max(cluster_hc$height)

# make a dendrogram
library(dendextend)

cluster_anno <- anno %>%
  select(primary_type_id, primary_type_label, primary_type_color) %>%
  unique()

dend <- as.dendrogram(cluster_hc)
dend_cluster_ids <- as.numeric(sub("primary_type_","",labels(dend)))
labels(dend) <- cluster_anno$primary_type_label[match(dend_cluster_ids, cluster_anno$primary_type_id)]
labels_colors(dend) <- cluster_anno$primary_type_color[match(dend_cluster_ids, cluster_anno$primary_type_id)]

dend_desc <- data.frame(base = "primary_type",
                        name = "Primary Type Hierarchy",
                        links_to = "primary_type_label")

# Save as rda so that we know what to expect
saveRDS(data, file = "rds/data.RData")
saveRDS(anno, file = "rds/anno.RData")
saveRDS(desc, file = "rds/desc.RData")
saveRDS(sparse_data, file = "rds/data_dgCMatrix.RData")
saveRDS(dend, file = "rds/dend.RData")

## Write feather files
library(feather)

fdata <- cbind(sample_name = colnames(data), as.data.frame(t(data)))

write_feather(fdata, "feather/data.feather")
write_feather(anno, "feather/anno.feather")
write_feather(desc, "feather/desc.feather")

# Write tome file
write_tome_data(exon_mat = as(data,"dgCMatrix"),
                tome = "tome/transcrip.tome",
                cols_are = "sample",
                overwrite = TRUE)
write_tome_anno(anno, "tome/transcrip.tome")
write_tome_anno_desc(desc, "tome/transcrip.tome")
write_tome_dend(dend, "primary_type" ,"tome/transcrip.tome")
write_tome_dend_desc(dend_desc, "tome/transcrip.tome")

# Fetch a small Loom file from loom.linnarssonlab.org

download.file("http://loom.linnarssonlab.org/dataset/cellmetadata/Mousebrain.org.level1/L1_Olfactory.agg.loom",
              destfile = "loom/L1_Olfactory.agg.loom")
