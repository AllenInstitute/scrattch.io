library(dplyr)
library(scrattch.io)
options(stringsAsFactors = FALSE)

# Get dataset from tasic2016data
library(tasic2016data)

select_cells <- tasic_2016_anno %>%
  filter(primary_type != "unclassified") %>%
  filter(grepl("Sst",primary_type)) %>%
  select(sample_id) %>%
  unlist()

# raw_anno will be used by formats like loom
raw_anno <- tasic_2016_anno %>%
  filter(sample_id %in% select_cells)

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

# desc will be used by tome and feather
desc <- data.frame(base = c("primary_type","secondary_type","core_intermediate","broad_type",
                            "mouse_line","cre_driver_1","cre_driver_2","cre_reporter",
                            "dissection","tdTomato","pass_qc_checks"),
                   name = c("Primary Cell Type","Secondary Cell Type","Core/Intermediate","Broad Cell Type",
                            "Donor Line","Donor Driver 1","Donor Driver 2","Donor Reporter Line",
                            "Dissected Layers","tdTomato expression","Passes QC Checks"),
                   type = rep("cat",11))

## Write feather files
library(feather)

fdata <- cbind(sample_id = colnames(data), as.data.frame(t(data)))

write_feather(fdata, "feather/data.feather")
write_feather(anno, "feather/anno.feather")
write_feather(desc, "feather/desc.feather")

# Write tome file
library(Matrix)

write_tome_data(exon_mat = as(data,"dgCMatrix"),
                tome = "tome/transcrip.tome",
                cols_are = "sample",
                overwrite = TRUE)
write_tome_anno(anno, "tome/transcrip.tome")
write_tome_anno_desc(desc, "tome/transcrip.tome")

