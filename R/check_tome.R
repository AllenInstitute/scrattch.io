check_tome_expected_objects <- data.frame(
  full_name = c(
    "/data",
    "/data/exon",
    "/data/exon/dims",
    "/data/exon/i",
    "/data/exon/p",
    "/data/exon/x",
    "/data/exon_t",
    "/data/exon_t/dims",
    "/data/exon_t/i",
    "/data/exon_t/p",
    "/data/exon_t/x",
    "/data/intron",
    "/data/intron/dims",
    "/data/intron/i",
    "/data/intron/p",
    "/data/intron/x",
    "/data/intron_t",
    "/data/intron_t/dims",
    "/data/intron_t/i",
    "/data/intron_t/p",
    "/data/intron_t/x",
    "/gene_names",
    "/sample_names",
    "/dend",
    "/dend/desc",
    "/sample_meta/anno",
    "/sample_meta/desc"
  ),
  cat_name = c(
    "data",
    "exon data",
    "exon dims",
    "exon indices",
    "exon pointers",
    "exon values",
    "transposed exon data",
    "transposed exon dims",
    "transposed exon indices",
    "transposed exon pointers",
    "transposed exon values",
    "intron data",
    "intron dims",
    "intron indices",
    "intron pointers",
    "intron values",
    "transposed intron data",
    "transposed intron dims",
    "transposed intron indices",
    "transposed intron pointers",
    "transposed intron values",
    "gene names",
    "sample names",
    "dendrograms",
    "dendrogram descriptions",
    "sample annotations",
    "sample annotation descriptions"
  )
)

check_tome <- function(tome_file) {

  tome_ls <- rhdf5::h5ls(tome_file) %>%
    mutate(full_name = paste0(group,"/",name)) %>%
    mutate(full_name = sub("//","/",full_name))

  obj_names <- tome_ls$full_name

  name_checks <- check_tome_expected_objects %>%
    dplyr::mutate(present = full_name %in% obj_names)

  name_checks

}
