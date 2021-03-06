context("Read .tome data.frame")
library(scrattch.io)

# Test data are available in inst/testdata/
# A tome file for testing is at inst/testdata/tome/transcrip.tome
# Reference files as .RData are in inst/testdata/rds/

tome_file <- system.file("testdata/tome",
                         "transcrip.tome",
                         package = "scrattch.io")

test_that(
  "read_tome_data.frame reads a compound object as a data.frame",
  {
    tome_desc <- read_tome_data.frame(tome = tome_file,
                                      df_name = "/sample_meta/desc",
                                      stored_as = "data.frame")

    expect_is(tome_desc, "data.frame")

    expect_equal_to_reference(
      tome_desc,
      system.file("testdata",
                  "rds/desc.RData",
                  package = "scrattch.io")
    )
  }
)

test_anno <- readRDS(system.file("testdata",
                                 "rds/anno.RData",
                                 package = "scrattch.io"))

test_that(
  "read_tome_data.frame can read vectorized objects as a data.frame",
  {

    tome_anno <- read_tome_data.frame(tome = tome_file,
                                      df_name = "/sample_meta/anno",
                                      stored_as = "vectors")
    tome_anno <- tome_anno[,names(test_anno)]

    expect_is(tome_anno, "data.frame")

    expect_equal_to_reference(
      tome_anno,
      system.file("testdata",
                  "rds/anno.RData",
                  package = "scrattch.io")
    )
  }
)

test_that(
  "read_tome_data.frame can read selected columns from a vectorized object as a data.frame",
  {
    tome_anno_select <- read_tome_data.frame(tome = tome_file,
                                             df_name = "/sample_meta/anno",
                                             columns = names(test_anno)[1:5],
                                             stored_as = "vectors")

    expect_is(tome_anno_select, "data.frame")

    expect_equal(tome_anno_select,
                 test_anno[,1:5])

    tome_anno_select_grep <- read_tome_data.frame(tome = tome_file,
                                                  df_name = "/sample_meta/anno",
                                                  columns = "primary_type",
                                                  match_type = "grep",
                                                  stored_as = "vectors")

    expect_equal(tome_anno_select_grep,
                 test_anno[,c("primary_type_color","primary_type_id","primary_type_label")])
  }
)

test_that(
  "read_tome_data.frame returns an error when no matching columns are requested",
  {

    expect_error(
      read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = "Notacolumn",
                           match_type = "exact",
                           stored_as = "vectors")
    )

    expect_error(
      read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = "Notacolumn",
                           match_type = "grep",
                           stored_as = "vectors"
                           )
    )

  }
)

test_that(
  "read_tome_data.frame returns a warning when matching and non-matching columns are requested",
  {

    expect_warning(
      read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = c("primary_type_label","Notacolumn"),
                           stored_as = "vectors",
                           match_type = "exact")
    )

    expect_warning(
      read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = c("primary_type","Notacolumn"),
                           stored_as = "vectors",
                           match_type = "grep")
    )

  }
)

test_that(
  "read_tome_data.frame arranges columns as expected when using get_all = TRUE",
  {

    primary_cols <- names(test_anno)[grepl("primary_type", names(test_anno))]

    tome_anno_expect <- test_anno[,c(primary_cols, setdiff(names(test_anno), primary_cols))]

    tome_anno_exact <- read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = primary_cols,
                           stored_as = "vectors",
                           match_type = "exact",
                           get_all = TRUE)

    expect_equal(tome_anno_exact, tome_anno_expect)

    tome_anno_grep <- read_tome_data.frame(tome = tome_file,
                           df_name = "/sample_meta/anno",
                           columns = "primary",
                           stored_as = "vectors",
                           match_type = "grep",
                           get_all = TRUE)

    expect_equal(tome_anno_grep, tome_anno_expect)


  }
)

context("Read .tome dgCMatrix")

test_that(
  "read_tome_dgCMatrix retrieves the entire matrix as a sparse matrix object",
  {
    tome_dgC <- read_tome_dgCMatrix(tome = tome_file,
                                    "/data/exon/")

    expect_equal_to_reference(tome_dgC,
                              system.file("testdata/rds",
                                          "data_dgCMatrix.RData",
                                          package = "scrattch.io"))

    tome_dgC_t <- read_tome_dgCMatrix(tome = tome_file,
                                      "/data/t_exon/")

    tome_dgC_t <- Matrix::t(tome_dgC_t)

    expect_equal_to_reference(tome_dgC_t,
                              system.file("testdata/rds",
                                          "data_dgCMatrix.RData",
                                          package = "scrattch.io"))
  }
)

context("Read .tome serialized")

test_that(
  "read_tome_serialized retrieves and unserializes an object.",
  {

    tome_serial <- read_tome_serialized(tome = tome_file,
                                        "/dend/primary_type")

    expect_equal_to_reference(tome_serial,
                              system.file("testdata/rds",
                                          "dend.RData",
                                          package = "scrattch.io"))

  }

)

context("Read .tome gene and sample data")

test_that(
  "read_tome_gene_data retrieves expression values for selected genes.",
  {
    tome_one_gene <- read_tome_gene_data(tome = tome_file,
                                         genes = "Sst",
                                         region = "exon",
                                         units = "counts",
                                         transform = "none",
                                         format = "data.frame")

    expect_is(tome_one_gene, "data.frame")
    expect_equal(ncol(tome_one_gene), 2)

    # Missing genes give an error.
    expect_error(read_tome_gene_data(tome = tome_file,
                                     genes = "Pvalb"))

    tome_multi_gene <- read_tome_gene_data(tome = tome_file,
                                           genes = c("Npy","Sst","Snap25"),
                                           region = "exon",
                                           units = "counts",
                                           transform = "none",
                                           format = "data.frame")

    expect_is(tome_multi_gene, "data.frame")
    expect_equal(ncol(tome_multi_gene), 4)

    tome_multi_gene_matrix <- read_tome_gene_data(tome = tome_file,
                                           genes = c("Npy","Sst","Snap25"),
                                           region = "exon",
                                           units = "counts",
                                           transform = "none",
                                           format = "matrix")

    expect_is(tome_multi_gene_matrix, "matrix")
    expect_equal(ncol(tome_multi_gene_matrix), 3)

    tome_multi_gene_dgCMatrix <- read_tome_gene_data(tome = tome_file,
                                                  genes = c("Npy","Sst","Snap25"),
                                                  region = "exon",
                                                  units = "counts",
                                                  transform = "none",
                                                  format = "dgCMatrix")

    expect_is(tome_multi_gene_dgCMatrix, "dgCMatrix")
    expect_equal(ncol(tome_multi_gene_dgCMatrix), 3)

  }
)

test_that(
  "read_tome_sample_data retrieves expression values for selected samples.",
  {
    tome_one_sample <- read_tome_sample_data(tome = tome_file,
                                         samples = "Gad2_tdTpositive_cell_68",
                                         region = "exon",
                                         units = "counts",
                                         transform = "none",
                                         format = "data.frame")

    expect_is(tome_one_sample, "data.frame")
    expect_equal(ncol(tome_one_sample), 2)

    # Missing samples give an error.
    expect_error(read_tome_sample_data(tome = tome_file,
                                     samples = "Pvalb"))

    samples_3 <- c("Gad2_tdTpositive_cell_68","Htr3a_tdTpositive_cell_84","Gad2_tdTpositive_cell_17")

    tome_multi_sample <- read_tome_sample_data(tome = tome_file,
                                           samples = samples_3,
                                           region = "exon",
                                           units = "counts",
                                           transform = "none",
                                           format = "data.frame")

    expect_is(tome_multi_sample, "data.frame")
    expect_equal(ncol(tome_multi_sample), 4)

    tome_multi_sample_matrix <- read_tome_sample_data(tome = tome_file,
                                                  samples = samples_3,
                                                  region = "exon",
                                                  units = "counts",
                                                  transform = "none",
                                                  format = "matrix")

    expect_is(tome_multi_sample_matrix, "matrix")
    expect_equal(ncol(tome_multi_sample_matrix), 3)

    tome_multi_sample_dgCMatrix <- read_tome_sample_data(tome = tome_file,
                                                     samples = samples_3,
                                                     region = "exon",
                                                     units = "counts",
                                                     transform = "none",
                                                     format = "dgCMatrix")

    expect_is(tome_multi_sample_dgCMatrix, "dgCMatrix")
    expect_equal(ncol(tome_multi_sample_dgCMatrix), 3)

  }
)
