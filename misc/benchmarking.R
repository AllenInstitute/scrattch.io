library(microbenchmark)
devtools::load_all()

tome_file <- system.file("testdata/tome",
                         "transcrip.tome",
                         package = "scrattch.io")

microbenchmark(
  read_tome_gene_data(tome_file, regions = "both", format = "dgCMatrix"),
  read_tome_dgCMatrix(tome_file, "/data/exon") + read_tome_dgCMatrix(tome_file, "/data/intron"),
  times = 5
)

genes <- read_tome_gene_names(tome_file)[sample(1:2000, 50)]
microbenchmark(
  read_tome_gene_data(tome_file, genes = genes, regions = "both", format = "dgCMatrix"),
  {
    g <- read_tome_gene_names(tome_file)
    g_col <- match(genes, g)
    ex <- read_tome_dgCMatrix(tome_file, "/data/exon")[,g_col]
    int <- read_tome_dgCMatrix(tome_file, "/data/intron")[,g_col]
    ex + int
    },
  times = 100
)
# About 10X faster to read the whole matrix in and then subset if the matrix is small.


index_list <- map2(exon_starts, exon_ends, function(x,y) {x:y})

test_read_1 <- function(tome, genes, regions = "exon") {
  gene_names <- read_tome_vector(tome,"/gene_names")

  if(is.null(genes)) {
    genes <- gene_names
  } else {

    # Check for genes not found in the dataset.
    missing_genes <- setdiff(genes, gene_names)

    if(length(missing_genes) > 0) {
      stop("Some genes not found in tome: ", paste(missing_genes, collapse = ", "))
    }

  }

  sample_names <- read_tome_vector(tome,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- read_tome_vector(tome, "data/exon/p")[gene_index] + 1
    exon_ends <- read_tome_vector(tome, "data/exon/p")[(gene_index + 1)]

    exon_values <- purrr::map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- rhdf5::h5read(tome, "data/exon/x", index = list(start:end))
      } else {
        values <- NA
      }
      as.vector(values)
    })

    exon_sample_indexes <- purrr::map2(exon_starts, exon_ends, function(start, end) {
      if(end > start) {
        values <- rhdf5::h5read(tome, "data/exon/i", index = list(start:end))
      } else {
        values <- NA
      }
      as.vector(values)
    })

    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = c(0, cumsum(map_int(exon_values, length))))


  }
  h5closeAll()
  exon
}

a <- test_read_1(tome_file, genes)

test_read_2 <- function(tome, genes, regions = "exon") {
  gene_names <- read_tome_vector(tome,"/gene_names")

  if(is.null(genes)) {
    genes <- gene_names
  } else {

    # Check for genes not found in the dataset.
    missing_genes <- setdiff(genes, gene_names)

    if(length(missing_genes) > 0) {
      stop("Some genes not found in tome: ", paste(missing_genes, collapse = ", "))
    }

  }

  sample_names <- read_tome_vector(tome,"/sample_names")

  gene_index <- match(genes, gene_names)

  ## Exon values
  if(regions == "exon" | regions == "both") {

    exon_starts <- read_tome_vector(tome, "data/exon/p")[gene_index] + 1
    exon_ends <- read_tome_vector(tome, "data/exon/p")[(gene_index + 1)]

    index_list <- unlist(map2(exon_starts, exon_ends, function(x,y) {x:y}))

    split_list <- unlist(map(1:length(exon_starts), function(x) { rep(x, length(exon_starts[x]:exon_ends[x])) }))

    exon_values <- rhdf5::h5read(tome, "data/exon/x", index = list(index_list))
    exon_values <- split(exon_values, split_list)
    names(exon_values) <- NULL

    exon_sample_indexes <- rhdf5::h5read(tome, "data/exon/i", index = list(index_list))
    exon_sample_indexes <- split(exon_sample_indexes, split_list)
    names(exon_sample_indexes) <- NULL
    # exon_values <- purrr::map2(exon_starts, exon_ends, function(start, end) {
    #   if(end > start) {
    #     values <- rhdf5::h5read(tome, "data/exon/x", index = list(start:end))
    #   } else {
    #     values <- NA
    #   }
    #   values
    # })

    # exon_sample_indexes <- purrr::map2(exon_starts, exon_ends, function(start, end) {
    #   if(end > start) {
    #     values <- rhdf5::h5read(tome, "data/exon/i", index = list(start:end))
    #   } else {
    #     values <- NA
    #   }
    #   values
    # })

    exon <- list(x = exon_values,
                 i = exon_sample_indexes,
                 p = c(0, cumsum(map_int(exon_values, length))))
    names(exon$p) <- NULL

  }
  h5closeAll()
  exon
}

b <- test_read_2(tome_file, genes)


genes <- read_tome_gene_names(tome_file)[sample(1:2000, 50)]
microbenchmark(
  #read_tome_gene_data(tome_file, genes = genes, regions = "both", format = "dgCMatrix"),
  test_read_1(tome_file, genes),
  test_read_2(tome_file, genes),
  times = 5
)


old_read_tome_genes_jagged <-
