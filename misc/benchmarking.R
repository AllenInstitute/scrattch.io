library(microbenchmark)
devtools::load_all()

tome_file <- tome_file <- system.file("testdata/tome",
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
  times = 5
)
# About 10X faster to read the whole matrix in and then subset if the matrix is small.
