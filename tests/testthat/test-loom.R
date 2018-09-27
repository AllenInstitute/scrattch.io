library(scrattch.io)

# Test data are available in inst/testdata/
# A loom file for testing is at inst/testdata/loom/Sympathetic.loom

context("Read Loom files")

loom_file <- system.file("testdata/loom",
                         "Sympathetic.loom",
                         package = "scrattch.io")

test_that(
  "read_loom_dgCMatrix can read a loom matrix.",
  {
    loom_mat <- read_loom_dgCMatrix(loom_file,
                                    row_names = "Gene",
                                    col_names = "Cell_id")
    expect_is(loom_mat, "dgCMatrix")
  }
)
