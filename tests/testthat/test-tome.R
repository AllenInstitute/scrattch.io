context("Read Tome Files")
library(scrattch.io)

test_that(
  "read_tome_data.frame reads a compound object as a data.frame",

  expect_equal_to_reference(
    read_tome_data.frame(tome = system.file("testdata",
                                            "tome/transcrip.tome",
                                            package = "scrattch.io"),
                                  df_name = "/sample_meta/desc",
                                  stored_as = "data.frame"),
    system.file("testdata",
                "rds/desc.RData",
                package = "scrattch.io")
  )
)
