library(scrattch.io)

context("annotation functions")

test_that(
  "cl_to_anno() converts from hicat's cl object to a data.frame",
  {
    cl <- readRDS(system.file("testdata/hicat", "tasic16_cl.RData", package = "scrattch.io"))

    result <- cl_to_anno(cl,
                         cl.df = NULL,
                         cl_col = "cl",
                         base = "cluster")

    expect_is(result, "data.frame")
    expect_equal(nrow(result), length(cl))
    expect_equal(names(result), c("sample_name","cluster_id","cluster_label","cluster_color"))
    expect_equal(sum(unique(result$cluster_label) %in% levels(cl)), length(levels(cl)))

  }
)

test_that(
  "cl_to_anno() uses cl.df to annotate hicat's cl object as a data.frame",
  {
    cl <- readRDS(system.file("testdata/hicat", "tasic16_cl.RData", package = "scrattch.io"))
    cl.df <- readRDS(system.file("testdata/hicat", "tasic16_cl.df.RData", package = "scrattch.io"))

    result <- cl_to_anno(cl,
                         cl.df = cl.df,
                         cl_col = "cl",
                         base = "cluster")

    expect_is(result, "data.frame")
    expect_equal(nrow(result), length(cl))
    expect_equal(names(result), c("sample_name","cluster_id","cluster_label","cluster_color"))

    result2 <- cl_to_anno(cl,
                          cl.df = cl.df[,names(cl.df) != "cluster_label"],
                          cl_col = "cl",
                          base = "cluster")

    expect_is(result2, "data.frame")
    expect_equal(nrow(result2), length(cl))
    expect_equal(names(result2), c("sample_name","cluster_id","cluster_label","cluster_color"))

    result3 <- cl_to_anno(cl,
                          cl.df = cl.df[,!names(cl.df) %in% c("cluster_label","cluster_color")],
                          cl_col = "cl",
                          base = "cluster")

    expect_is(result3, "data.frame")
    expect_equal(nrow(result3), length(cl))
    expect_equal(names(result3), c("sample_name","cluster_id","cluster_label","cluster_color"))


    result4 <- cl_to_anno(cl,
                          cl.df = cl.df[,!names(cl.df) %in% c("cluster_id","cluster_color")],
                          cl_col = "cl",
                          base = "cluster")

    expect_is(result4, "data.frame")
    expect_equal(nrow(result4), length(cl))
    expect_equal(names(result4), c("sample_name","cluster_id","cluster_label","cluster_color"))

    cl.df2 <- cl.df
    rownames(cl.df2) <- NULL

    expect_error(cl_to_anno(cl,
                            cl.df = cl.df2,
                            cl_col = "wrong",
                            base = "cluster"))


    cl.df3 <- cl.df
    rownames(cl.df3) <- letters[1:nrow(cl.df)]

    expect_error(cl_to_anno(cl,
                            cl.df = cl.df3,
                            cl_col = "wrong",
                            base = "cluster"))

    }
)
