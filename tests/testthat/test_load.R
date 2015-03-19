library(peprr)
context("testing load functions")

test_that("init_peprDB creates database",
          expect_equal(class(init_peprDB(data_dir = "data_test")),
                       list("src_sqlite", "src_sql","src"))
)

test_that("init_peprDB loads database",
          expect_equal(class(init_peprDB(db_path = "data_test")),
                       list("src_sqlite", "src_sql","src"))
)

## for testing ...
# peprDB <- init_peprDB(data_dir = "?? test data??",create = TRUE)
# qc_stats_dir <- "/Users/nolson/Desktop/micro_rm/pepr-data/MG001/Re.Order.RC.CFSAN008157.HGAP_miseq_qc_stats"
