## 跨语言对拍：R 侧结果必须与 Python 引擎逐位相同。
##
## fixture 由 data-raw/make_fixtures.py 生成：GSE63703 的 47 个谱系 marker
## 加 3000 个随机背景基因，两边从同一份矩阵跑到最终分数。期望值来自
## 引擎本身（sample_types.layer_consistency + headline_score），不是手填。

fix_dir <- function(...) test_path("fixtures", ...)

read_fixture_matrix <- function(name) {
  d <- utils::read.delim(fix_dir(name), row.names = 1, check.names = FALSE)
  as.matrix(d)
}

test_that("purity reproduces the engine bit for bit", {
  skip_if_not(file.exists(fix_dir("purity_expected.csv")),
              "fixtures not generated; run data-raw/make_fixtures.py")
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  exp <- utils::read.csv(fix_dir("purity_expected.csv"), stringsAsFactors = FALSE)

  got <- ery_purity(m, layer = "transcriptome", level = "high")

  expect_equal(nrow(got), nrow(exp))
  expect_equal(got$sample, exp$sample)

  ## 主读数与一致性分：这是最终交付给人的数字，容差收到 0.001
  ## （两边都已 round 到一位小数，0.001 只吸收浮点末位）。
  expect_equal(got$headline, exp$headline, tolerance = 1e-6)
  expect_equal(got$ery_consistency, exp$ery_consistency, tolerance = 1e-6)
  expect_equal(got$consistency_rho, exp$consistency_rho, tolerance = 1e-6)
  expect_equal(got$prior, exp$prior, tolerance = 1e-6)
  expect_equal(got$ceiling, exp$ceiling, tolerance = 1e-6)

  ## 中间量：秩、AUC、logistic 分。容差 0.002 = 两个末位，
  ## 吸收 round() 在 .5 边界上的浮点抖动。
  expect_equal(got$ery_rank, exp$ery_rank, tolerance = 2e-3)
  expect_equal(got$leu_rank, exp$leu_rank, tolerance = 2e-3)
  expect_equal(got$delta, exp$delta, tolerance = 2e-3)
  expect_equal(got$auc, exp$auc, tolerance = 2e-3)
  expect_equal(got$rna_identity_score, exp$rna_identity_score, tolerance = 2e-3)
  expect_equal(got$n_markers, exp$n_markers)
  expect_equal(got$conflict, as.logical(exp$conflict))
})

test_that("the primary reading separates pure RBC from PBMC", {
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_purity(m, layer = "transcriptome", level = "high")
  rbc <- got$headline[grepl("^RBC", got$sample)]
  pbmc <- got$headline[grepl("^PBMC", got$sample)]
  expect_true(all(rbc > 85))
  expect_true(all(pbmc < 45))
  ## 主读数封顶在 p95 换算值，永不到 100
  expect_true(all(got$headline <= got$ceiling))
  expect_true(all(got$ceiling < 100))
})

test_that("no observation means abstention, never the prior alone", {
  ## unverifiable 权重为 0：不把库内先验当成本样本读数报出去。
  ## 这是最容易滑向「动不动 90 分」的路径，单独测一次。
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_purity(m, layer = "transcriptome", level = "unverifiable")
  expect_true(all(is.na(got$headline)))
  ## 但先验与天花板照常给出，便于解释为什么弃权
  expect_true(all(!is.na(got$prior)))
  expect_true(all(got$prior > 0))
})

test_that("a matrix without lineage markers fails loudly", {
  m <- matrix(abs(stats::rnorm(200)), nrow = 20,
              dimnames = list(paste0("GENE", 1:20), paste0("S", 1:10)))
  expect_error(ery_purity(m), "erythroid markers")
})

test_that("compute primitives match the engine definitions", {
  ## 这些是纯粹数学，一旦与引擎分叉，所有模块都会静默偏移。
  v <- c(3, 1, 4, 1, 5, NA, 9, 2, 6)
  ## 0-based 平均秩，NA 保持
  r <- ErythroWorkbench:::.rank_profile(v)
  expect_equal(r[which(!is.na(v))], sort(rank(v, ties.method = "average") - 1)[order(order(v[!is.na(v)]))])
  expect_true(is.na(r[6]))

  ## ddof = 0：R 的 sd() 是 ddof = 1，直接换会差 sqrt(n/(n-1))
  x <- c(1, 2, 3, 4, NA)
  expect_equal(ErythroWorkbench:::.nansd0(x), sqrt(mean((c(1, 2, 3, 4) - 2.5)^2)))
  expect_false(isTRUE(all.equal(ErythroWorkbench:::.nansd0(x), stats::sd(x, na.rm = TRUE))))

  ## 尺度判定：p99 > 50 才取 log2
  lin <- matrix(exp(seq(0, 12, length.out = 100)), nrow = 10)
  expect_equal(ery_as_log2(lin)$scale, "linear")
  expect_equal(ery_as_log2(log2(lin + 1))$scale, "log2")
})

test_that("metabolite names normalise the same way as the engine", {
  expect_equal(ery_norm_metabolite("Glucose-6-phosphate"), "glucose6phosphate")
  expect_equal(ery_norm_metabolite(" glucose 6 phosphate "), "glucose6phosphate")
  expect_true(is.na(ery_norm_metabolite("")))
})
