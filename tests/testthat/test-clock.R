## 时钟对拍。两条路径分别验：
##   同层模型（miRNA / 代谢组）与跨层模块桥（代谢组 → 任意层）。
## 注意 fixture 用的数据集本身就是各自时钟的训练集，所以误差偏乐观——
## 这里只验证 R 与引擎在同样输入下给出同样天数，不验证准确度。

fix_dir <- function(...) test_path("fixtures", ...)

read_fixture_matrix <- function(name) {
  d <- utils::read.delim(fix_dir(name), row.names = 1, check.names = FALSE)
  as.matrix(d)
}

cmp_clock <- function(got, exp) {
  expect_equal(got$sample, exp$sample)
  expect_equal(got$day, exp$day, tolerance = 1e-6)
  ## lo/hi 是「预测 ± 分位残差」后舍入到一位小数，边界上 R 与 numpy 的
  ## round() 可能差一个末位（银行家舍入遇上浮点表示），放 0.11 容忍一个末位。
  expect_equal(got$lo, exp$lo, tolerance = 0.11)
  expect_equal(got$hi, exp$hi, tolerance = 0.11)
  expect_equal(got$raw, exp$raw, tolerance = 1e-6)
  expect_equal(got$clipped, as.logical(exp$clipped))
}

test_that("miRNA clock reproduces the engine", {
  skip_if_not(file.exists(fix_dir("clock_mirna_expected.csv")),
              "fixtures not generated; run data-raw/make_fixtures.py")
  m <- read_fixture_matrix("gse114990_mirna.tsv.gz")
  exp <- utils::read.csv(fix_dir("clock_mirna_expected.csv"), stringsAsFactors = FALSE)
  got <- ery_clock(m, space = "mirna")
  cmp_clock(got, exp)
  expect_equal(attr(got, "mode"), "reference")
  expect_equal(attr(got, "model"), "mirna")
})

test_that("metabolite clock reproduces the engine", {
  skip_if_not(file.exists(fix_dir("clock_metabolite_expected.csv")),
              "fixtures not generated")
  m <- read_fixture_matrix("st004148_metabolome.tsv.gz")
  exp <- utils::read.csv(fix_dir("clock_metabolite_expected.csv"), stringsAsFactors = FALSE)
  got <- ery_clock(m, space = "metabolite")
  cmp_clock(got, exp)
  expect_equal(attr(got, "model"), "metabolite")
  expect_equal(attr(got, "day_range"), c(1, 49))
})

test_that("the cross-layer module bridge reproduces the engine", {
  skip_if_not(file.exists(fix_dir("clock_module_expected.csv")),
              "fixtures not generated")
  m <- read_fixture_matrix("st004148_metabolome.tsv.gz")
  exp <- utils::read.csv(fix_dir("clock_module_expected.csv"), stringsAsFactors = FALSE)
  got <- ery_clock(m, space = "metabolite", path = "module")
  cmp_clock(got, exp)
  expect_equal(attr(got, "mode"), "extrapolated")
  expect_true(attr(got, "matched_modules") >= 3)
})

test_that("true days are used to validate, never to fit", {
  m <- read_fixture_matrix("st004148_metabolome.tsv.gz")
  base <- ery_clock(m, space = "metabolite")
  ## 给一份"真实天数"：应该只改变验证块，不改变预测本身
  fake <- stats::setNames(rep(7, nrow(base)), base$sample)
  v <- ery_clock(m, space = "metabolite", days = fake)
  expect_equal(v$day, base$day)
  expect_equal(attr(v, "mode"), "reference_verified")
  expect_true(isTRUE(attr(v, "validation")$ok))
  ## 全部真值相同 → 偏倚等于平均误差，但只在给出天数的样本上算
  okk <- is.finite(v$day)
  expect_true(sum(okk) >= 4)
  expect_equal(attr(v, "validation")$bias_days,
               round(mean(v$day[okk] - 7), 2), tolerance = 0.01)
  expect_equal(attr(v, "validation")$n, sum(okk))
})

test_that("storage days are clipped to the trained range and flagged", {
  m <- read_fixture_matrix("st004148_metabolome.tsv.gz")
  got <- ery_clock(m, space = "metabolite")
  rng <- attr(got, "day_range")
  expect_true(all(got$day >= rng[1] & got$day <= rng[2], na.rm = TRUE))
  ## 被截断的样本必须标出来，否则看不出这是外推
  expect_true(is.logical(got$clipped))
  expect_equal(got$clipped, abs(got$raw - got$day) > 0.5)
})

test_that("robust row z falls back to SD only when MAD is zero", {
  ## 中位 + 1.4826*MAD；MAD 为 0 时退到标准差；两者都退化则该特征无信息。
  ## 兜底比分母设成 1 会让该特征保留原始量级，模块分数被推到 1e5。
  x <- rbind(a = c(1, 2, 3, 4, 5), b = rep(2, 5), c = c(1, 1, 1, 1, 5))
  rz <- ErythroWorkbench:::.robust_row_z(x)
  expect_equal(unname(rz$ok[1]), TRUE)
  ## b 行 MAD = 0 且 SD = 0 → 无信息
  expect_equal(unname(rz$ok[2]), FALSE)
  ## c 行 MAD = 0 但 SD > 0 → 退到 SD
  expect_equal(unname(rz$ok[3]), TRUE)
  expect_true(max(abs(rz$z[3, ])) <= 5)
})
