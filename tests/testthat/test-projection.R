## 投影对拍。注意这份 fixture 的样本本身就在参考图谱里（GSE63703 是图谱
## 来源之一），所以 best_rho 会到 1.0——这是"自投影"，引擎 CAVEATS 里明确
## 警告过：准确率必须看留一**数据集**验证，不能看这种样本内相似度。
## 这里只验证 R 与引擎在同样输入下给出同样的标签与弃权理由。

fix_dir <- function(...) test_path("fixtures", ...)

read_fixture_matrix <- function(name) {
  d <- utils::read.delim(fix_dir(name), row.names = 1, check.names = FALSE)
  as.matrix(d)
}

test_that("projection reproduces the engine's labels and abstentions", {
  skip_if_not(file.exists(fix_dir("projection_expected.csv")),
              "fixtures not generated; run data-raw/make_fixtures.py")
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  exp <- utils::read.csv(fix_dir("projection_expected.csv"), stringsAsFactors = FALSE)

  got <- ery_transfer(m)

  expect_equal(got$sample, exp$sample)
  ## 标签：空串在 fixture 里表示弃权
  expect_equal(ifelse(is.na(got$label), "", got$label), exp$label)
  expect_equal(ifelse(is.na(got$abstain), "", got$abstain), exp$abstain)
  expect_equal(got$n_markers, exp$n_markers)
  expect_equal(got$best_rho, exp$best_rho, tolerance = 2e-3)

  ## 概率与熵只在给出标签时有值
  keep <- nzchar(exp$label)
  expect_equal(got$probability[keep], exp$probability[keep], tolerance = 2e-3)
})

test_that("both gates come from measured nulls, not constants", {
  atlas <- ery_atlas()
  ## OOD 闸门 = 图谱自比对秩相似度的 p5
  expect_equal(unname(atlas$ood_gate$rank_p5), 0.6549, tolerance = 1e-4)
  ## margin 闸门 = 留一数据集投票 margin 的 p10；逐对重排秩才对得上 0.4334
  expect_equal(unname(atlas$margin_gate_rank), 0.4334, tolerance = 5e-4)
  ## 这两个数必须来自数据，写死成 0.6 / 0.2 之类的常数会让弃权率失控
  expect_gt(atlas$margin_gate_rank, 0)
  expect_lt(atlas$margin_gate_rank, 1)
})

test_that("the projection is a projection: abstention is explicit", {
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_transfer(m)
  ## 弃权的样本不给标签，而且理由是 OOD / ambiguous / no_overlap 之一
  ab <- got[is.na(got$label), ]
  if (nrow(ab)) {
    expect_true(all(ab$abstain %in% c("ood", "ambiguous", "no_overlap",
                                      "no_positive_similarity")))
  }
  ## 覆盖率与回答数应自洽
  expect_equal(attr(got, "n_answered"), sum(!is.na(got$label)))
  expect_equal(attr(got, "coverage"), round(mean(!is.na(got$label)), 3))
})

test_that("the atlas rank matrix keeps missing markers in place", {
  ## 曾经的 bug：用 unlist() 读 JSON，null 会被整个丢掉，后面的值左移、
  ## 末尾补 NA。缺失总数一样，位置全错，相似度算在错位向量上且不报错。
  atlas <- ery_atlas()
  r1 <- atlas$rank[1, ]
  ## GSE108378 RBC1 在引擎里的缺失位置是第 3、14 位与第 25-46 位（0-based）
  expect_equal(unname(which(is.na(r1))) - 1L, c(2L, 13L, 24:45))
  expect_equal(unname(atlas$rank[1, 1]), 1)
  expect_equal(unname(atlas$rank[1, 2]), 0.9942, tolerance = 1e-3)
  expect_true(is.na(atlas$rank[1, 3]))
})
