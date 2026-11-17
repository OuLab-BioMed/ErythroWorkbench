## 通路指纹对拍。within 与 cohort 两种归一的分数必须分别对上，
## 因为它们量的不是同一件事（样本内相对 vs 相对本批）。

fix_dir <- function(...) test_path("fixtures", ...)
read_fixture_matrix <- function(name) {
  as.matrix(utils::read.delim(fix_dir(name), row.names = 1, check.names = FALSE))
}

test_that("pathway fingerprint reproduces the engine in both modes", {
  skip_if_not(file.exists(fix_dir("pathway_within_expected.csv")),
              "fixtures not generated")
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")

  for (md in c("within", "cohort")) {
    exp <- utils::read.csv(fix_dir(sprintf("pathway_%s_expected.csv", md)),
                           stringsAsFactors = FALSE, check.names = FALSE)
    got <- ery_pathways(m, mode = md)
    tab <- got$table
    expect_equal(tab$pathway_id, exp$pathway_id, info = md)
    expect_equal(tab$n_measured, exp$n_measured, info = md)
    expect_equal(tab$consistency, exp$consistency, tolerance = 1e-6, info = md)
    smp <- setdiff(names(exp), c("pathway_id", "n_members", "n_measured", "consistency"))
    for (s in smp) {
      expect_equal(tab[[s]], exp[[s]], tolerance = 1e-6,
                   info = paste(md, s))
    }
  }
})

test_that("abundance carriers are excluded from scoring but kept in curation", {
  ## HBB 这类丰度载体不参与打分：单基因就能把整条通路的均值 z 抬到 26.9。
  ## 但策展成员数（n_members）里仍然算它，供机制注释使用。
  scorable <- ErythroWorkbench:::.pathway_members("gene", scoring_only = TRUE)
  allm <- ErythroWorkbench:::.pathway_members("gene", scoring_only = FALSE)
  ## HBB / ALB / CA1 都在 no_metabolism，role 标着 abundance-excluded
  nem_all <- allm[["no_metabolism"]]
  nem_sc <- scorable[["no_metabolism"]]
  expect_true("HBB" %in% nem_all)
  expect_false("HBB" %in% nem_sc)
  expect_true(length(nem_all) > length(nem_sc))
})

test_that("pathway ranking multiplies by sqrt(n_measured)", {
  ## 不乘 sqrt(n) 时，成员最少的通路在纯噪声下也有一半概率排第一。
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_pathways(m, mode = "within")
  tab <- got$table
  smp <- got$samples
  key <- apply(tab[, smp, drop = FALSE], 1, function(v) max(abs(v))) *
    sqrt(tab$n_measured)
  expect_equal(order(-key), seq_len(nrow(tab)))
})

test_that("unscorable pathways are listed, never silently dropped", {
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_pathways(m, mode = "within")
  allm <- ErythroWorkbench:::.pathway_members("gene", scoring_only = TRUE)
  expect_equal(sort(c(got$table$pathway_id, got$excluded$pathway_id)),
               sort(names(allm)))
  ## 缺席必须给出可核对的理由
  if (nrow(got$excluded)) {
    expect_true(all(got$excluded$reason %in%
                      c("too_few_scorable_members", "not_measured_in_this_matrix")))
  }
})

test_that("cohort mode centres each feature, so a pathway mean across samples is zero", {
  ## 这是通路指纹最容易读错的地方：cohort 模式下每条通路跨样本平均恒为 0，
  ## 画「平均偏移」柱图必然是一排 0——那不代表没有变化。
  m <- read_fixture_matrix("gse63703_subset.tsv.gz")
  got <- ery_pathways(m, mode = "cohort")
  tab <- got$table
  smp <- got$samples
  means <- rowMeans(tab[, smp, drop = FALSE])
  expect_true(all(abs(means) < 0.05))
})
