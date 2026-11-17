## 条件谱匹配：与引擎逐位对拍。
##
## 三件事都要对：主排序（best_rho 而非条件均值）、未舍入的 rho、以及置换 p 值。
## p 值只有在喂同一组置换时才可比——numpy 的 PCG64 与 R 的 MT 不可能同流，
## 所以 python 侧把实际用到的置换导出来（attribution_perms_*.csv.gz）。

fx <- function(...) testthat::test_path("fixtures", ...)

## CSV 里的中文是 UTF-8 字节，但 read.csv 在非 UTF-8 locale 下会把它们标成
## "unknown"；与包装载后的字符串（标记 UTF-8）按 == 比较会失败，尽管字节完全
## 相同。强制标记即可——不要 enc2utf8()，那会把已经是 UTF-8 的字节再转一次。
u8 <- function(x) { Encoding(x) <- "UTF-8"; x }

## CSV 里的中文是 UTF-8 字节，但 read.csv 在非 UTF-8 locale 下会把它们标成
## "unknown"；与包装载后的字符串（标记 UTF-8）按 == 比较会失败，尽管字节完全
## 相同。强制标记即可——不要 enc2utf8()，那会把已经是 UTF-8 的字节再转一次。
u8 <- function(x) { Encoding(x) <- "UTF-8"; x }
rd <- function(f) as.matrix(utils::read.delim(fx(f), row.names = 1, check.names = FALSE))
rperm <- function(f) as.matrix(utils::read.csv(gzfile(fx(f)), header = FALSE))

m <- rd("gse63703_subset.tsv.gz")

attr_run <- function() {
  ## 样本 1（RBC_1）与样本 4（PBMC_1）各带一组置换，其余用 R 自己的 RNG
  pl <- vector("list", ncol(m))
  pl[[1]] <- rperm("attribution_perms_0.csv.gz")
  pl[[4]] <- rperm("attribution_perms_3.csv.gz")
  ery_attribution(m, permutations = pl)
}

test_that("attribution reports the reference library shape it actually used", {
  a <- attr_run()
  expect_equal(a$space, "gene")
  expect_equal(a$n_samples, 8L)
  expect_equal(a$n_centroids, 34L)          # 本 id 空间里共享特征够多的红系画像
  expect_equal(a$n_reference_centroids, 1L) # 细胞类型参照只有 1 张
  expect_equal(a$n_centroids_excluded, 7L)  # mixed 角色：两边证据都不够，不进判别
  expect_equal(a$n_conditions, 15L)
  expect_equal(length(a$results), 8L)
  expect_equal(length(a$caveats$zh), 6L)
  expect_equal(length(a$caveats$en), 6L)
})

test_that("per-sample condition call, rounded and raw, matches the engine", {
  a <- attr_run()
  exp <- utils::read.csv(fx("attribution_expected.csv"), stringsAsFactors = FALSE)
  exp$top_label <- u8(exp$top_label)
  exp$runner_label <- u8(exp$runner_label)
  raw <- attr(a, "raw")
  expect_equal(raw$sample, exp$sample)

  for (i in seq_len(nrow(exp))) {
    r <- a$results[[i]]
    expect_true(isTRUE(r$ok), info = exp$sample[i])
    expect_equal(r$top$label, exp$top_label[i], info = exp$sample[i])
    expect_equal(r$top$best_rho, exp$best_rho[i], info = exp$sample[i])
    expect_equal(r$top$mean_rho, exp$mean_rho[i], info = exp$sample[i])
    ## 中位数落在 .0005 边界时，numpy 与 R 的 round 会差一个末位（0.2875 在
    ## 二进制里是 0.2874999…）。raw 那一组断言才是严格比对，这里放宽一个末位。
    ## expect_equal 的 tolerance 是**相对**容差；这里是绝对值差一个末位，
    ## 用 expect_lt 才表达得准
    expect_lt(abs(r$background_median_rho - exp$background_median_rho[i]), 0.0011)
    expect_equal(r$n_cohorts_compared, exp$n_compared[i], info = exp$sample[i])
    expect_equal(r$n_conditions, exp$n_conditions[i], info = exp$sample[i])
    expect_equal(length(r$conditions), min(8L, exp$n_conditions[i]),
                 info = exp$sample[i])
    expect_equal(r$runner_up$label, exp$runner_label[i], info = exp$sample[i])
    expect_equal(r$runner_up$best_rho, exp$runner_best_rho[i], info = exp$sample[i])
    expect_equal(r$top$consistency, exp$consistency[i], info = exp$sample[i])
    expect_equal(r$identity$compartment, exp$identity[i], info = exp$sample[i])

    ## 未舍入值：这是能把"浮点末位"和"算法真的不一样"分开的唯一办法
    expect_equal(raw$best_rho_raw[i], exp$best_rho_raw[i], tolerance = 1e-9,
                 info = exp$sample[i])
    expect_equal(raw$mean_rho_raw[i], exp$mean_rho_raw[i], tolerance = 1e-9,
                 info = exp$sample[i])
    expect_equal(raw$bg_raw[i], exp$bg_raw[i], tolerance = 1e-9,
                 info = exp$sample[i])
  }
})

test_that("permutation p-values match when the same permutations are supplied", {
  a <- attr_run()
  exp <- utils::read.csv(fx("attribution_expected.csv"), stringsAsFactors = FALSE)
  ## 只有样本 1 与 4 用了导出的置换，其余样本的 p 依赖各自 RNG，不作断言
  for (i in c(1L, 4L)) {
    expect_equal(a$results[[i]]$top$p, exp$p_value[i], tolerance = 1e-9,
                 info = exp$sample[i])
    expect_equal(a$results[[i]]$top$null_q95, exp$null_q95[i], tolerance = 1e-9,
                 info = exp$sample[i])
  }
  expect_equal(a$results[[1]]$top$margin, exp$margin[1], tolerance = 1e-9)
})

test_that("conditions are ranked by the strongest cohort, not by the class mean", {
  ## 这条是模块最容易改错的地方：GSE301441 的健康样本与自己组的 rho 是 0.97，
  ## 但被 0.63 的终末组拉到均值 0.61，按均值排会输给 MM 的 0.83。
  a <- attr_run()
  for (r in a$results) {
    if (!isTRUE(r$ok)) next
    best <- vapply(r$conditions, `[[`, 0, "best_rho")
    expect_equal(best, sort(best, decreasing = TRUE), info = r$sample)
    expect_true(r$conditions[[1]]$best_rho >= r$conditions[[1]]$mean_rho)
    ## 一致性 = 均值 / 最强，只能是 [0,1] 之间的数
    expect_true(r$top$consistency >= 0 && r$top$consistency <= 1)
  }
})

test_that("a non-erythroid sample carries the identity warning", {
  a <- attr_run()
  pbmc <- Filter(function(r) grepl("^PBMC", r$sample), a$results)
  expect_length(pbmc, 3L)
  for (r in pbmc) {
    expect_equal(r$identity$compartment, "leukocyte")
    expect_false(is.null(r$identity_warning))
    expect_true(grepl("\u9519\u4f4d\u63d0\u95ee", r$identity_warning))
  }
  ## 红细胞与储存样本不该被加这条警告
  rbc <- Filter(function(r) grepl("^RBC", r$sample), a$results)
  for (r in rbc) expect_true(is.null(r$identity_warning))
})

test_that("attribution abstains instead of guessing when nothing shares features", {
  ## 只有 3 个基因、且都不在画像库的 id 空间里 → 没有任何 centroid 可用
  tiny <- matrix(1:6, nrow = 3,
                 dimnames = list(c("ZZZ1", "ZZZ2", "ZZZ3"), c("s1", "s2")))
  a <- ery_attribution(tiny, n_perm = 10)
  r <- a$results[[1]]
  expect_false(isTRUE(r$ok))
  expect_equal(r$reason, "no_centroid_with_enough_shared_features")
})

test_that("metabolite input stays on the metabolite space", {
  ## 跨空间不外推：基因侧画像不能拿来解释代谢物谱。
  ## 命名要按库里代谢物的实际写法（化学名），否则 id 空间启发式会判成基因。
  m2 <- matrix(c(100, 200, 50, 300, 120, 80), nrow = 3,
               dimnames = list(c("L-alanine", "L-arginine", "L-glutamine"),
                               c("s1", "s2")))
  a <- ery_attribution(m2, n_perm = 10)
  expect_equal(a$space, "metabolite")
  ## 3 个代谢物达不到 MIN_SHARED（30），所以没有任何画像可比
  expect_equal(a$n_centroids, 0L)
  expect_true(all(vapply(a$results, function(r) !isTRUE(r$ok), TRUE)))
})
