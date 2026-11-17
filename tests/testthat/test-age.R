## 红系衰老：三条路径都要对拍。
##
##   A 档（脂质分子种直接投影）—— 主读数，附留一阶段精度；
##   C 档（基因层走先验外推）—— 只在 A 档弃权时给，unverifiable；
##   拦截（身份不是红细胞）—— 连分数都不给。
##
## 置换 p 值同 attribution：numpy 与 R 的 RNG 不同流，所以 python 侧把实际
## 用到的排列导出（age_perms.csv.gz / age_feature_perms.csv.gz）。

fx <- function(...) testthat::test_path("fixtures", ...)
rd <- function(f) as.matrix(utils::read.delim(fx(f), row.names = 1, check.names = FALSE))
rperm <- function(f) as.matrix(utils::read.csv(gzfile(fx(f)), header = FALSE))
rjson <- function(f) jsonlite_free(f)

## 不引 jsonlite：只用它读一个扁平 JSON，用正则取标量就够了
jsonlite_free <- function(f) {
  txt <- paste(readLines(fx(f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  m <- regmatches(txt, gregexpr("\"([A-Za-z0-9_]+)\"\\s*:\\s*(\"[^\"]*\"|-?[0-9.]+|true|false|null)",
                                txt, perl = TRUE))[[1]]
  out <- list()
  for (h in m) {
    kv <- regmatches(h, regexec("\"([A-Za-z0-9_]+)\"\\s*:\\s*(.+)", h, perl = TRUE))[[1]]
    v <- kv[3]
    if (grepl("^\"", v)) v <- gsub("^\"|\"$", "", v)
    out[[kv[2]]] <- suppressWarnings(if (v %in% c("true", "false")) {
      v == "true"
    } else if (v == "null") {
      NULL
    } else if (!is.na(suppressWarnings(as.numeric(v)))) as.numeric(v) else v)
  }
  out
}

lipid_mat <- rd("age_lipid_input.tsv.gz")
age_perms <- rperm("age_perms.csv.gz")

test_that("the five reference stages land exactly on 0-4", {
  a <- ery_age(lipid_mat, permutations = age_perms)
  exp <- utils::read.csv(fx("age_expected.csv"), stringsAsFactors = FALSE)
  expect_true(a$ok)
  expect_equal(vapply(a$samples, `[[`, "", "sample"), exp$sample)

  for (i in seq_len(nrow(exp))) {
    s <- a$samples[[i]]
    expect_true(isTRUE(s$ok), info = exp$sample[i])
    expect_equal(s$position, exp$position[i], info = exp$sample[i])
    expect_equal(s$rho, exp$rho[i], tolerance = 1e-9, info = exp$sample[i])
    expect_equal(s$p_perm, exp$p_perm[i], tolerance = 1e-9, info = exp$sample[i])
    expect_equal(s$band[1], exp$band_lo[i], info = exp$sample[i])
    expect_equal(s$band[2], exp$band_hi[i], info = exp$sample[i])
    expect_equal(s$n_matched, exp$n_matched[i], info = exp$sample[i])
    expect_equal(s$n_within, exp$n_within[i], tolerance = 1e-9, info = exp$sample[i])
    expect_equal(s$p_value, exp$p_binom[i], tolerance = 1e-9, info = exp$sample[i])
    expect_equal(s$zone$key, exp$zone_key[i], info = exp$sample[i])
    expect_equal(s$tiers$A, exp$tier_A[i], info = exp$sample[i])
    expect_equal(s$tiers$B, exp$tier_B[i], info = exp$sample[i])
    ## 0–100 衰老评分必须与位置严格对应
    expect_equal(s$senescence_score, round(s$position / 4 * 100, 1))
  }
})

test_that("coverage, summary and the leave-one-stage error match the engine", {
  a <- ery_age(lipid_mat, permutations = age_perms)
  e <- rjson("age_model_expected.json")

  expect_equal(a$coverage$n_matched, e$n_matched)
  expect_equal(a$coverage$n_anchor, e$n_anchor)
  expect_equal(a$coverage$tier_A, e$tier_A)
  expect_equal(a$coverage$tier_B, e$tier_B)

  expect_equal(a$summary$n_ok, e$n_ok)
  expect_equal(a$summary$n_refused, e$n_refused)
  expect_equal(a$summary$median_position, e$median_position)
  expect_equal(a$summary$median_zone$key, e$median_zone)
  expect_equal(a$summary$n_young, e$n_young)
  expect_equal(a$summary$n_old, e$n_old)
  expect_equal(a$summary$n_between, e$n_between)
  expect_equal(a$summary$median_senescence_score, round(e$median_position / 4 * 100, 1))

  ## 留一阶段是唯一诚实的精度数字：自投影 1.00 不能当精度
  expect_equal(a$validation$loso_mae, e$loso_mae)
  expect_equal(a$validation$loso_max, e$loso_max)
  expect_equal(a$validation$loso_mae_pct, e$loso_mae_pct)
  expect_equal(a$validation$n_stages, 5L)
})

test_that("too few lipids gives insufficient_overlap, not a low score", {
  short <- lipid_mat[1:6, , drop = FALSE]
  a <- ery_age(short, permutations = age_perms)
  e <- rjson("age_refusal_expected.json")
  expect_false(a$ok)
  expect_equal(a$reason_code, "insufficient_overlap")
  expect_equal(a$reason_code, e$first_reason)
  expect_equal(length(a$samples), e$n_samples)
  ## 弃权样本不得带任何位置
  for (s in a$samples) expect_true(is.null(s$position))
  expect_true(grepl("12", a$reason))
})

test_that("a gene matrix abstains on the axis and falls back to priors", {
  g <- rd("gse63703_subset.tsv.gz")
  a <- ery_age(g, permutations = age_perms)
  e <- rjson("age_prior_markers_expected.json")
  expect_false(a$ok)
  expect_equal(a$reason_code, "no_lipid_names")
  expect_false(is.null(a$prior_extrapolation))
  expect_false(isTRUE(a$prior_extrapolation$ok))
  expect_equal(a$prior_extrapolation$reason_code, e$reason_code)
  expect_equal(a$prior_extrapolation$n_directional, e$n_directional)
  ## 弃权也要说清"命中几条"，否则用户不知道下一步补什么数据
  ## 先验表规模是数据驱动的（97 行），只断言它确实读到了
  expect_gte(a$prior_extrapolation$n_prior_rows, 50L)
})

test_that("the cross-layer tier reports the confidence it actually has", {
  p <- rd("age_prior_input.tsv.gz")
  a <- ery_age(p, permutations = age_perms)
  e <- rjson("age_prior_expected.json")
  pe <- a$prior_extrapolation
  expect_true(isTRUE(pe$ok))
  expect_true(isTRUE(pe$unverifiable))
  expect_equal(pe$grade, e$grade)
  expect_equal(pe$coverage, e$coverage, tolerance = 1e-6)
  expect_equal(pe$direction_agreement, e$direction_agreement, tolerance = 1e-6)
  expect_equal(pe$n_directional, e$n_directional)
  expect_equal(pe$n_hit, e$n_hit)
  expect_equal(pe$n_context_only, e$n_context_only)
  expect_equal(pe$confidence, e$confidence, tolerance = 1e-6)
  expect_equal(pe$feature_loo_r2, e$feature_loo_r2, tolerance = 1e-6)
  ## 跨层没有配对数据可验证 → 最高只到 low，绝不给 moderate/high
  expect_true(e$grade %in% c("low", "exploratory"))
  expect_true(all(nchar(pe$caveats$zh) > 10))
  ## 不给 p 值、不给区间、不给 MAE：这三个字段必须是 NULL
  expect_true(is.null(pe$p_value))
  expect_true(is.null(pe$ci))
  expect_true(is.null(pe$accuracy))
})

test_that("prior-extrapolated positions and scores match sample by sample", {
  p <- rd("age_prior_input.tsv.gz")
  a <- ery_age(p, permutations = age_perms)
  exp <- utils::read.csv(fx("age_prior_samples_expected.csv"), stringsAsFactors = FALSE)
  got <- a$prior_extrapolation$samples
  expect_equal(vapply(got, `[[`, "", "sample"), exp$sample)
  for (i in seq_along(got)) {
    expect_equal(got[[i]]$position, exp$position[i], tolerance = 1e-9)
    expect_equal(got[[i]]$score, exp$score[i], tolerance = 1e-6)
    expect_equal(got[[i]]$senescence_score, exp$senescence_score[i], tolerance = 1e-9)
    ## 期望文件里 zone 存成 Python dict 的字面量，只取其中的 key
    zk <- sub("^.*'key': '([^']+)'.*$", "\\1", exp$zone[i])
    expect_equal(got[[i]]$zone$key, zk)
  }
})

test_that("the feature model behind the extrapolation matches the engine", {
  perms <- rperm("age_feature_perms.csv.gz")
  fm <- ErythroWorkbench:::ery_age_feature_model(permutations = perms)
  e <- rjson("age_prior_expected.json")

  expect_equal(fm$n_species, 45L)
  expect_equal(fm$r2, 0.4355, tolerance = 1e-6)
  expect_equal(fm$loo_r2, 0.3185, tolerance = 1e-6)
  ## 只用特征权重能否把五个阶段排回正确顺序：这是"结构确实在描述成熟"的证据
  expect_equal(fm$stage_recovery_rho, 1, tolerance = 1e-6)
  ## 置换 p 用导出的同一组排列才可比
  ## 1000 次置换一次都没超过观测 R² → (1+0)/(1+1000)
  expect_equal(fm$permutation_p, 1 / 1001, tolerance = 1e-6)
  ## 1000 次置换里偶尔有落在 ss_tot=0 上的，会被剔除，所以用下界
  expect_gte(fm$n_perm, 990L)
  ## 碳数与双键数在 PC/PE 里高度共线，这个数要如实报出来
  expect_true(fm$chain_unsat_collinearity$PC > 0.8)
  expect_true(fm$chain_unsat_collinearity$PE > 0.8)
  ## directions 除 stage_marker 外都由面板实测，且被归一化到 [-1,1]
  expect_equal(fm$directions$stage_marker, 1)
  expect_true(all(abs(unlist(fm$directions)) <= 1 + 1e-9))
})

test_that("the identity gate refuses non-erythroid samples outright", {
  m <- lipid_mat[, seq_len(2), drop = FALSE]
  colnames(m) <- c("PBMC_1", "PBMC_2")
  a <- ery_age(m, permutations = age_perms)
  expect_false(a$ok)
  expect_equal(a$reason_code, "no_erythroid_samples")
  expect_equal(length(a$samples), 0L)
  expect_true(is.null(a$summary))
  ## 被拦的样本连分数都不给——给个低分是最坏的做法
  expect_true(grepl("\u975e\u7ea2\u7ec6\u80de", a$reason) ||
                grepl("not red cells", a$reason_en))
})

test_that("non-positive values are refused instead of silently closed", {
  bad <- lipid_mat
  bad[1, 1] <- -1
  a <- ery_age(bad, permutations = age_perms)
  ## RY 那一列含负值 → 该样本弃权并说明原因是尺度问题；其余样本照常
  s1 <- Filter(function(s) identical(s$sample, "RY"), a$samples)[[1]]
  expect_false(isTRUE(s1$ok))
  expect_equal(s1$reason, "requires_positive_values")
  expect_true(length(Filter(function(s) isTRUE(s$ok), a$samples)) >= 4L)
})
