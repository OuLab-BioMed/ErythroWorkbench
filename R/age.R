## ---------------------------------------------------------------------------
## 红系衰老：整条链路的对外入口
##
## 三档，按证据强度从高到低，**不并列显示**：
##
##   A 档（脂质分子种直接投影）—— 主读数，0–100 衰老评分，附留一阶段精度；
##   B 档（转录/蛋白层走「脂质特征先验」外推）—— 只在 A 档弃权时给，
##        标注 unverifiable，不给 p 值、不给区间、最高只到 low；
##   拦截（身份不是红细胞）—— 连分数都不给。给个低分是最坏的选项：
##        它会被人当成"年轻的红细胞"引用。
## ---------------------------------------------------------------------------

.ery_ery_fraction_min <- 0.80

## ---- 身份闸门 -------------------------------------------------------------
## 脂质层没有 RNA，算不了红系占比，所以用注释 / 样本名的细胞身份判定。
## 引擎还有一条「红系 RNA 占比 ≥ 0.80」的分支，它需要端点混合模型的结果；
## 本包未移植 fraction 端点模型，因此 gene 空间的输入在这里会明确落到
## unverifiable，而不是假装核验过。
.age_gate <- function(mat, annotation = NULL) {
  id <- ery_identify(mat, annotation)
  comp <- id$compartment
  if (any(comp != "unknown")) {
    ok <- !(comp %in% .si_block)
    details <- lapply(seq_along(comp), function(i) {
      list(sample = id$sample[i], compartment = comp[i], passed = ok[i])
    })
    return(list(
      threshold = .ery_ery_fraction_min, level = "identity",
      n_input = length(comp), n_after = sum(ok),
      passed = id$sample[ok], blocked = id$sample[!ok], details = details,
      note = "\u672c\u5c42\u6ca1\u6709 RNA\uff0c\u65e0\u6cd5\u7b97\u7ea2\u7cfb\u5360\u6bd4\uff1b\u6539\u7528\u6ce8\u91ca / \u6837\u672c\u540d\u7684\u7ec6\u80de\u8eab\u4efd\u5224\u5b9a\uff0c\u660e\u786e\u4e3a\u767d\u7ec6\u80de / \u8840\u6d46 / \u5176\u4ed6\u7ec4\u7ec7\u7684\u6837\u672c\u88ab\u62e6\u622a\u3002",
      note_en = paste("This layer has no RNA, so an erythroid fraction cannot be computed;",
                      "annotation/sample-name identity is used instead and leukocyte, plasma",
                      "and other-tissue samples are blocked.")
    ))
  }
  list(
    threshold = .ery_ery_fraction_min, level = "unverifiable",
    n_input = length(comp), n_after = length(comp),
    passed = id$sample, blocked = character(0), details = list(),
    note = paste("\u65e2\u65e0\u53ef\u7528\u7684\u7ea2\u7cfb\u5360\u6bd4\uff0c\u4e5f\u6ca1\u6709\u53ef\u5224\u8bfb\u7684\u7ec6\u80de\u8eab\u4efd\u6ce8\u91ca\uff1a\u672c\u6279\u6837\u672c\u8eab\u4efd**\u672a\u6838\u9a8c**\uff0c",
                 "\u7ed3\u679c\u53ea\u80fd\u5f53\u4f5c\u63a2\u7d22\u6027\u8bfb\u6570\u3002\u8102\u8d28\u5206\u5b50\u79cd\u9762\u677f\u672c\u8eab\u7684\u8f68\u8ff9\u76f8\u5173\u68c0\u9a8c\u662f\u53e6\u4e00\u9053\u72ec\u7acb\u95f8\u95e8\u3002"),
    note_en = paste("Neither an erythroid fraction nor readable identity annotation is",
                    "available: identity is UNVERIFIED and the result is exploratory only.",
                    "The trajectory-correlation test on the lipid panel is a separate gate.")
  )
}

## ---- 跨层先验外推（C 档）--------------------------------------------------
## 先验表：基因 → 它推动哪一类脂质特征（链长 / 不饱和度 / 类别占比 / 阶段标志）。
## effect=0 或 tier=C 的条目只贡献覆盖，不贡献方向。
.age_prior_table <- function() {
  if (!is.null(.ewb_cache$age_priors)) return(.ewb_cache$age_priors)
  p <- .ewb_extdata("lipid_regulators.tsv")
  tab <- utils::read.delim(p, stringsAsFactors = FALSE, quote = "",
                           check.names = FALSE, fileEncoding = "UTF-8-BOM")
  tab$effect <- as.integer(tab$effect)
  tab$tier_w <- unname(c(A = 1.0, B = 0.6, C = 0.0)[tab$tier])
  tab$tier_w[is.na(tab$tier_w)] <- 0
  .ewb_cache$age_priors <- tab
  tab
}

.age_prior_sources <- function() {
  p <- .ewb_extdata("lipid_regulator_sources.tsv")
  utils::read.delim(p, stringsAsFactors = FALSE, quote = "",
                    check.names = FALSE, fileEncoding = "UTF-8-BOM")
}

## 最小二乘：用 SVD 伪逆，与 numpy.linalg.lstsq(rcond = None) 同一个解。
## 不用 qr.solve：矩阵奇异时它会直接报错，而 numpy 给最小范数解——两边行为
## 不一致会让 r2 在退化输入上分叉。
.lstsq <- function(X, y) {
  s <- svd(X)
  tol <- .Machine$double.eps * max(dim(X)) * s$d[1]
  dinv <- ifelse(s$d > tol, 1 / s$d, 0)
  as.numeric(s$v %*% (dinv * crossprod(s$u, y)))
}

## 设计矩阵：截距 / 归一化碳数 / 归一化双键数 / PE 的碳数交互 / PC 的双键交互。
## 交互项正是「PE 34/36 升、38/40 降」与「PC 多不饱和升」允许有类别特异方向。
.age_design <- function(keys) {
  cn <- (keys$c - 36) / 4
  dn <- (keys$db - 2) / 2
  is_pe <- as.numeric(keys$cls == "PE")
  is_pc <- as.numeric(keys$cls == "PC")
  cbind(1, cn, dn, cn * is_pe, dn * is_pc)
}

#' Feature model behind the cross-layer extrapolation
#'
#' Measures, on the 45-species five-stage panel, which structural features move
#' with maturation: chain length, unsaturation and class share. The three are
#' highly collinear inside PC and PE (r(c, db) about 0.85), so the per-feature
#' directions are reported with that collinearity attached rather than
#' pretending the two can be separated by a regression coefficient.
#'
#' @param permutations Optional matrix of 1-based permutations of the anchor
#'   rows used for the null distribution. Supplying it makes `permutation_p`
#'   reproducible across implementations, which is how the parity test compares
#'   this against the Python engine.
#' @return A list with `beta`, `r2`, `loo_r2`, `permutation_p`,
#'   `stage_recovery_rho`, `class_direction`, `class_feature_direction`,
#'   `chain_unsat_collinearity`, `groups` and `directions`.
#' @keywords internal
ery_age_feature_model <- function(permutations = NULL) {
  if (is.null(permutations) && !is.null(.ewb_cache$age_feature_model)) {
    return(.ewb_cache$age_feature_model)
  }
  m <- ery_lipid_model()
  keys <- m$keys
  Z <- m$z
  P <- m$pct
  stage <- seq_len(length(m$stages)) - 1

  d <- vapply(seq_len(nrow(keys)), function(i) .pearson(Z[i, ], stage), 0)
  X <- .age_design(keys)

  beta <- .lstsq(X, d)
  fitted <- as.numeric(X %*% beta)
  ss_res <- sum((d - fitted)^2)
  ss_tot <- sum((d - mean(d))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NaN

  ## 留一分子种：拟合优度不能只报 in-sample R²
  loo <- numeric(length(d))
  for (i in seq_along(d)) {
    keep <- seq_along(d) != i
    b <- .lstsq(X[keep, , drop = FALSE], d[keep])
    loo[i] <- sum(X[i, ] * b)
  }
  loo_r2 <- if (ss_tot > 0) 1 - sum((d - loo)^2) / ss_tot else NaN

  ## 置换：把 (类别, 碳数, 双键) 在分子种间打乱，检验结构是否优于随机分配。
  ## 这个 p 只检验脂质面板内部的特征结构，不能当跨层预测的显著性。
  n_perm <- 1000L
  perm_r2 <- numeric(n_perm)
  if (is.null(permutations)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else NULL
    set.seed(20261112L)
  }
  for (t in seq_len(n_perm)) {
    idx <- if (is.null(permutations)) sample.int(length(d)) else permutations[t, ]
    b <- .lstsq(X[idx, , drop = FALSE], d)
    f <- as.numeric(X[idx, , drop = FALSE] %*% b)
    ## 必须显式写 NaN：只写 if(ss_tot>0) 的话 else 分支会保留初始的 0，
    ## 0 是一个合法的 R²，会被算进置换分布里，把 p 值往下压
    ## 必须显式写 NaN：只写 if(ss_tot>0) 的话 else 分支会保留初始的 0，
    ## 0 是一个合法的 R²，会被算进置换分布里，把 p 值往下压
    perm_r2[t] <- if (ss_tot > 0) 1 - sum((d - f)^2) / ss_tot else NaN
  }
  if (is.null(permutations) && exists(".Random.seed", envir = globalenv())) {
    if (is.null(old_seed)) rm(".Random.seed", envir = globalenv())
    else assign(".Random.seed", old_seed, envir = globalenv())
  }
  perm_r2 <- perm_r2[is.finite(perm_r2)]
  perm_p <- if (length(perm_r2)) (1 + sum(perm_r2 >= r2)) / (1 + length(perm_r2)) else NULL

  ## 只用特征权重（不用逐分子种方向）能否把五个阶段排回正确顺序
  feat_score <- vapply(seq_along(stage), function(t) sum(fitted * Z[, t]), 0)
  recovery <- .spearman(stage, feat_score)

  cls_dir <- list(); cls_delta <- list()
  for (cl in c("PC", "PE", "SM")) {
    idx <- which(keys$cls == cl)
    if (!length(idx)) next
    tot <- colSums(P[idx, , drop = FALSE])
    cls_dir[[cl]] <- .pearson(tot, stage)
    cls_delta[[cl]] <- round(tot[length(tot)] - tot[1], 3)
  }

  ## 用户提出的两条规律，用最直白的方式实测（描述性，5 个合并样本）
  group_delta <- function(sel) {
    v <- P[sel, ncol(P)] - P[sel, 1]
    list(n = length(v),
         mean_delta_pct = if (length(v)) round(mean(v), 3) else NULL,
         n_up = sum(v > 0), n_down = sum(v < 0))
  }
  groups <- list(
    pc_db_le1 = group_delta(keys$cls == "PC" & keys$db <= 1),
    pc_db_ge2 = group_delta(keys$cls == "PC" & keys$db >= 2),
    pe_c34_36 = group_delta(keys$cls == "PE" & keys$c <= 36),
    pe_c38_40 = group_delta(keys$cls == "PE" & keys$c >= 38)
  )

  ## 按类别实测「高链长组 − 低链长组」「高双键组 − 低双键组」的成熟方向差。
  ## 不用回归系数：PC/PE 面板里碳数与双键数高度共线，最小二乘会把同一段信号
  ## 随意分给两个共线列。分组对比至少不假装能把两者分离。
  class_feature <- list(); collinearity <- list()
  for (cl in c("PC", "PE", "SM")) {
    idx <- which(keys$cls == cl)
    if (!length(idx)) next
    cc <- keys$c[idx]; db <- keys$db[idx]; dd <- d[idx]
    class_feature[[cl]] <- list(
      chain = mean(dd[cc >= stats::median(cc)]) - mean(dd[cc < stats::median(cc)]),
      unsat = mean(dd[db >= stats::median(db)]) - mean(dd[db < stats::median(db)])
    )
    collinearity[[cl]] <- round(.pearson(cc, db), 3)
  }

  n_cls <- vapply(c("PC", "PE", "SM"),
                  function(cl) sum(keys$cls == cl), 0L)
  tot_n <- sum(n_cls)
  chain_net <- sum(n_cls * vapply(names(n_cls), function(cl) class_feature[[cl]]$chain, 0)) / tot_n
  unsat_net <- sum(n_cls * vapply(names(n_cls), function(cl) class_feature[[cl]]$unsat, 0)) / tot_n

  raw <- list(
    chain = chain_net, unsat = unsat_net,
    pc_frac = cls_dir$PC %||% 0, pe_frac = cls_dir$PE %||% 0,
    sm_frac = cls_dir$SM %||% 0,
    ## 单位标尺：该 feature 的数值就是成熟方向本身，先验表用 effect=-1 表示
    ## VPS13A / CD71 高表达推向年轻端
    stage_marker = 1.0
  )
  fin <- vapply(raw, function(v) is.finite(v), TRUE)
  mx <- if (any(fin)) max(abs(unlist(raw[fin]))) else 1
  if (!is.finite(mx) || mx == 0) mx <- 1
  dirs <- lapply(raw, function(v) if (is.finite(v)) as.numeric(v / mx) else 0)
  dirs$context <- 0

  out <- list(
    beta = round(as.numeric(beta), 4),
    beta_labels = c("intercept", "chain", "unsat", "chain_x_PE", "unsat_x_PC"),
    r2 = round(r2, 4), loo_r2 = round(loo_r2, 4),
    permutation_p = perm_p, n_perm = length(perm_r2),
    stage_recovery_rho = round(recovery, 4),
    class_direction = cls_dir, class_delta_pct = cls_delta,
    class_feature_direction = class_feature,
    chain_unsat_collinearity = collinearity,
    groups = groups, directions = dirs, n_species = nrow(keys),
    note = paste("directions \u9664 stage_marker \u5916\u5168\u90e8\u7531 45 \u5206\u5b50\u79cd\u4e94\u9636\u6bb5\u9762\u677f\u5b9e\u6d4b\uff1b",
                 "stage_marker=+1 \u53ea\u662f\u5355\u4f4d\u6807\u5c3a\uff0c\u5148\u9a8c\u8868\u7528 effect=-1 \u8868\u8fbe VPS13A/CD71",
                 "\u968f\u6210\u719f\u4e0b\u964d\uff08\u9ad8\u8868\u8fbe\u63a8\u5411\u5e74\u8f7b\u7aef\uff09\u3002PC/PE \u7684\u78b3\u6570\u4e0e\u53cc\u952e\u6570\u9ad8\u5ea6\u5171\u7ebf\uff0c",
                 "\u4e24\u4e2a feature \u7684\u65b9\u5411\u4e0d\u53ef\u72ec\u7acb\u89e3\u91ca\u3002")
  )
  if (is.null(permutations)) .ewb_cache$age_feature_model <- out
  out
}

.age_prior_caveats <- function(abstain) {
  zh <- c(
    "\u673a\u5236\u8fb9\u754c\uff1aMinetti 2025\uff08PMID 40016214\uff09\u6307\u51fa\u7f51\u7ec7\u7ea2\u4e0e\u6210\u719f\u7ea2\u7ec6\u80de\u4e0d\u80fd\u4ece\u5934\u5408\u6210\u78f7\u8102\uff0c\u819c\u8102\u91cd\u5851\u6765\u81ea\u9009\u62e9\u6027\u53bb\u9664\u4e0e\u8840\u6d46\u8102\u8d28\u4ea4\u6362\uff08VPS13A \u53ef\u80fd\u53c2\u4e0e\uff09\u3002\u628a ELOVL/FADS/LPCAT \u7b49\u9176\u7684\u8868\u8fbe\u63a5\u5230\u7ea2\u7ec6\u80de\u819c\u8102\u94fe\u957f\u4e0a\uff0c\u662f**\u95f4\u63a5\u5047\u8bf4**\uff0c\u4e0d\u662f\u5df2\u8bc1\u5b9e\u7684\u8c03\u63a7\u673a\u5236\u3002",
    "0\u2013100 \u8870\u8001\u8bc4\u5206\u662f**\u672c\u6279\u6837\u672c\u5185\u7684\u76f8\u5bf9\u79e9\u6295\u5f71**\uff08\u4e94\u9636\u6bb5\u8f68\u8ff9\u53ea\u662f\u523b\u5ea6\uff09\uff0c\u4e0d\u662f\u7edd\u5bf9\u50a8\u5b58\u5929\u6570\uff0c\u4e5f\u4e0d\u662f\u53ef\u8de8\u6279\u6b21\u6bd4\u8f83\u7684\u5e74\u9f84\u5750\u6807\u3002\u56e0\u4e3a\u662f\u79e9\u6295\u5f71\uff0c\u4e24\u7aef\u6837\u672c\u5fc5\u7136\u843d\u5728\u63a5\u8fd1 0 \u548c 100 \u7684\u5206\u6570\u4e0a\uff0c**\u8fd9\u4e0d\u4ee3\u8868\u8fd9\u6279\u6837\u672c\u771f\u7684\u6a2a\u8de8\u4e86 RY\u2192EO \u5168\u7a0b**\uff1b\u6362\u4e00\u6279\u6837\u672c\uff0c\u540c\u4e00\u4e2a\u6837\u672c\u7684\u5206\u6570\u4f1a\u53d8\u3002",
    "\u672c\u7ed3\u679c unverifiable\uff1a\u6ca1\u6709\u4e00\u6279\u6837\u672c\u540c\u65f6\u6d4b\u8fc7\u8102\u8d28\u5206\u5b50\u79cd\u4e0e\u8f6c\u5f55/\u86cb\u767d\uff0c\u65e0\u6cd5\u8bad\u7ec3\u4e5f\u65e0\u6cd5\u9a8c\u8bc1\u8868\u8fbe\u2192\u4f4d\u7f6e\u7684\u6620\u5c04\uff1b\u56e0\u6b64\u4e0d\u7ed9 p \u503c\u3001\u4e0d\u7ed9\u7f6e\u4fe1\u533a\u95f4\u3001\u4e0d\u7ed9 MAE\u3002",
    "\u7f6e\u4fe1\u5ea6\u662f\u8bc1\u636e\u8986\u76d6\u3001\u5148\u9a8c\u65b9\u5411\u4e00\u81f4\u6027\u3001\u7279\u5f81\u6a21\u578b\u7559\u4e00\u89e3\u91ca\u5ea6\u4e0e\u5206\u6570\u8de8\u5ea6\u7684\u63cf\u8ff0\u6027\u5408\u6210\uff0c**\u4e0d\u662f\u6821\u51c6\u6982\u7387\uff0c\u4e5f\u4e0d\u4ee3\u8868\u9884\u6d4b\u6b63\u786e\u7387**\uff1b\u6863\u4f4d\u6700\u9ad8\u53ea\u5230 low\uff0c\u4e0d\u8bbe moderate/high\u3002",
    "\u8102\u8d28\u5c42\u7684 A \u6863\u7ed3\u679c\uff08\u7559\u4e00\u9636\u6bb5 MAE 0.24\uff09\u4e0e\u672c\u6863\u4e0d\u53ef\u6bd4\u8f83\uff1b\u80fd\u7528\u8102\u8d28\u5206\u5b50\u79cd\u65f6\u4e0d\u8981\u770b\u8fd9\u4e00\u6863\u3002"
  )
  en <- c(
    "Mechanistic boundary: Minetti 2025 (PMID 40016214) report that reticulocytes and erythrocytes cannot synthesise phospholipids de novo; remodelling reflects selective removal and exchange with plasma (possibly via VPS13A). Linking ELOVL/FADS/LPCAT expression to red-cell acyl-chain length is therefore an INDIRECT HYPOTHESIS.",
    "The score is a WITHIN-BATCH relative-rank projection onto the 0-100 scale (the five stages are only ticks), not an absolute storage day or a cross-batch comparable age coordinate. Because it is a rank projection, the extremes necessarily land near 0 and 100; this does NOT mean the batch spans the whole RY-to-EO range.",
    "This output is UNVERIFIABLE: no cohort measures both the 45 lipid species and transcript/protein abundance, so an expression-to-position map can be neither trained nor validated. No p-value, interval or MAE is supplied.",
    "Confidence is a descriptive combination of evidence coverage, prior-direction agreement, leave-one-out fit of the feature model and score spread. It is NOT a calibrated probability; the highest grade is 'low'.",
    "The lipid-only tier-A result (leave-one-stage-out MAE 0.24) is not comparable with this tier; ignore this block whenever species-level lipidomics is available."
  )
  if (abstain) {
    zh <- c(zh, "\u672c\u6b21\u5f03\u6743\uff1a\u547d\u4e2d\u5148\u9a8c\u4e0d\u8db3 / \u6837\u672c\u6570\u4e0d\u8db3 / \u5148\u9a8c\u65b9\u5411\u4e92\u76f8\u62b5\u6d88 / \u6279\u5185\u5206\u6570\u65e0\u5206\u8fa8\u529b\u3002")
    en <- c(en, paste("Abstained: too few informative priors, too few samples,",
                      "cancelling prior directions, or no spread of the batch scores."))
  }
  list(zh = zh, en = en)
}

## 基因/蛋白矩阵 → 先验外推位置 + 描述性置信度
.age_project_prior <- function(x, permutations = NULL) {
  fm <- ery_age_feature_model(permutations)
  P <- .age_prior_table()
  cl <- .clean_matrix(x)
  id <- ery_id_space(cl$features)
  base <- list(ok = FALSE, unverifiable = TRUE, p_value = NULL, ci = NULL,
               accuracy = NULL, basis = "within-batch relative rank",
               diagnostics = fm,
               sources = .age_prior_sources(),
               n_prior_rows = nrow(P))

  if (!identical(id$verdict, "gene")) {
    base$reason_code <- "non_gene_space"
    base$reason <- sprintf("\u7279\u5f81\u540d\u4e0d\u5728\u57fa\u56e0\u7b26\u53f7\u7a7a\u95f4\uff08%s\uff09\uff0c\u8102\u8d28\u8c03\u63a7\u5148\u9a8c\u65e0\u6cd5\u5bf9\u9f50\u3002",
                           id$verdict)
    base$reason_en <- sprintf(paste("Feature identifiers are not in gene-symbol space",
                                    "(%s); lipid-regulator priors cannot be aligned."),
                              id$verdict)
    base$caveats <- .age_prior_caveats(TRUE)
    return(base)
  }
  sp <- ery_to_space(cl$x, cl$features, "gene")
  sym <- sp$x
  rownames(sym) <- sp$features

  lookup <- stats::setNames(seq_len(nrow(sym)), toupper(rownames(sym)))
  hit <- list()
  for (i in seq_len(nrow(P))) {
    lab <- lookup[toupper(P$gene[i])]
    if (is.na(lab)) next
    row <- as.numeric(sym[lab, ])
    if (!any(is.finite(row))) next
    w <- P$tier_w[i] * (fm$directions[[P$feature[i]]] %||% 0) * P$effect[i]
    hit[[length(hit) + 1L]] <- list(
      gene = P$gene[i], process = P$process[i], feature = P$feature[i],
      effect = P$effect[i], tier = P$tier[i], source = P$source[i],
      label = rownames(sym)[lab], values = row, weight = as.numeric(w)
    )
  }
  directional <- Filter(function(h) abs(h$weight) > 1e-9, hit)
  context <- Filter(function(h) abs(h$weight) <= 1e-9, hit)

  n_samples <- ncol(sym)
  n_dir_prior <- sum(P$effect != 0 & P$tier != "C")
  coverage <- length(directional) / max(1, n_dir_prior)
  w <- vapply(directional, `[[`, 0, "weight")
  agreement <- if (length(w)) abs(sum(w)) / sum(abs(w)) else 0

  reason <- NULL
  if (length(directional) < 8L) reason <- "too_few_priors"
  else if (n_samples < 3L) reason <- "too_few_samples"
  else if (agreement < 0.05) reason <- "prior_directions_cancel"

  score <- NULL
  if (is.null(reason)) {
    ## 每个基因在批内取秩百分位；NaN 用该基因中位数补齐，避免整行丢弃
    R <- matrix(NA_real_, nrow = length(directional), ncol = n_samples)
    for (k in seq_along(directional)) {
      v <- directional[[k]]$values
      med <- if (any(is.finite(v))) stats::median(v[is.finite(v)]) else 0
      v[!is.finite(v)] <- med
      R[k, ] <- (.rank_profile(v) + 0.5) / length(v)
    }
    score <- as.numeric(crossprod(w, R - 0.5)) / sum(abs(w))
    if (stats::sd(score) < 1e-6) reason <- "no_score_spread"
  }

  base <- utils::modifyList(base, list(
    n_hit = length(hit), n_directional = length(directional),
    n_context_only = length(context), n_samples = n_samples,
    id_space = list(before = nrow(cl$x), after = nrow(sym), verdict = id$verdict),
    coverage = round(coverage, 4), direction_agreement = round(agreement, 4),
    n_tier_a = sum(vapply(hit, `[[`, "", "tier") == "A"),
    n_tier_b = sum(vapply(hit, `[[`, "", "tier") == "B"),
    feature_loo_r2 = fm$loo_r2
  ))

  if (!is.null(reason)) {
    base$reason_code <- reason
    base$reason <- switch(
      reason,
      too_few_priors = sprintf("\u547d\u4e2d\u7684\u5e26\u65b9\u5411\u5148\u9a8c\u53ea\u6709 %d \u6761\uff08\u9608\u503c 8\uff09\uff0c\u4e0d\u8db3\u4ee5\u7ed9\u51fa\u4e00\u4e2a\u65b9\u5411\u3002",
                               length(directional)),
      too_few_samples = sprintf("\u672c\u6279\u53ea\u6709 %d \u4e2a\u6837\u672c\uff0c\u6279\u5185\u76f8\u5bf9\u987a\u5e8f\u6ca1\u6709\u610f\u4e49\u3002", n_samples),
      prior_directions_cancel = sprintf("\u547d\u4e2d\u5148\u9a8c\u7684\u6b63\u8d1f\u65b9\u5411\u51e0\u4e4e\u5b8c\u5168\u62b5\u6d88\uff08\u4e00\u81f4\u6027 %.3f\uff09\uff0c\u4e0d\u7ed9\u4f4d\u7f6e\u3002",
                                        agreement),
      no_score_spread = "\u6279\u5185\u5206\u6570\u6ca1\u6709\u5206\u8fa8\u529b\uff0c\u65e0\u6cd5\u6392\u5e8f\u3002")
    base$caveats <- .age_prior_caveats(TRUE)
    return(base)
  }

  pct <- (rank(score, ties.method = "average") - 0.5) / n_samples
  conf <- 0.35 * min(1, coverage) + 0.25 * agreement +
    0.25 * max(0, min(1, fm$loo_r2)) + 0.15 * min(1, n_samples / 20)
  ## 最高只到 low：跨层没有配对数据可验证，覆盖再好也不能叫 moderate/high
  grade <- if (coverage < 0.35 || agreement < 0.20 || fm$loo_r2 < 0.10) {
    "exploratory"
  } else "low"

  rows <- lapply(seq_len(n_samples), function(j) {
    pos <- 4 * pct[j]
    list(sample = colnames(sym)[j], position = round(pos, 3),
         senescence_score = round(pos / 4 * 100, 1),
         score = round(score[j], 5), zone = .ery_zone(pos))
  })
  detail <- lapply(directional, function(h) {
    list(gene = h$gene, process = h$process, feature = h$feature,
         effect = h$effect, tier = h$tier, source = h$source,
         weight = round(h$weight, 4))
  })

  utils::modifyList(base, list(
    ok = TRUE, samples = rows, prior_rows = detail,
    confidence = round(conf, 4), grade = grade,
    confidence_parts = list(
      coverage = round(coverage, 4),
      direction_agreement = round(agreement, 4),
      feature_loo_r2 = fm$loo_r2,
      score_sd = round(stats::sd(score), 5),
      n_samples = n_samples),
    confidence_formula = paste("0.35*min(1,coverage) + 0.25*direction_agreement",
                               "+ 0.25*clamp(loo_R2,0,1) + 0.15*min(1,n_samples/20)"),
    caveats = .age_prior_caveats(FALSE)
  ))
}

#' Place samples on the reticulocyte-to-senescent-red-cell maturation axis
#'
#' Answers "how old is this red cell, relative to the maturation trajectory".
#' The primary reading is a 0-100 senescence score derived from a one-dimensional
#' projection onto the membrane-lipid trajectory built from five pooled stages
#' (reticulocyte young to erythrocyte old). Alongside it: a permutation p-value
#' for the profile shape, a per-lipid agreement rate with its own binomial test
#' (does the *whole* profile move together, which matters for compositional
#' data), and the leave-one-stage reconstruction error as the only honest
#' accuracy figure.
#'
#' Identity is checked before anything is projected. A PBMC or plasma sample has
#' no biological meaning on a red-cell lipid axis, so blocked samples receive no
#' score at all rather than a low one.
#'
#' When the input is not lipidomics (no lipid names resolve), the function
#' abstains and — only then — attaches a cross-layer extrapolation from
#' lipid-regulator priors. That tier is unverifiable by construction and is
#' never shown next to a measured result.
#'
#' @param x Numeric matrix, lipids in rows (names such as `PC:34:1 (16:0-18:1w9)`
#'   or `PC(34:1)`), samples in columns. Values must be positive (intensities or
#'   percentages): compositional closure on log-scale values would produce a
#'   ratio with no meaning, so the function abstains instead.
#' @param annotation Optional data.frame used by the identity gate.
#' @param prior_input Optional gene/protein matrix for the cross-layer tier.
#'   Defaults to `x`.
#' @param n_perm Permutations for the shape test (default 200).
#' @param seed Seed for the permutation RNG.
#' @param permutations Optional matrix of 1-based permutations (length = matched
#'   lipids). Supplying it makes the p-value reproducible across implementations.
#' @return A list with `ok`, `samples`, `summary`, `coverage`, `validation`,
#'   `gate`, `knowledge`-free `prior_extrapolation`, `note` and `caveats`.
#' @examples
#' \donttest{
#' m <- read.delim(system.file("extdata", "PMC11868425_lipids.tsv.gz",
#'                             package = "ErythroWorkbench"), row.names = 1)
#' ery_age(m)$summary$median_senescence_score
#' }
#' @export
ery_age <- function(x, annotation = NULL, prior_input = NULL, n_perm = 200L,
                    seed = 20261112L, permutations = NULL) {
  cl <- .clean_matrix(x)
  mat <- cl$x
  rownames(mat) <- cl$features

  gate <- .age_gate(mat, annotation)
  keep <- gate$passed
  if (length(keep) && length(keep) < ncol(mat)) {
    mat <- mat[, colnames(mat) %in% keep, drop = FALSE]
  }
  gate$n_after <- ncol(mat)

  if (!length(keep)) {
    ## 一个都没通过：直接弃权，绝不能让投影拿原矩阵又算出一份分数，
    ## 否则页面会同时显示「全部拦截」和一张结果表，自相矛盾
    gate$n_after <- 0
    return(list(
      ok = FALSE, reason_code = "no_erythroid_samples",
      reason = "\u5168\u90e8\u6837\u672c\u88ab\u7ea2\u7cfb\u8eab\u4efd\u95f8\u95e8\u62e6\u622a\uff1a\u8fd9\u662f\u4e00\u4efd\u975e\u7ea2\u7ec6\u80de\u6837\u672c\uff08\u767d\u7ec6\u80de / \u8840\u6d46 / \u5176\u4ed6\u7ec4\u7ec7 / \u7ec6\u80de\u7cfb\uff09\uff0c\u5e74\u9f84\u8f74\u5728\u5b83\u4eec\u8eab\u4e0a\u6ca1\u6709\u5b9a\u4e49\u3002",
      reason_en = paste("Every sample was blocked by the erythroid identity gate: the input",
                        "is not red cells, and the axis is not defined for them."),
      gate = gate, samples = list(), summary = NULL, coverage = list(),
      prior_extrapolation = NULL,
      caveats = list(
        zh = "\u5e74\u9f84\u8f74\u53ea\u5728\u7ea2\u7ec6\u80de\u819c\u8102\u8f68\u8ff9\u4e0a\u5b9a\u4e49\uff1b\u88ab\u95f8\u95e8\u62e6\u622a\u7684\u6837\u672c\u4e0d\u8f93\u51fa\u4efb\u4f55\u5e74\u9f84\u8bfb\u6570\u3002",
        en = paste("The axis is defined on red-cell membrane lipids; blocked samples",
                   "receive no age score of any kind."))
    ))
  }

  res <- .lipid_project_matrix(mat, permutations)
  res$gate <- gate

  ## 无论 A 档还是 C 档，主读数统一改成 0–100 的「衰老评分」；0–4 位置退到
  ## 诊断 / 下载里当参考轨迹坐标，不再拿 RY/RM/EY/EM/EO 当样本标签
  for (i in seq_along(res$samples)) {
    s <- res$samples[[i]]
    if (!is.null(s$position) && is.numeric(s$position)) {
      s$senescence_score <- round(s$position / 4 * 100, 1)
      res$samples[[i]] <- s
    }
  }
  if (is.list(res$summary) && !is.null(res$summary$median_position)) {
    res$summary$median_senescence_score <- round(res$summary$median_position / 4 * 100, 1)
  }

  ## C 档只在 A 档定位不了时给，避免把不可验证的数字和实测结果并列
  if (!isTRUE(res$ok)) {
    pe <- tryCatch(
      .age_project_prior(if (is.null(prior_input)) x else prior_input),
      error = function(e) list(ok = FALSE, unverifiable = TRUE,
                               reason_code = "prior_engine_error",
                               reason = sprintf("\u5148\u9a8c\u5916\u63a8\u5f15\u64ce\u5f02\u5e38\uff1a%s", conditionMessage(e)),
                               caveats = .age_prior_caveats(TRUE)))
    res$prior_extrapolation <- pe

    code <- res$reason %||% {
      rr <- res$refusal_reasons
      if (length(rr)) names(rr)[which.max(unlist(rr))] else "no_lipid_names"
    }
    res$reason_code <- code
    res$reason <- switch(
      code,
      insufficient_overlap = sprintf("\u547d\u4e2d\u8f68\u8ff9\u5206\u5b50\u79cd\u4e0d\u8db3\uff08\u9700\u8981 \u2265%d \u4e2a\uff0c\u5b9e\u9645 %d \u4e2a\uff09\uff1a\u8fd9\u4efd\u77e9\u9635\u7684\u8102\u8d28\u547d\u540d\u6216\u9762\u677f\u8986\u76d6\u4e0e\u53c2\u8003\u8f68\u8ff9\u6ca1\u6709\u4ea4\u96c6\u3002",
                                     .ery_min_lipids, res$coverage$n_matched %||% 0),
      off_trajectory = "\u8c31\u5f62\u4e0e\u8f68\u8ff9\u76f8\u5173\u8fc7\u4f4e\uff1a\u8fd9\u662f\u4e00\u4e2a\u8102\u8d28\u7ec4\u77e9\u9635\uff0c\u4f46\u6574\u6761\u8c31\u4e0d\u6cbf\u8fd9\u6761\u6210\u719f\u8f68\u8ff9\u79fb\u52a8\uff08\u8de8\u5b9e\u9a8c\u5ba4\u9762\u677f\u5e38\u89c1\uff09\u3002",
      requires_positive_values = "\u77e9\u9635\u542b\u975e\u6b63\u503c\uff08log2 \u5f3a\u5ea6\u6216 log \u6bd4\u503c\uff09\uff1a\u7ec4\u6210\u578b\u95ed\u5408\u9700\u8981\u539f\u59cb\u5f3a\u5ea6\u6216\u767e\u5206\u6bd4\uff0c\u786c\u95ed\u5408\u4f1a\u9020\u51fa\u6ca1\u6709\u610f\u4e49\u7684\u6bd4\u4f8b\u3002",
      insufficient_varying_lipids = sprintf("\u547d\u4e2d\u5206\u5b50\u79cd\u91cc\u53ef\u7528\u7684\u4e0d\u5230 %d \u4e2a\uff08\u9636\u6bb5\u95f4\u65e0\u53d8\u5f02\uff09\uff0c\u65e0\u6cd5\u5b9a\u4f4d\u3002",
                                            .ery_min_lipids),
      no_lipid_names = "\u884c\u540d\u4e0d\u662f\u53ef\u89e3\u6790\u7684\u8102\u8d28\u5206\u5b50\u79cd\uff08\u4f8b\u5982 PC:34:1 / SM(d42:1) / TG(18:2_17:1_18:2)\uff09\uff0c\u65e0\u6cd5\u6295\u5f71\u3002",
      sprintf("\u65e0\u6cd5\u5b9a\u4f4d\uff1a%s", code))
    res$reason_en <- switch(
      code,
      insufficient_overlap = sprintf("Too few trajectory species matched (need >= %d, got %d).",
                                     .ery_min_lipids, res$coverage$n_matched %||% 0),
      off_trajectory = paste("The profile does not move along this maturation trajectory",
                             "(common for cross-laboratory panels)."),
      requires_positive_values = paste("Non-positive values (log2 intensity or log ratio):",
                                       "compositional closure needs raw intensities or percentages."),
      insufficient_varying_lipids = sprintf("Fewer than %d usable species after the variance filter.",
                                            .ery_min_lipids),
      no_lipid_names = "Row names are not parseable lipid species.",
      sprintf("Not placeable: %s", code))
    if (length(res$coverage)) {
      res$note <- sprintf("\u547d\u4e2d %d / %d \u4e2a\u8f68\u8ff9\u5206\u5b50\u79cd\uff08A \u7ea7 %d\u3001B \u7ea7 %d\uff09\u3002",
                          res$coverage$n_matched, res$coverage$n_anchor,
                          res$coverage$tier_A, res$coverage$tier_B)
      res$note_en <- sprintf("%d of %d trajectory species matched (tier A %d, tier B %d).",
                             res$coverage$n_matched, res$coverage$n_anchor,
                             res$coverage$tier_A, res$coverage$tier_B)
    }
    return(res)
  }

  res$prior_extrapolation <- NULL
  summ <- res$summary
  res$note <- sprintf("%d \u4e2a\u6837\u672c\u7ed9\u51fa\u4f4d\u7f6e\uff0c%d \u4e2a\u5f03\u6743\uff1b\u4e2d\u4f4d\u4f4d\u7f6e %.2f\uff08%s\uff09\u3002",
                      summ$n_ok, summ$n_refused, summ$median_position,
                      summ$median_zone$zh)
  res$note_en <- sprintf("%d samples placed, %d refused; median position %.2f (%s).",
                         summ$n_ok, summ$n_refused, summ$median_position,
                         summ$median_zone$en)
  res$caveats <- list(
    zh = c(
      "\u4f4d\u7f6e\u662f**\u5728\u53c2\u8003\u8f68\u8ff9\u4e0a\u7684\u4e00\u7ef4\u6295\u5f71**\uff0c\u4e0d\u662f\u5206\u7c7b\u5668\u8f93\u51fa\u7684\u6982\u7387\u3002\u53c2\u8003\u7814\u7a76\u6bcf\u4e2a\u9636\u6bb5\u53ea\u6709\u4e00\u4e2a\u5408\u5e76\u6837\u672c\uff084 \u540d\u4f9b\u8005\u5408\u5e76\uff09\uff0c\u65e0\u6cd5\u505a\u7559\u4e00\u9636\u6bb5\u7684\u4e2a\u4f53\u7ea7\u9a8c\u8bc1\uff1b\u7559\u4e00\u9636\u6bb5\u91cd\u5efa\u66f2\u7ebf\u7684\u7cbe\u5ea6\u89c1 validation\u3002",
      "\u53ea\u80fd\u8bfb\u300c\u76f8\u5bf9\u8fd9\u6761\u8f68\u8ff9\u504f\u5e74\u8f7b/\u504f\u8870\u8001\u300d\uff0c\u4e0d\u80fd\u8bfb\u6210\u7edd\u5bf9\u7ec6\u80de\u5bff\u547d\uff08\u5929\uff09\u3002\u819c\u8102\u8d28\u4e0d\u662f\u65f6\u949f\uff0c\u5b83\u540c\u65f6\u53d7\u50a8\u5b58\u3001\u6c27\u5316\u548c\u4f9b\u8005\u5dee\u5f02\u5f71\u54cd\u3002",
      "\u5360\u6bd4\u662f\u7ec4\u6210\u578b\u6570\u636e\uff1a\u4e00\u4e2a\u5206\u5b50\u79cd\u5347\u9ad8\u5fc5\u7136\u6324\u5360\u522b\u7684\u5206\u5b50\u79cd\uff0c\u56e0\u6b64\u5355\u770b\u67d0\u4e2a\u8102\u8d28\u7684\u9ad8\u4f4e\u6ca1\u6709\u610f\u4e49\uff0c\u8981\u770b\u6574\u6761\u8c31\u662f\u5426\u4e00\u8d77\u79fb\u52a8\uff08\u9010\u8102\u8d28\u4e00\u81f4\u7387\u5c31\u662f\u8fd9\u4e2a\u68c0\u9a8c\uff09\u3002",
      "A \u7ea7\u547d\u4e2d\uff1d\u7c7b\u522b+\u603b\u78b3+\u603b\u53cc\u952e+\u919a\u952e/\u78b1\u57fa\u5168\u540c\uff1bB \u7ea7\u547d\u4e2d\uff1d\u53ea\u5bf9\u5230\u7c7b\u522b+\u603b\u78b3+\u603b\u53cc\u952e\uff0c\u919a\u8102\u4e0e t/d \u7cfb\u5217\u53ef\u80fd\u88ab\u5e76\u5230\u4e00\u8d77\uff0c\u8bfb\u7684\u65f6\u5019\u770b B \u7ea7\u5360\u6bd4\u3002"),
    en = c(
      "The position is a one-dimensional projection onto the reference trajectory, not a classifier probability. Each reference stage is a single pooled sample (4 donors), so individual-level cross-validation is impossible; see validation for the leave-one-stage reconstruction error.",
      "It reads as younger/older relative to this trajectory, not as an absolute cell lifespan in days. Membrane lipids are not a clock: storage, oxidation and donor variation all move them.",
      "Percentages are compositional: one species rising necessarily takes share from others, so a single lipid level means little; what matters is whether the whole profile shifts together.",
      "Tier A = same class, total carbons, total double bonds, ether bond and base; tier B = class + carbons + double bonds only, so ether and t/d series may be merged.")
  )
  res
}
