## ---------------------------------------------------------------------------
## 条件谱匹配：这个样本的谱最接近库里**哪一类条件**
##
## 移植自 app/backend/engines/toolkit_tools.py 的 run_attribution()。三件事按
## 顺序回答，顺序不能颠倒：
##
##   1. 它是不是红细胞（身份 / 细胞类型参照）——不是的话，后面两个问题都错位；
##   2. 像哪一类条件（疾病 / 表型），按**该类里最像的那个队列**排序；
##   3. 这个相似度是不是特征身份带来的（打乱查询侧秩序求零分布）。
##
## 排序口径是最容易被改错的一处：用「条件内平均 ρ」会把同一疾病下分化阶段
## 不同的组一起平均掉。GSE301441 的健康样本与自己组的 ρ 是 0.97，但被 0.63 的
## 终末组拉到均值 0.61，反而输给 MM 的 0.83——健康对照被判成多发性骨髓瘤就是
## 这么来的。所以主排序用 best_rho，mean_rho 只作跨队列一致性参考。
## ---------------------------------------------------------------------------

## 引擎的 MIN_SHARED：共享特征少于这个数就不做相关（30 个名字以下的相关不稳）
.ATTR_MIN_SHARED <- 30L

.attribution_centroids <- function() {
  tr <- .ewb_read_rds("toolkit_refs.rds")
  cent <- tr$attribution$centroids
  if (is.null(cent) || !length(cent)) {
    stop("Reference bundle has no attribution centroids.", call. = FALSE)
  }
  cent
}

## 把 centroid 列表按查询侧的 id 空间对齐：只保留共享特征够多的那些
.attr_align <- function(pool, qpos) {
  out <- list()
  for (c in pool) {
    feats <- as.character(unlist(c$features, use.names = FALSE))
    prof <- .num_vec(c$profile)
    qi <- qpos[feats]
    keep <- !is.na(qi)
    if (sum(keep) < .ATTR_MIN_SHARED) next
    out[[length(out) + 1L]] <- list(
      dataset_id = as.character(c$dataset_id),
      group = as.character(c$group %||% ""),
      disease = as.character(c$disease %||% ""),
      layer = as.character(c$layer %||% ""),
      n = as.integer(c$n %||% NA_integer_),
      role = as.character(c$role %||% "erythroid"),
      shared = sum(keep),
      qi = as.integer(qi[keep]),
      cv = prof[keep]
    )
  }
  out
}

## 按行求秩（并列取平均秩），与 Python 的 _rank_rows 等价
.rank_rows <- function(A) {
  if (!length(A)) return(A)
  t(apply(A, 1L, function(v) rank(v, ties.method = "average") - 1))
}

## 置换零分布：打乱查询侧秩序，看同一条件内最强的相关还能到多高。
## 统计量仍是「同一条件内取最强队列」，非有限值按 -Inf 参与 max。
##
## 随机数：Python 用 numpy 的 PCG64，R 用 Mersenne-Twister，两边的置换序列
## 不可能相同。因此 `permutations` 可以外部传入——测试时喂同一组置换，p 值
## 才能逐位对拍；日常使用留 NULL，各用各的 RNG，n_perm 足够时 p 的抽样误差
## 是 1/(n_perm+1) 量级。
.null_max_rho <- function(q, members, n_perm, permutations = NULL) {
  n_feat <- length(q)
  if (!all(is.finite(q))) return(NULL)
  rb_list <- list()
  for (m in members) {
    cv <- m$cv
    if (length(cv) < 3L || !all(is.finite(cv))) return(NULL)
    rb <- .rank_profile(cv)
    rb <- rb - mean(rb)
    nb <- sqrt(sum(rb * rb))
    if (!is.finite(nb) || nb <= 0) return(NULL)
    rb_list[[length(rb_list) + 1L]] <- list(rb = rb, nb = nb)
  }
  if (!length(rb_list)) return(NULL)

  perms <- permutations
  if (is.null(perms)) {
    perms <- matrix(0L, n_perm, n_feat)
    for (t in seq_len(n_perm)) perms[t, ] <- sample.int(n_feat)
  }
  perms <- as.matrix(perms)
  Qp <- matrix(q[perms], nrow = nrow(perms))

  fat <- matrix(NA_real_, nrow = nrow(perms), ncol = length(members))
  for (k in seq_along(members)) {
    A <- Qp[, members[[k]]$qi, drop = FALSE]
    ra <- .rank_rows(A)
    ra <- ra - rowMeans(ra)
    na <- sqrt(rowSums(ra * ra))
    den <- na * rb_list[[k]]$nb
    rho <- rep(-Inf, nrow(ra))
    okd <- is.finite(den) & den > 0
    if (any(okd)) {
      rho[okd] <- as.numeric((ra[okd, , drop = FALSE] %*% rb_list[[k]]$rb) / den[okd])
    }
    rho[!is.finite(rho)] <- -Inf
    fat[, k] <- rho
  }
  best <- apply(fat, 1L, max)
  best[is.finite(best)]
}

#' Match a sample's profile against the reference condition library
#'
#' Answers, in this order: is the sample erythroid at all, which condition class
#' does it most resemble, and is that similarity carried by feature identity
#' rather than chance. The order matters: running a condition match on a
#' leukocyte or plasma sample is the wrong question, so the identity verdict is
#' reported first and a warning is attached when the sample is not erythroid.
#'
#' Conditions are ranked by the **strongest cohort inside the class** (`best_rho`),
#' not by the class mean. Averaging within a condition mixes differentiation
#' stages: a healthy sample scoring rho = 0.97 against its own group can fall
#' below another class once a poorly matching terminal group is averaged in.
#' `mean_rho` is still reported, as a consistency check across cohorts.
#'
#' @param x Numeric matrix, features in rows and samples in columns. Duplicate
#'   features are collapsed onto the standard identifier space first.
#' @param annotation Optional data.frame used for the annotation prior (a known
#'   label, reported separately from the predicted condition) and for identity.
#' @param top Number of conditions kept per sample (default 8, as on the site).
#' @param n_perm Permutations for the null distribution of the top condition.
#' @param seed Seed for the permutation RNG. Ignored when `permutations` is given.
#' @param permutations Optional matrix (or list of matrices, one per sample) of
#'   1-based permutations used for the null distribution. Supplying them makes
#'   the p-value reproducible across implementations, which is how the parity
#'   tests compare this function with the Python engine.
#' @param quiet Suppress the per-sample progress messages.
#' @return A list with `space`, `n_samples`, `n_centroids`, `n_conditions`,
#'   `results` (one entry per sample) and `caveats`.
#' @examples
#' \donttest{
#' m <- read.delim(system.file("extdata", "GSE63703_matrix.tsv.gz",
#'                             package = "ErythroWorkbench"), row.names = 1)
#' a <- ery_attribution(m)
#' a$results[[1]]$top$label
#' }
#' @export
ery_attribution <- function(x, annotation = NULL, top = 8L, n_perm = 200L,
                            seed = 20261110L, permutations = NULL,
                            quiet = TRUE) {
  cl <- .clean_matrix(x)
  id <- ery_id_space(cl$features)
  space <- switch(id$verdict,
                  gene = "gene",
                  metabolite = "metabolite",
                  mirna = "mirna",
                  NA_character_)
  if (is.na(space)) {
    stop("Feature identifiers are not on a supported space (detected: ",
         id$verdict, ").", call. = FALSE)
  }
  sp <- ery_to_space(cl$x, cl$features, space)
  X <- sp$x
  names_ <- sp$features

  ## 与引擎 _prepare 一致：缺失用行中位数填补（秩归一后全有限，除数才是 n-1）
  finite <- is.finite(X)
  for (i in seq_len(nrow(X))) {
    row <- X[i, ]
    if (any(!is.finite(row))) {
      fill <- stats::median(row, na.rm = TRUE)
      row[!is.finite(row)] <- if (is.finite(fill)) fill else 0
      X[i, ] <- row
    }
  }
  Q <- apply(X, 2L, .rank01)
  if (is.null(dim(Q))) Q <- matrix(Q, nrow = nrow(X))
  qpos <- stats::setNames(seq_along(names_), names_)

  cent <- .attribution_centroids()
  cent <- Filter(function(c) identical(as.character(c$space), space), cent)
  if (!length(cent)) {
    stop("No reference profile in the ", space, " space.", call. = FALSE)
  }
  ref_cent <- Filter(function(c) identical(as.character(c$role %||% ""),
                                           "celltype_reference"), cent)
  ery_cent <- Filter(function(c) identical(as.character(c$role %||% "erythroid"),
                                           "erythroid"), cent)
  if (!length(ery_cent)) {
    stop("No erythroid reference profile; condition matching needs one.",
         call. = FALSE)
  }
  aligned <- .attr_align(ery_cent, qpos)
  aligned_ref <- .attr_align(ref_cent, qpos)

  ## 查询样本自己的身份：注释能认就用注释，认不出看样本名
  ident_tab <- NULL
  cols <- colnames(cl$x)
  if (!is.null(cols)) {
    ident_tab <- tryCatch(ery_identify(cl$x, annotation), error = function(e) NULL)
  }
  ident_of <- function(sample) {
    if (is.null(ident_tab)) return(NULL)
    i <- match(sample, ident_tab$sample)
    if (is.na(i)) return(NULL)
    list(sample = ident_tab$sample[i], compartment = ident_tab$compartment[i],
         in_scope = ident_tab$in_scope[i], confidence = ident_tab$confidence[i],
         evidence = if (nzchar(ident_tab$evidence[i]))
           strsplit(ident_tab$evidence[i], ";", fixed = TRUE)[[1]] else character(0))
  }

  ann_rows <- NULL
  if (!is.null(annotation) && is.data.frame(annotation) && nrow(annotation)) {
    al <- .si_row_for(cols, annotation)
    ann_rows <- al$rows
  }
  ann_prior_of <- function(sample, i) {
    if (is.null(ann_rows)) return(NULL)
    j <- ann_rows[i]
    if (is.na(j)) return(NULL)
    row <- annotation[j, , drop = FALSE]
    for (pc in c("disease", "diagnosis")) {
      k <- which(tolower(trimws(names(row))) == pc)
      if (!length(k)) next
      v <- trimws(as.character(row[[k[1]]])[1])
      if (nzchar(v) && !tolower(v) %in% c("nan", "none", "unknown", "na")) return(v)
      return(NULL)
    }
    NULL
  }

  perms_for <- function(j) {
    if (is.null(permutations)) return(NULL)
    if (!is.list(permutations)) return(permutations)
    ## 列表比样本数短时不能直接 [[j]]：越界会抛 "subscript out of bounds"，
    ## 而语义上"这个样本没有指定置换"就是 NULL（用调用方的 RNG）
    if (j > length(permutations)) return(NULL)
    permutations[[j]]
  }

  results <- vector("list", length(cols))
  for (j in seq_along(cols)) {
    sample <- cols[j]
    q <- Q[, j]
    ident <- ident_of(sample)

    hits <- list()
    for (a in aligned) {
      rho <- .spearman(q[a$qi], a$cv)
      if (is.finite(rho)) {
        hits[[length(hits) + 1L]] <- c(a[c("dataset_id", "group", "disease",
                                           "layer", "n", "role", "shared")],
                                       list(rho = round(rho, 3), rho_raw = rho, a = a))
      }
    }
    if (!length(hits)) {
      results[[j]] <- list(sample = sample, ok = FALSE,
                           reason = "no_centroid_with_enough_shared_features",
                           space = space, identity = ident)
      next
    }

    ## 非红系参照里最像的那一个（回答「像不像红细胞」）
    ref_best <- NULL
    for (a in aligned_ref) {
      rho <- .spearman(q[a$qi], a$cv)
      if (is.finite(rho) && (is.null(ref_best) || rho > ref_best$rho)) {
        ref_best <- list(dataset_id = a$dataset_id, group = a$group,
                         disease = a$disease, role = a$role, rho = round(rho, 3))
      }
    }

    ## 条件级聚合：label 用 disease 字段
    labels <- vapply(hits, `[[`, "", "disease")
    ## 按**首次出现顺序**分组，与 Python dict 的插入序一致。默认 split 按字母序，
    ## best_rho 并列时会让两边的「第一名」不同——这类差异只在少数样本上出现，
    ## 排查成本极高，所以在源头对齐。
    by_cond <- split(seq_along(hits), factor(labels, levels = unique(labels)))
    conds <- list()
    for (label in names(by_cond)) {
      hs <- hits[by_cond[[label]]]
      rhos <- vapply(hs, `[[`, 0, "rho")
      rhos_raw <- vapply(hs, `[[`, 0, "rho_raw")
      best <- hs[[which.max(rhos)]]
      ord <- order(-vapply(hs, `[[`, 0, "rho"))
      conds[[length(conds) + 1L]] <- list(
        label = label,
        n_cohorts = length(hs),
        mean_rho = round(mean(rhos), 3),
        best_rho = round(max(rhos), 3),
        mean_rho_raw = mean(rhos_raw), best_rho_raw = max(rhos_raw),
        best_dataset = best$dataset_id,
        best_group = best$group,
        layers = sort(unique(vapply(hs, `[[`, "", "layer"))),
        cohorts = lapply(hs[ord], function(h) {
          list(dataset_id = h$dataset_id, group = h$group, rho = h$rho,
               n = h$n, shared = h$shared)
        })
      )
    }
    ## 排序看条件内最强队列，不看平均（见文件头注释里的 GSE301441 反例）
    o <- order(-vapply(conds, `[[`, 0, "best_rho"),
               -vapply(conds, `[[`, 0, "mean_rho"))
    conds <- conds[o]
    topc <- conds[[1]]

    members <- lapply(hits[by_cond[[topc$label]]], `[[`, "a")
    null <- .null_max_rho(q, members, n_perm, perms_for(j))
    if (is.null(null)) {
      ## 有非有限值 / 秩退化：退回逐次实现，语义相同只是慢
      null <- numeric(0)
      nf <- length(q)
      for (t in seq_len(n_perm)) {
        pm <- perms_for(j)
        idx <- if (is.null(pm)) sample.int(nf) else pm[t, ]
        qq <- q[idx]
        vals <- vapply(members, function(m) .spearman(qq[m$qi], m$cv), 0)
        vals <- vals[is.finite(vals)]
        if (length(vals)) null <- c(null, max(vals))
      }
    }
    p <- if (length(null)) {
      (1 + sum(null >= topc$best_rho)) / (1 + length(null))
    } else NULL
    bg <- stats::median(vapply(hits, `[[`, 0, "rho"))

    prior <- ann_prior_of(sample, j)
    agree <- NULL
    if (!is.null(prior) && nzchar(prior)) {
      agree <- grepl(prior, topc$label, fixed = TRUE) ||
        grepl(topc$label, prior, fixed = TRUE)
    }

    ## 与「非红系参照」比与任何红系条件都更像 → 它根本不是红细胞，
    ## 这条必须排在条件名前面
    nonery_dominant <- !is.null(ref_best) &&
      identical(ref_best$role, "celltype_reference") &&
      ref_best$rho >= topc$best_rho

    warning <- NULL
    if (!is.null(ident) && ident$compartment %in% c("leukocyte", "plasma", "other_tissue")) {
      ev <- if (length(ident$evidence)) paste(ident$evidence, collapse = "\uff1b") else "\u6309\u540d\u79f0"
      warning <- sprintf(
        "\u6837\u672c\u8eab\u4efd\u5224\u5b9a\u4e3a\u300c%s\u300d\uff08%s\uff09\uff0c\u4e0d\u5c5e\u4e8e\u7ea2\u7cfb\u6837\u672c\uff1b\u5bf9\u5b83\u505a\u6761\u4ef6\u5f52\u56e0\u662f\u9519\u4f4d\u63d0\u95ee\uff0c\u4e0b\u9762\u7684\u6761\u4ef6\u540d\u53ea\u80fd\u5f53\u4f5c\u300c\u5b83\u6700\u4e0d\u50cf\u7ea2\u7ec6\u80de\u300d\u7684\u65c1\u8bc1\u3002",
        ident$compartment, ev)
    } else if (nonery_dominant) {
      warning <- sprintf(
        "\u8be5\u6837\u672c\u4e0e\u7ec6\u80de\u7c7b\u578b\u53c2\u7167\uff08%s / %s\uff0c\u03c1=%s\uff09\u6bd4\u4efb\u4f55\u7ea2\u7cfb\u6761\u4ef6\u90fd\u66f4\u50cf\uff0c\u8bf4\u660e\u5b83\u66f4\u50cf\u975e\u7ea2\u7cfb\u7ec6\u80de\u800c\u4e0d\u662f\u7ea2\u7ec6\u80de\uff1b\u6761\u4ef6\u5f52\u56e0\u7684\u7ed3\u8bba\u4e0d\u53ef\u7528\u3002",
        ref_best$dataset_id, ref_best$group, ref_best$rho)
    }

    consistency <- if (topc$best_rho != 0) round(topc$mean_rho / topc$best_rho, 3) else NULL
    reading <- sprintf(
      "\u6700\u63a5\u8fd1\u300c%s\u300d\uff1a\u8be5\u7c7b\u91cc\u6700\u50cf\u7684\u961f\u5217 %s%s\uff0c\u03c1=%s\uff08\u5168\u90e8\u961f\u5217\u4e2d\u4f4d \u03c1=%s\uff0c\u9ad8\u51fa %s\uff09\uff1b\u7b2c\u4e8c\u540d\u300c%s\u300d\u03c1=%s\u3002\u8de8\u961f\u5217\u4e00\u81f4\u6027 %s\uff08\u5e73\u5747 \u03c1=%s\uff0c%d \u4e2a\u961f\u5217\uff09\u3002",
      topc$label,
      topc$best_dataset,
      if (nzchar(topc$best_group)) paste0(" / ", topc$best_group) else "",
      topc$best_rho, round(bg, 3), round(topc$best_rho - bg, 3),
      if (length(conds) > 1) conds[[2]]$label else "\u2014",
      if (length(conds) > 1) conds[[2]]$best_rho else "\u2014",
      if (!is.null(consistency)) round(topc$mean_rho / topc$best_rho, 2) else "\u2014",
      topc$mean_rho, topc$n_cohorts)

    results[[j]] <- list(
      sample = sample, ok = TRUE, space = space,
      identity = ident, identity_warning = warning,
      annotation_prior = prior, prior_agrees = agree,
      nonery_reference = ref_best, nonery_dominant = nonery_dominant,
      conditions = utils::head(conds, top),
      ## conditions 只给前 8 条；总数单独记，否则「有几类条件被比过」和
      ## 「返回了几条」会被当成同一个数
      n_conditions = length(conds),
      ## conditions 只给前 8 条；总数单独记，否则「有几类条件被比过」和
      ## 「返回了几条」会被当成同一个数
      n_conditions = length(conds),
      top = list(
        label = topc$label, mean_rho = topc$mean_rho, best_rho = topc$best_rho,
        mean_rho_raw = topc$mean_rho_raw, best_rho_raw = topc$best_rho_raw,
        n_cohorts = topc$n_cohorts,
        agg = "best_rho\uff08\u6761\u4ef6\u5185\u6700\u5f3a\u961f\u5217\uff1bmean_rho \u4f5c\u8de8\u961f\u5217\u4e00\u81f4\u6027\u53c2\u8003\uff09",
        consistency = consistency,
        p = if (!is.null(p)) round(p, 4) else NULL,
        null_q95 = if (length(null)) round(as.numeric(stats::quantile(null, 0.95, names = FALSE)), 3) else NULL,
        margin = round(topc$best_rho - bg, 3)
      ),
      background_median_rho = round(bg, 3),
      background_median_rho_raw = bg,
      runner_up = if (length(conds) > 1) {
        list(label = conds[[2]]$label, best_rho = conds[[2]]$best_rho)
      } else NULL,
      n_cohorts_compared = length(hits),
      reading = reading
    )
    if (!quiet) {
      message(sprintf("  %-16s -> %s (\u03c1=%s)", sample, topc$label, topc$best_rho))
    }
  }

  n_excluded <- length(cent) - length(ery_cent)
  ## 未舍入的主读数单独挂在 attr 上：对外仍给 3 位小数，但跨语言对拍要的是
  ## 原始浮点——舍入值在 .0005 边界上两边会差一个末位，那是噪声不是差异。
  raw <- data.frame(
    sample = vapply(results, function(r) as.character(r$sample), ""),
    best_rho_raw = vapply(results, function(r) {
      if (isTRUE(r$ok)) r$top$best_rho_raw else NA_real_ }, 0),
    mean_rho_raw = vapply(results, function(r) {
      if (isTRUE(r$ok)) r$top$mean_rho_raw else NA_real_ }, 0),
    bg_raw = vapply(results, function(r) {
      if (isTRUE(r$ok)) r$background_median_rho_raw else NA_real_ }, 0),
    stringsAsFactors = FALSE)
  out <- list(
    space = space, n_samples = length(results),
    n_centroids = length(aligned),
    n_reference_centroids = length(aligned_ref),
    n_centroids_excluded = n_excluded,
    n_conditions = length(unique(vapply(ery_cent, function(c) as.character(c$disease), ""))),
    n_perm = as.integer(n_perm),
    results = utils::head(results, 50),
    caveats = list(
      zh = c(
        "\u8f93\u51fa\u662f\u300c\u6700\u50cf\u54ea\u4e00\u7c7b\u6761\u4ef6\u300d\uff0c\u4e0d\u662f\u8bca\u65ad\uff1a\u4efb\u4f55\u6837\u672c\u90fd\u80fd\u627e\u5230\u6700\u50cf\u7684\u4e00\u7c7b\u3002",
        sprintf("\u53c2\u7167\u753b\u50cf\u5148\u6309**\u6837\u672c\u8eab\u4efd**\u7b5b\u8fc7\uff1aPBMC / \u8840\u6d46 / \u5176\u4ed6\u7ec4\u7ec7\u7684\u6837\u672c\u4e0d\u8fdb\u7ea2\u7cfb\u6761\u4ef6\u753b\u50cf\uff08\u672c\u6b21\u6392\u9664 %d \u5f20\uff09\uff0c\u5b83\u4eec\u53ea\u5355\u72ec\u7528\u4f5c\u300c\u50cf\u4e0d\u50cf\u7ea2\u7ec6\u80de\u300d\u7684\u53c2\u7167\u3002", n_excluded),
        "p \u503c\u53ea\u56de\u7b54\u300c\u8fd9\u4e2a\u76f8\u4f3c\u5ea6\u662f\u4e0d\u662f\u7279\u5f81\u8eab\u4efd\u5e26\u6765\u7684\u300d\uff0c\u4e0d\u56de\u7b54\u300c\u5c5e\u4e8e\u54ea\u4e00\u7c7b\u300d\u3002",
        "\u53ea\u6bd4\u8f83\u540c\u4e00\u6807\u8bc6\u7a7a\u95f4\u7684\u753b\u50cf\uff1b\u8de8\u7a7a\u95f4\uff08\u57fa\u56e0 vs \u4ee3\u8c22\u7269\uff09\u4e0d\u505a\u5916\u63a8\u3002",
        "\u540c\u4e00\u6761\u4ef6\u4e0b\u591a\u4e2a\u961f\u5217\u7684 \u03c1 \u53d6\u7b49\u6743\u5747\u503c\uff1b\u961f\u5217\u6570\u5c11\uff081\u20132 \u4e2a\uff09\u7684\u6761\u4ef6\uff0c\u5747\u503c\u5c31\u662f\u5355\u961f\u5217\u7684\u6570\uff0c\u6ce2\u52a8\u5f88\u5927\u3002",
        "\u8de8\u5e73\u53f0\u53ea\u6bd4\u79e9\u76f8\u5173\uff0c\u5e73\u53f0\u5dee\u5f02\u88ab\u5438\u6536\uff0c\u771f\u5b9e\u7684\u6279\u6b21\u6548\u5e94\u4e5f\u4e00\u8d77\u88ab\u5438\u6536\u6389\u3002"),
      en = c(
        "The output names the closest condition, not a diagnosis: every sample has a closest condition.",
        sprintf("Reference profiles are filtered by sample identity first: PBMC / plasma / other-tissue samples never enter erythroid condition profiles (%d excluded); they serve only as the yes-or-no erythroid reference.", n_excluded),
        "The p-value only covers whether the similarity is carried by feature identity, not which class the sample belongs to.",
        "Only centroids in the same id space are compared; cross-space (gene vs metabolite) is never extrapolated.",
        "Cohorts under one condition are averaged with equal weight; conditions with 1-2 cohorts are single-cohort numbers and move a lot.",
        "Cross-platform comparison uses ranks only, which absorbs platform effects but also real batch effects.")
    )
  )
  attr(out, "raw") <- raw
  out
}
