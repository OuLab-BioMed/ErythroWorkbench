## purity.R — 红细胞身份与混入：主读数是「与本层库内红细胞参照的秩一致性」。
##
## 对应引擎：app/backend/engines/toolkit_tools.py 的 run_purity() 与
## _whole_rank_evidence()，加上 sample_types.py 的 layer_consistency() /
## headline_score()。
##
## 三条不许改口径的地方：
##   1. 秩向量用**全基因**背景下的秩，不是 47 个 marker 内部的相对秩。
##      只用 marker 内部排序时，纯 RBC 与纯 PBMC 的一致性实测都在 67-68，
##      主读数会失去判别力。
##   2. 主读数封顶在 min(95, 50*(p95+1))，所以永远不会出现 100。
##   3. 没有数据观测就弃权，绝不把库内先验当成本样本的读数报出去。

#' Erythroid identity and admixture
#'
#' Scores how much of a sample's signal is erythroid. The primary reading is
#' the rank agreement between the sample's 47 lineage markers (positions within
#' the whole transcriptome) and the per-layer consensus measured over in-house
#' red-cell samples, rescaled to 0-100, then fused with the in-house prior by
#' confidence and capped at the 95th percentile of in-house self-agreement.
#'
#' This is **not** cell-count purity. Mature red cells are enucleated, so what
#' the matrix measures is erythroid-derived RNA.
#'
#' @param x Numeric matrix or data.frame, features in rows and samples in
#'   columns. Row names are gene symbols. Probe-level input should be collapsed
#'   to symbols first (duplicates are averaged).
#' @param layer Omics layer used to pick the reference consensus: `"transcriptome"`
#'   (default) or `"proteome"`. Passed explicitly, because the choice of
#'   reference frame must not be inferred silently from a handful of row names.
#'   The Chinese layer name is available as `attr(result, "layer_label")`; it is
#'   deliberately not a key, see [ery_layer_reference()].
#' @param level Confidence grade used to weight the observation against the
#'   prior: `"high"` (w = 4), `"medium"` (w = 2), `"low"` (w = 1) or
#'   `"unverifiable"` (w = 0, abstain). Use [ery_confidence()] to derive it
#'   from the data instead of guessing.
#'
#' @return A data.frame with one row per sample and the columns
#'   `sample`, `ery_rank`, `leu_rank`, `delta`, `auc`, `q_leu_rna`,
#'   `rna_identity_score`, `ery_consistency`, `consistency_rho`,
#'   `n_markers`, `identity`, `conflict`, `headline`, `prior`, `ceiling`.
#'   `headline` is the primary reading; `NA` means the module abstained.
#'
#' @examples
#' \dontrun{
#' m <- read.delim("matrix.tsv", row.names = 1)
#' ery_purity(m)
#' }
#' @export
ery_purity <- function(x, layer = c("transcriptome", "proteome"), level = "high") {
  layer <- match.arg(layer)
  refs <- .ewb_read_rds("layer_references.rds")
  trefs <- .ewb_read_rds("toolkit_refs.rds")$purity

  cl <- .clean_matrix(x)
  ## 先归一到标准符号空间：别名、带版本号的 Ensembl、UniProt 外壳都走同一套
  ## canonical 规则，多探针保留离散度最大的一行（不平均）。
  sp <- ery_to_space(cl$x, cl$features, "gene")
  id <- ery_id_space(sp$features)
  if (id$verdict != "gene") {
    stop("Input is not on the gene-symbol space (detected: ", id$verdict,
         "). The panel needs symbols such as HBB / GYPA / SLC4A1.", call. = FALSE)
  }

  lg <- ery_as_log2(sp$x)
  L <- lg$x
  features <- sp$features

  ery_all <- unlist(refs$markers$ery, use.names = FALSE)
  leu_all <- unlist(refs$markers$leu, use.names = FALSE)
  has_ery <- sum(ery_all %in% features)
  has_leu <- sum(leu_all %in% features)
  if (has_ery < 3) {
    stop(sprintf(paste0("Only %d of %d erythroid markers found, and %d of %d ",
                        "leukocyte markers. This is not a sample problem: the ",
                        "matrix carries too little lineage evidence. A targeted ",
                        "panel that measured only metabolic enzymes will do this."),
                 has_ery, length(ery_all), has_leu, length(leu_all)), call. = FALSE)
  }

  ## 逐样本：全基因秩归一（列内），这是唯一与引擎一致的口径
  R <- apply(L, 2, function(v) .rank_profile(v) / max(1, sum(is.finite(v)) - 1))
  if (is.null(dim(R))) R <- matrix(R, nrow = 1)
  rownames(R) <- features

  cons_ref <- .find_layer(refs, layer)
  if (is.null(cons_ref)) {
    stop("No in-house reference for layer '", layer, "'.", call. = FALSE)
  }
  ## 逐元素转数值：unlist() 会吞掉 null，让 47 维向量左移
  consensus <- .num_vec(cons_ref$consensus)
  if (length(consensus) != length(ery_all) + length(leu_all)) {
    stop("Reference consensus length does not match the marker panels.", call. = FALSE)
  }

  ident <- trefs$identity
  cp <- trefs$cell_prior
  k_slope <- as.numeric(ident$k)
  center <- as.numeric(ident$center)
  fit_min <- as.numeric(ident$fit_r2_min %||% 0.5)
  auc_leu_max <- as.numeric(cp$auc_leu_endpoint_max %||% NA_real_)
  n_leu_min <- 8L

  ## AUC 曲线：auc -> q（白细胞直判的把握），按 auc 升序插值
  curve <- cp$auc_curve
  if (length(curve)) {
    ax <- vapply(curve, function(z) as.numeric(z[[2]]), 0)   # auc
    qy <- vapply(curve, function(z) as.numeric(z[[1]]), 0)   # q
    o <- order(ax)
    ax <- ax[o]; qy <- qy[o]
  } else {
    ax <- qy <- NULL
  }

  samples <- colnames(L)
  out <- vector("list", length(samples))
  for (j in seq_along(samples)) {
    rj <- R[, j]
    e_full <- ifelse(ery_all %in% features, rj[match(ery_all, features)], NA_real_)
    l_full <- ifelse(leu_all %in% features, rj[match(leu_all, features)], NA_real_)
    full_vec <- c(e_full, l_full)

    er <- .nanmean0(e_full)
    lr <- .nanmean0(l_full)
    d <- er - if (is.finite(lr)) lr else 0.5
    rna_score <- 1 / (1 + exp(-max(-30, min(30, k_slope * (d - center)))))

    a <- e_full[is.finite(e_full)]
    b <- l_full[is.finite(l_full)]
    auc <- if (length(a) && length(b)) {
      mean(outer(a, b, ">")) + 0.5 * mean(outer(a, b, "=="))
    } else {
      NA_real_
    }
    q <- if (is.finite(auc) && !is.null(ax)) {
      max(0, min(1, stats::approx(ax, qy, xout = auc, rule = 2)$y))
    } else {
      NA_real_
    }

    cons <- .layer_consistency(full_vec, consensus, cons_ref)

    strong_leu <- isTRUE(length(b) >= n_leu_min && is.finite(auc) &&
                          is.finite(auc_leu_max) && auc <= auc_leu_max + 0.005)
    conflict <- isTRUE(is.finite(auc) && auc < 0.5) && !strong_leu

    call <- if (strong_leu) {
      "leukocyte_sample"
    } else if (!is.finite(cons$score)) {
      "unverifiable"
    } else if (cons$score >= 80) {
      "erythroid_sample"
    } else if (cons$score >= 60) {
      "erythroid_dominant"
    } else {
      "mixed"
    }

    h <- .headline_score(cons$score, cons_ref, level)
    out[[j]] <- data.frame(
      sample = samples[j],
      ery_rank = round(er, 3), leu_rank = round(lr, 3), delta = round(d, 3),
      auc = round(auc, 3), q_leu_rna = round(q, 3),
      rna_identity_score = round(rna_score, 3),
      ery_consistency = if (is.null(cons$score)) NA_real_ else cons$score,
      consistency_rho = if (is.null(cons$rho)) NA_real_ else cons$rho,
      n_markers = cons$n_markers,
      identity = call, conflict = conflict,
      headline = h$score, prior = h$prior, ceiling = h$ceiling,
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)

  attr(res, "layer") <- layer
  attr(res, "layer_label") <- cons_ref$key_label
  attr(res, "scale") <- lg$scale
  attr(res, "level") <- level
  attr(res, "n_ref_samples") <- as.integer(cons_ref$n_ref_samples)
  attr(res, "n_ref_datasets") <- as.integer(cons_ref$n_ref_datasets)
  attr(res, "n_ery_markers") <- has_ery
  attr(res, "n_leu_markers") <- has_leu
  attr(res, "reference_provenance") <- ery_provenance()$engine_version
  res
}

## 与 sample_types.layer_consistency() 逐位相同：先 _rank01 两侧，再 Pearson。
.layer_consistency <- function(ery_vec, consensus, cons_ref) {
  ok <- is.finite(ery_vec) & is.finite(consensus)
  if (sum(ok) < 8) {
    return(list(ok = FALSE, score = NULL, rho = NULL, n_markers = sum(ok)))
  }
  x <- .rank01(ery_vec[ok])
  y <- .rank01(consensus[ok])
  if (stats::sd(x) <= 0 || stats::sd(y) <= 0) {
    return(list(ok = FALSE, score = NULL, rho = NULL, n_markers = sum(ok)))
  }
  rho <- max(-1, min(1, mean((x - mean(x)) * (y - mean(y))) /
                       (sqrt(mean((x - mean(x))^2)) * sqrt(mean((y - mean(y))^2)))))
  list(ok = TRUE, score = round(50 * (rho + 1), 1), rho = round(rho, 3),
       n_markers = sum(ok),
       n_ref_samples = cons_ref$n_ref_samples,
       n_ref_datasets = cons_ref$n_ref_datasets)
}

## 与 sample_types.headline_score() 逐位相同。
##   prior   = 50*(p50+1)，同层库内样本彼此一致性的中位数（实算）
##   ceiling = min(95, 50*(p95+1))，所以不会出现 100
##   w       = 置信度给权重，观测与先验按 w:1 融合
## 观测缺失或 unverifiable 时**弃权**（score = NA），不给一个纯先验的数。
.headline_score <- function(observed, cons_ref, level) {
  p50 <- as.numeric(cons_ref$self_rho_p50 %||% 0.5)
  p95 <- as.numeric(cons_ref$self_rho_p95 %||% 0.9)
  prior <- round(50 * (p50 + 1), 1)
  ceiling <- round(min(95, 50 * (p95 + 1)), 1)
  w <- switch(level, high = 4, medium = 2, low = 1, unverifiable = 0, 1)
  score <- if (is.null(observed) || !is.finite(observed) || w == 0) {
    NA_real_
  } else {
    round(min((w * observed + prior) / (w + 1), ceiling), 1)
  }
  list(score = score, prior = prior, ceiling = ceiling, weight = w)
}

#' Confidence grade for the erythroid identity reading
#'
#' The grade is set only by measured quantities: how many in-house samples sit
#' in the same layer, how many leukocyte markers the input actually carries,
#' and whether the two endpoints disagree. It never comes from the sample's
#' declared type alone.
#'
#' @param n_ref_samples Number of in-house red-cell reference samples in the
#'   same layer.
#' @param n_leu_markers Leukocyte markers found in the input.
#' @param conflict Whether the endpoint evidence conflicts.
#' @return `"high"`, `"medium"`, `"low"` or `"unverifiable"`.
#' @export
ery_confidence <- function(n_ref_samples = 0, n_leu_markers = 0, conflict = FALSE) {
  lv <- if (n_ref_samples >= 100) {
    "high"
  } else if (n_ref_samples >= 50) {
    "medium"
  } else if (n_ref_samples >= 10) {
    "low"
  } else {
    "unverifiable"
  }
  order_lv <- c("unverifiable", "low", "medium", "high")
  i <- match(lv, order_lv)
  if (isTRUE(n_leu_markers < 3)) i <- max(1L, i - 1L)
  if (isTRUE(conflict)) i <- max(1L, i - 1L)
  order_lv[i]
}

`%||%` <- function(a, b) if (is.null(a)) b else a
