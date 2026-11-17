## projection.R — 把样本投到固定参考图谱做身份标签转移。
##
## 对应引擎：app/backend/engines/reference_projection.py。
##
## 两条硬约束（引擎 CAVEATS 里写死，这里照搬）：
##   1. 这是**投影不是分类**：几何上任何样本都会落到最近的群，所以 OOD 与
##      投票歧义两道闸门必须同时过，任一不过就弃权、不给标签。
##   2. 验证是留一**数据集**，不是留一样本：同批样本的相似度会被当成精度。
##
## 白细胞参照只有 3 例且全在同一数据集，留掉该数据集后没有白细胞参照，
## 因此对白细胞制品的识别**未经验证**。这个边界必须随结果一起返回。

#' Transfer reference identity by projection onto a fixed atlas
#'
#' Projects each sample onto a fixed atlas of in-house profiled samples and
#' transfers the identity of its nearest neighbours. Two gates decide whether a
#' label is allowed: an out-of-distribution gate (the best similarity must
#' clear the 5th percentile of the atlas's own null distribution) and a margin
#' gate (the winning compartment must beat the runner-up by more than the 10th
#' percentile of margins observed when each dataset is held out).
#'
#' This is a projection, not a classifier: geometrically every sample lands next
#' to some reference group, so a sample that passes neither gate gets **no
#' label**, not a low-confidence one.
#'
#' @param x Numeric matrix or data.frame, features in rows and samples in
#'   columns; row names are gene symbols.
#' @param k Number of neighbours used in the weighted vote. Defaults to the
#'   value baked into the atlas (7) and should not normally be changed.
#'
#' @return A data.frame with one row per sample: `sample`, `label` (`NA` when
#'   abstained), `probability`, `best_rho`, `entropy`, `n_markers`, `abstain`
#'   (reason), and the top-5 neighbours condensed into `top_dataset` /
#'   `top_condition`. Attributes carry the gates and the validation caveats.
#'
#' @examples
#' \dontrun{
#' ery_transfer(read.delim("matrix.tsv", row.names = 1))
#' }
#' @export
ery_transfer <- function(x, k = NULL) {
  atlas <- .ewb_read_rds("reference_atlas.rds")
  cl <- .clean_matrix(x)
  sp <- ery_to_space(cl$x, cl$features, "gene")
  id <- ery_id_space(sp$features)
  if (id$verdict != "gene") {
    stop("Input is not on the gene-symbol space (detected: ", id$verdict, ").",
         call. = FALSE)
  }

  K <- if (is.null(k)) as.integer(atlas$k %||% 7L) else as.integer(k)
  gate <- as.numeric(atlas$ood_gate$rank_p5)
  mgate <- as.numeric(atlas$margin_gate_rank)

  lg <- ery_as_log2(sp$x)
  L <- lg$x
  ## 与引擎一致：秩取**全基因**背景下算，再抽 47 个 marker
  R <- apply(L, 2, function(v) .rank_profile(v) / max(1, sum(is.finite(v)) - 1))

  markers <- atlas$markers
  pos <- match(markers, sp$features)
  ref <- atlas$rank
  meta <- atlas$meta
  lab_idx <- which(meta$labeled)

  samples <- colnames(L)
  out <- vector("list", length(samples))
  for (j in seq_along(samples)) {
    q <- ifelse(!is.na(pos), R[pos, j], NA_real_)
    ## 与全部带标签参考逐对求相似（共同子集上重排秩）
    sims <- rep(NA_real_, nrow(ref))
    for (i in lab_idx) sims[i] <- .sim_rank(q, ref[i, ])

    r <- .vote(sims, lab_idx, meta, gate, mgate, K)
    nb <- lab_idx[order(-ifelse(is.finite(sims[lab_idx]), sims[lab_idx], -9))][1:5]
    nb <- nb[is.finite(sims[nb])]

    out[[j]] <- data.frame(
      sample = samples[j],
      label = r$label, probability = r$prob, best_rho = r$best,
      entropy = r$entropy, n_markers = sum(is.finite(q)),
      abstain = r$why,
      top_dataset = paste(meta$dataset_id[nb], collapse = ";"),
      top_condition = paste(meta$condition[nb], collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  res <- do.call(rbind, out)

  ok <- !is.na(res$label)
  attr(res, "k") <- K
  attr(res, "ood_gate") <- gate
  attr(res, "margin_gate") <- mgate
  attr(res, "n_labeled") <- length(lab_idx)
  attr(res, "n_datasets") <- atlas$n_datasets
  attr(res, "n_answered") <- sum(ok)
  attr(res, "coverage") <- if (nrow(res)) round(mean(ok), 3) else NA_real_
  attr(res, "validation") <- .projection_validation()
  attr(res, "caveats") <- .projection_caveats()
  res
}

## 与 reference_projection._sim 逐位相同：共同非缺失子集上重排秩，再 Pearson。
.sim_rank <- function(a, b, min_common = 8L) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_common) return(NA_real_)
  x <- .rank_profile(a[ok])
  y <- .rank_profile(b[ok])
  sdx <- sqrt(mean((x - mean(x))^2))
  sdy <- sqrt(mean((y - mean(y))^2))
  if (!is.finite(sdx) || !is.finite(sdy) || sdx <= 0 || sdy <= 0) return(NA_real_)
  mean((x - mean(x)) * (y - mean(y))) / (sdx * sdy)
}

## 与 reference_projection._vote 逐位相同。
.vote <- function(sims, ref_idx, meta, gate, mgate, K) {
  ok <- ref_idx[is.finite(sims[ref_idx])]
  if (!length(ok)) {
    return(list(label = NA_character_, prob = NA_real_, best = NA_real_,
                entropy = NA_real_, why = "no_overlap"))
  }
  best <- sims[ok[which.max(sims[ok])]]
  if (best < gate) {
    return(list(label = NA_character_, prob = NA_real_, best = best,
                entropy = NA_real_, why = "ood"))
  }
  top <- ok[order(-sims[ok])][seq_len(min(K, length(ok)))]
  w <- pmax(sims[top], 0)
  if (sum(w) <= 0) {
    return(list(label = NA_character_, prob = NA_real_, best = best,
                entropy = NA_real_, why = "no_positive_similarity"))
  }
  tally <- tapply(w, meta$compartment[top], sum)
  ranked <- sort(as.numeric(tally), decreasing = TRUE)
  lab <- names(tally)[order(-as.numeric(tally))][1]
  s1 <- ranked[1]
  s2 <- if (length(ranked) > 1) ranked[2] else 0
  margin <- (s1 - s2) / sum(w)
  prob <- s1 / sum(w)
  ent <- if (length(tally) > 1) {
    p <- as.numeric(tally) / sum(w)
    -sum(p * log(pmax(p, 1e-12))) / log(length(tally))
  } else 0
  if (margin < mgate) {
    return(list(label = NA_character_, prob = NA_real_, best = best,
                entropy = round(ent, 3), why = "ambiguous"))
  }
  list(label = lab, prob = round(prob, 3), best = round(best, 3),
       entropy = round(ent, 3), why = NA_character_)
}

## 留一数据集验证的实测数字，随包分发（与引擎用的是同一份 JSON）。
.projection_validation <- function() {
  p <- tryCatch(.ewb_extdata("reference_projection_validation.rds"),
                error = function(e) NULL)
  if (!is.null(p) && file.exists(p)) return(readRDS(p))
  NULL
}

.projection_caveats <- function() {
  list(
    zh = c("\u8fd9\u662f**\u6295\u5f71\u4e0d\u662f\u5206\u7c7b**\uff1a\u51e0\u4f55\u4e0a\u4efb\u4f55\u6837\u672c\u90fd\u4f1a\u843d\u5230\u6700\u8fd1\u7684\u7fa4\uff0c\u6240\u4ee5\u7ed3\u679c\u5fc5\u987b\u540c\u65f6\u770b OOD \u4e0e\u6295\u7968\u6b67\u4e49\uff1b\u4e24\u8005\u4efb\u4e00\u4e0d\u8fc7\u5173\u5c31\u5f03\u6743\uff0c\u4e0d\u7ed9\u6807\u7b7e\u3002",
           "\u7559\u4e00\u6570\u636e\u96c6\u9a8c\u8bc1\uff08\u4e0d\u662f\u7559\u4e00\u6837\u672c\uff09\uff1a\u540c\u6279\u6837\u672c\u7684\u76f8\u4f3c\u5ea6\u4f1a\u88ab\u5f53\u6210\u7cbe\u5ea6\u3002",
           "\u53c2\u8003\u91cc\u767d\u7ec6\u80de\u6837\u672c\u53ea\u6709 3 \u4f8b\u4e14\u5168\u5728\u540c\u4e00\u4e2a\u6570\u636e\u96c6\uff0c\u56e0\u6b64\u5bf9\u767d\u7ec6\u80de\u5236\u54c1\u7684\u8bc6\u522b\u80fd\u529b**\u6ca1\u6709\u7ecf\u8fc7\u6709\u6548\u9a8c\u8bc1**\u3002",
           "\u79e9\u7a7a\u95f4\u5bf9\u4efb\u610f\u5355\u8c03\u5e73\u53f0\u504f\u5dee\u514d\u75ab\uff0c\u6545\u6821\u6b63\u5728\u8be5\u7a7a\u95f4\u4e0d\u6539\u53d8\u4efb\u4f55\u8bfb\u6570\u3002"),
    en = c("This is a **projection, not a classifier**: geometrically every sample lands next to some reference group, so OOD and vote ambiguity are mandatory; a label is abstained if either fails.",
           "Validation is leave-one-**dataset**-out, not leave-one-sample-out.",
           "Only 3 leukocyte reference samples exist, all in one dataset, so recognition of leukocyte preparations is **not validated**.",
           "Rank space is immune to monotone platform shifts, so correction cannot change it.")
  )
}
