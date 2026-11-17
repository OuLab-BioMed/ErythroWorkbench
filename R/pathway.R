## pathway.R — 单样本通路扰动指纹。
##
## 对应引擎：toolkit_tools.run_pathway_fingerprint + toolkit_core.pathway_*。
##
## 三条口径：
##   1. within 是**样本内跨特征**标准化（高=相对其他通路偏高），cohort 是
##      **跨样本按特征**标准化（高=高于本批均值）。两者不能混读。
##   2. 打分只用量得到、且真正反映通路活性的成员；role 标了
##      abundance-excluded 的丰度载体（HBB / ALB / CA1 / CA2）留作机制注释，
##      不参与打分——HBB 单基因就能把整条通路的均值 z 抬到 26.9。
##   3. 排序按 max|分数| × √(n_measured)。不乘 √n 时，成员最少的通路在纯
##      随机噪声下也有一半概率排第一。

.pathway_kb <- function() {
  if (!is.null(.ewb_cache[["pw_kb"]])) return(.ewb_cache[["pw_kb"]])
  pw <- .ewb_read_tsv("pathways.tsv")
  pg <- .ewb_read_tsv("pathway_genes.tsv")
  pm <- .ewb_read_tsv("pathway_metabolites.tsv")
  kb <- list(pathways = pw, genes = pg, mets = pm)
  .ewb_cache[["pw_kb"]] <- kb
  kb
}

## 返回 {pathway_id: [成员名...]}，与 toolkit_core.pathway_members 一致。
.pathway_members <- function(space, scoring_only = TRUE) {
  kb <- .pathway_kb()
  if (identical(space, "metabolite")) {
    nm <- ery_norm_metabolite(kb$mets$metabolite_name)
    keep <- !is.na(nm)
    m2p <- split(as.character(kb$mets$pathway_id)[keep], nm[keep])
  } else {
    sym <- toupper(trimws(as.character(kb$genes$gene_symbol)))
    role <- as.character(kb$genes$role)
    keep <- !is.na(sym) & nzchar(sym)
    if (scoring_only) {
      keep <- keep & !grepl("abundance-excluded", ifelse(is.na(role), "", role), fixed = TRUE)
    }
    m2p <- split(as.character(kb$genes$pathway_id)[keep], sym[keep])
  }
  ## m2p 是「成员 → 通路集合」的反查表，这里必须真反一次再返回。
  ## 直接把反查表返回去，会让所有通路的成员数静默变成 0。
  pids <- unlist(m2p, use.names = FALSE)
  mems <- rep(names(m2p), lengths(m2p))
  out <- split(mems, pids)
  lapply(out, function(v) sort(unique(v)))
}

#' Single-sample pathway perturbation fingerprint
#'
#' Scores each curated pathway as the mean z of its measured, scoreable members.
#' In `"within"` mode every feature is standardised inside each sample, so a
#' positive value means the pathway sits above the other pathways in that same
#' sample. In `"cohort"` mode every feature is standardised across samples, so a
#' positive value means this sample is above the batch mean.
#'
#' Pathways whose curated members include abundance carriers (HBB, ALB, CA1,
#' CA2 and similar) are scored **without** those carriers: one highly abundant
#' transcript would otherwise drag the whole pathway's mean z to absurd values.
#'
#' @param x Numeric matrix or data.frame, features in rows and samples in
#'   columns.
#' @param mode `"within"` (default, works on a single sample) or `"cohort"`.
#' @param space `"auto"` (default), `"gene"` or `"metabolite"`.
#'
#' @return A list with `table` (data.frame: `pathway_id`, `name_cn`, `name_en`,
#'   `n_members`, `n_measured`, `consistency` and one column per sample) and
#'   `excluded` (pathways that could not be scored, with the reason).
#'   Pathways are ordered by `max|score| * sqrt(n_measured)`.
#'
#' @examples
#' \dontrun{
#' res <- ery_pathways(read.delim("matrix.tsv", row.names = 1))
#' head(res$table)
#' }
#' @export
ery_pathways <- function(x, mode = c("within", "cohort"), space = "auto") {
  mode <- match.arg(mode)
  cl <- .clean_matrix(x)
  if (identical(space, "auto")) {
    id <- ery_id_space(cl$features)
    space <- switch(id$verdict, metabolite = "metabolite", gene = "gene", "gene")
  }
  sp <- ery_to_space(cl$x, cl$features, space)
  lg <- ery_as_log2(sp$x)
  L <- lg$x
  ## 通路打分**要**做 as_log2：HBB 这类高丰度基因在原始强度下会把列内 sd
  ## 撑大，99.5% 的特征 z 被压到 |z|<0.05，同时自己拿到 z=134。
  Z <- if (identical(mode, "within")) .within_sample_z(L) else .zscore_rows(L)
  Z[!is.finite(Z)] <- 0

  members <- .pathway_members(space, scoring_only = TRUE)
  all_members <- .pathway_members(space, scoring_only = FALSE)
  kb <- .pathway_kb()
  pw <- kb$pathways

  feat <- sp$features
  if (identical(space, "metabolite")) {
    fk <- ery_norm_metabolite(feat)
    mk <- function(v) ery_norm_metabolite(v)
  } else {
    fk <- toupper(feat)
    mk <- function(v) toupper(v)
  }
  pos <- match(fk, fk)                       # 已去重，位置即名
  names(pos) <- fk

  samples <- colnames(sp$x)
  rows <- list()
  for (pid in names(members)) {
    ms <- members[[pid]]
    idx <- unname(pos[mk(ms)])
    idx <- idx[!is.na(idx)]
    if (length(idx) < 2) next
    sub <- Z[idx, , drop = FALSE]
    sc <- colMeans(sub)
    ## 通路内一致性：每个样本上，成员朝同一方向的占比，取中位数
    cons <- apply(sub, 2, function(v) max(sum(v > 0), sum(v < 0)) / length(v))
    nm <- pw[pw$pathway_id == pid, , drop = FALSE]
    rows[[pid]] <- data.frame(
      pathway_id = pid,
      name_cn = if (nrow(nm)) nm$name_cn[1] else pid,
      name_en = if (nrow(nm)) nm$name_en[1] else pid,
      n_members = length(ms), n_measured = length(idx),
      consistency = round(stats::median(cons), 3),
      stringsAsFactors = FALSE
    )
    rows[[pid]][, samples] <- round(as.numeric(sc), 3)
  }
  if (!length(rows)) {
    stop("None of the curated pathways has two or more measured members in ",
         "this matrix.", call. = FALSE)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  ## 排序：max|分数| × √(n_measured)，抵消小通路均值的抽样方差
  mx <- apply(tab[, samples, drop = FALSE], 1, function(v) max(abs(v)))
  ord <- order(-(mx * sqrt(tab$n_measured)))
  tab <- tab[ord, , drop = FALSE]

  ## 缺席的通路要显式列出：读者有权知道「没出现」是可打分成员不足，
  ## 而不是「没有变化」
  scored <- tab$pathway_id
  excluded <- lapply(setdiff(names(all_members), scored), function(pid) {
    ms <- all_members[[pid]]
    sc <- members[[pid]]
    if (is.null(sc)) sc <- character(0)
    hit <- sum(!is.na(unname(pos[mk(sc)])))
    nm <- pw[pw$pathway_id == pid, , drop = FALSE]
    data.frame(
      pathway_id = pid,
      name_cn = if (nrow(nm)) nm$name_cn[1] else pid,
      name_en = if (nrow(nm)) nm$name_en[1] else pid,
      n_curated = length(ms), n_scorable = length(sc), n_measured = hit,
      reason = if (length(sc) < 2) "too_few_scorable_members" else
        "not_measured_in_this_matrix",
      stringsAsFactors = FALSE
    )
  })
  excluded <- if (length(excluded)) do.call(rbind, excluded) else NULL
  if (!is.null(excluded)) excluded <- excluded[order(excluded$pathway_id), , drop = FALSE]

  list(table = tab, excluded = excluded, mode = mode, space = space,
       scale = lg$scale, n_features = length(feat),
       n_pathways = nrow(tab), samples = samples)
}
