## profile.R — 一次调用把样本表征的各个模块跑完。
##
## 对应引擎的 /api/toolkit/analyze_all。设计上和引擎一样：**一个模块失败不
## 影响其他模块**，失败的那个把原因写进结果，而不是让整次调用报错——读者
## 需要知道"这一项为什么没出来"。

#' Profile a sample across the erythroid modules
#'
#' Runs the available sample-profiling modules in sequence on one matrix and
#' returns each block, or the reason a block does not apply. A module that
#' cannot run is recorded with its message instead of aborting the whole call:
#' knowing *why* the storage clock does not apply is part of the answer.
#'
#' @param x Numeric matrix or data.frame, features in rows and samples in
#'   columns.
#' @param layer Omics layer for the identity reference: `"transcriptome"`
#'   (default) or `"proteome"`.
#' @param space Feature space for the clock and pathway modules: `"auto"`
#'   (default), `"gene"`, `"mirna"` or `"metabolite"`.
#' @param days Optional named vector of true storage days, used only to
#'   validate the clock, never to fit it.
#'
#' @return A named list with `purity`, `transfer`, `clock`, `pathways`,
#'   `provenance` and `caveats`. Each module block carries its own attributes;
#'   a module that failed has `$error` set and nothing else.
#'
#' @examples
#' \dontrun{
#' p <- ery_profile(read.delim("matrix.tsv", row.names = 1))
#' p$purity
#' p$clock
#' }
#' @export
ery_profile <- function(x, layer = c("transcriptome", "proteome"),
                        space = "auto", days = NULL) {
  layer <- match.arg(layer)
  out <- list()

  ## 身份与混入：这是入口模块，失败要如实记录而不是静默跳过
  out$purity <- .safe_run(function() {
    r <- ery_purity(x, layer = layer, level = "high")
    ## 置信度按同层参照量重算一次，high 只是初值
    lv <- ery_confidence(attr(r, "n_ref_samples"), attr(r, "n_leu_markers"),
                         any(r$conflict))
    if (!identical(lv, "high")) {
      r <- ery_purity(x, layer = layer, level = lv)
    }
    r
  })

  out$transfer <- .safe_run(function() ery_transfer(x))
  out$clock <- .safe_run(function() ery_clock(x, days = days, space = space))
  out$pathways <- .safe_run(function() ery_pathways(x, mode = "within", space = space))
  ## 条件谱匹配与红系衰老：前者按 id 空间自动选画像库（脂质矩阵会因没有共享
  ## 特征而自然弃权），后者只在行名能解析成脂质分子种时才给 A 档读数。
  out$attribution <- .safe_run(function() ery_attribution(x, n_perm = 100L))
  out$age <- .safe_run(function() ery_age(x, n_perm = 200L))

  out$provenance <- ery_provenance()
  out$caveats <- list(
    zh = c("\u6bcf\u4e2a\u6a21\u5757\u7684\u4e3b\u8bfb\u6570\u90fd\u5728\u5b83\u81ea\u5df1\u7684\u5c5e\u6027\u91cc\uff08attr(result, \"mode\") / attr(result, \"level\") \u7b49\uff09\uff0c\u4e0d\u8981\u53ea\u770b\u6570\u503c\u3002",
           "\u4e3b\u8bfb\u6570\u662f\u79e9\u4e00\u81f4\u6027\u8bc4\u5206\uff0c\u4e0d\u662f\u7ec6\u80de\u8ba1\u6570\u7eaf\u5ea6\uff1a\u6210\u719f\u7ea2\u7ec6\u80de\u65e0\u6838\uff0c\u77e9\u9635\u6d4b\u5230\u7684\u662f\u7ea2\u7cfb\u6765\u6e90 RNA\u3002",
           "\u8eab\u4efd\u8f6c\u79fb\u662f\u6295\u5f71\u4e0d\u662f\u5206\u7c7b\uff0cOOD \u4e0e\u6b67\u4e49\u4e24\u9053\u95f8\u95e8\u4efb\u4e00\u4e0d\u8fc7\u5c31\u5f03\u6743\uff1b\u767d\u7ec6\u80de\u53c2\u7167\u53ea\u6709 3 \u4f8b\u4e14\u540c\u5c5e\u4e00\u4e2a\u6570\u636e\u96c6\u3002",
           "\u50a8\u5b58\u65f6\u949f\u5728\u8f6c\u5f55\u7ec4/\u86cb\u767d\u7ec4\u4e0a\u8d70\u7684\u662f\u8de8\u5c42\u6a21\u5757\u6865\uff0c\u5c5e\u5916\u63a8\uff1b\u5929\u6570\u622a\u65ad\u5728\u8bad\u7ec3\u533a\u95f4\u5185\u3002",
           "\u901a\u8def\u6307\u7eb9 within \u6a21\u5f0f\u53ea\u56de\u7b54\u300c\u6837\u672c\u5185\u76f8\u5bf9\u300d\uff0c\u4e0d\u56de\u7b54\u300c\u8de8\u6837\u672c\u9ad8\u4f4e\u300d\u3002",
           "\u6761\u4ef6\u5f52\u56e0\u662f\u300c\u6700\u50cf\u54ea\u4e00\u7c7b\u6761\u4ef6\u300d\u4e0d\u662f\u8bca\u65ad\uff1b\u975e\u7ea2\u7cfb\u6837\u672c\u4f1a\u5148\u88ab\u6807\u8bb0\u4e3a\u9519\u4f4d\u63d0\u95ee\u3002",
           "\u7ea2\u7cfb\u8870\u8001\u7684 0\u2013100 \u662f\u672c\u6279\u6837\u672c\u5185\u7684\u76f8\u5bf9\u79e9\u6295\u5f71\uff0c\u4e0d\u662f\u7edd\u5bf9\u5b58\u50a8\u5929\u6570\uff1b\u8de8\u5c42\u5148\u9a8c\u5916\u63a8\u6807\u4e3a unverifiable\u3002"),
    en = c("Each module keeps its primary reading in its own attributes (attr(result, \"mode\") etc.); do not read the numbers alone.",
           "The primary reading is a rank-agreement score, not cell-count purity: mature red cells are enucleated.",
           "Label transfer is a projection, not a classifier; a label is abstained if either gate fails, and only 3 leukocyte references exist.",
           "On gene/protein input the storage clock is a cross-layer module bridge, i.e. an extrapolation, clipped to the trained range.",
           "The within-sample pathway fingerprint answers only 'high relative to other pathways in this sample'.",
           "Condition matching names the closest condition class, not a diagnosis; non-erythroid samples are flagged as the wrong question first.",
           "The 0-100 ageing score is a within-batch rank projection, not an absolute storage age; the cross-layer prior tier is marked unverifiable.")
  )
  attr(out, "engine_version") <- out$provenance$engine_version
  ## 设 class，否则 print() 走的是 list 打印，print.ery_profile 永远不触发
  class(out) <- c("ery_profile", "list")
  out
}

.safe_run <- function(f) {
  tryCatch(f(), error = function(e) list(error = conditionMessage(e)))
}

#' Print a profiling result
#'
#' @param x An object returned by [ery_profile()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.ery_profile <- function(x, ...) {
  cat("<ErythroWorkbench profile>\n")
  cat("  engine:", attr(x, "engine_version") %||% "unknown", "\n\n")

  p <- x$purity
  if (!is.null(p$error)) {
    cat("identity : failed -", p$error, "\n")
  } else {
    cat(sprintf("identity : %d samples, layer %s (%s)\n", nrow(p),
                attr(p, "layer"), attr(p, "layer_label")))
    hd <- p[, c("sample", "headline", "ery_consistency", "identity")]
    print(utils::head(hd, 10), row.names = FALSE)
  }

  tr <- x$transfer
  cat("\nlabel transfer:\n")
  if (!is.null(tr$error)) {
    cat("  failed -", tr$error, "\n")
  } else {
    cat(sprintf("  coverage %.1f%%, k = %d, gates: ood %.3f / margin %.3f\n",
                100 * attr(tr, "coverage"), attr(tr, "k"),
                attr(tr, "ood_gate"), attr(tr, "margin_gate")))
    print(utils::head(tr[, c("sample", "label", "abstain")], 10), row.names = FALSE)
  }

  cl <- x$clock
  cat("\nstorage clock:\n")
  if (!is.null(cl$error)) {
    cat("  not applicable -", cl$error, "\n")
  } else {
    cat(sprintf("  mode %s, model %s, matched %s features\n",
                attr(cl, "mode"), attr(cl, "model"),
                attr(cl, "matched_features") %||% attr(cl, "matched_modules")))
    print(utils::head(cl[, c("sample", "day", "lo", "hi")], 10), row.names = FALSE)
  }

  pw <- x$pathways
  cat("\npathways:\n")
  if (!is.null(pw$error)) {
    cat("  failed -", pw$error, "\n")
  } else {
    cat(sprintf("  %d scored, %d excluded, mode %s\n",
                pw$n_pathways, if (is.null(pw$excluded)) 0L else nrow(pw$excluded),
                pw$mode))
    print(utils::head(pw$table[, c("pathway_id", "n_measured", "consistency")], 6),
          row.names = FALSE)
  }
  at <- x$attribution
  cat("\ncondition match:\n")
  if (!is.null(at$error)) {
    cat("  not applicable -", at$error, "\n")
  } else {
    ok <- Filter(function(r) isTRUE(r$ok), at$results)
    cat(sprintf("  %s space, %d centroids, %d conditions; %d/%d samples placed\n",
                at$space, at$n_centroids, at$n_conditions, length(ok), at$n_samples))
    for (r in utils::head(ok, 6)) {
      cat(sprintf("    %-12s %-28s rho=%.3f p=%s\n", r$sample, r$top$label,
                  r$top$best_rho,
                  if (is.null(r$top$p)) "NA" else sprintf("%.4f", r$top$p)))
    }
  }

  ag <- x$age
  cat("\nred-cell ageing:\n")
  if (!is.null(ag$error)) {
    cat("  not applicable -", ag$error, "\n")
  } else if (isTRUE(ag$ok)) {
    cat(sprintf("  %d placed, %d refused; median senescence %.1f (%s); LOSO MAE %.2f\n",
                ag$summary$n_ok, ag$summary$n_refused,
                ag$summary$median_senescence_score, ag$summary$median_zone$en,
                ag$validation$loso_mae))
  } else {
    cat("  abstained -", ag$reason_code, "\n")
    pe <- ag$prior_extrapolation
    if (!is.null(pe) && isTRUE(pe$ok)) {
      cat(sprintf("  cross-layer tier: %d samples, grade %s, confidence %.3f (unverifiable)\n",
                  length(pe$samples), pe$grade, pe$confidence))
    }
  }

  invisible(x)
}
