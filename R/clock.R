## clock.R — 储存时钟：一个单位偏离新鲜状态多少天。
##
## 对应引擎：toolkit_tools.run_clock / _clock_reference / _clock_module /
## _attach_validation。
##
## 四条口径：
##   1. 一律**先用预训练模型预测**，注释里的天数只用来做外部验证。用标签
##      训练再用留一误差自评，报出来的准确度是自己评自己。
##   2. miRNA 与代谢组各有同层模型；转录组/蛋白组没有任何天数标签，只能用
##      「模块分数 → 天数」的跨层外推，属低置信路径。
##   3. 平台整体平移**逐样本**校正（中位 z 超过 1 才校），不是整批一个数。
##   4. 预测截断在训练集覆盖的天数范围内，标记 clipped；外推没有依据。

## 逐特征稳健 z：中位 + 1.4826×MAD，退化到标准差，再退化则该特征无信息。
## 与 toolkit_core.robust_row_z 逐位相同。不能用 IQR 兜底成 1：代谢组里
## 很多特征四分位距正好为 0，兜底成 1 会让该特征保留原始浓度量级，模块分数
## 被推到 1e5（训练时又 z 了一次，看不出来；跨层一用全部归零）。
.robust_row_z <- function(x, clip = 5) {
  ctr <- apply(x, 1, stats::median, na.rm = TRUE)
  dev <- abs(sweep(x, 1, ctr, "-"))
  madv <- apply(dev, 1, stats::median, na.rm = TRUE)
  s <- 1.4826 * madv
  sdv <- apply(x, 1, .nansd0)
  s[!is.finite(s) | s <= 0] <- sdv[!is.finite(s) | s <= 0]
  ok <- is.finite(s) & s > 0 & is.finite(ctr)
  s[!ok | !is.finite(s) | s <= 0] <- 1
  z <- sweep(sweep(x, 1, ctr, "-"), 1, s, "/")
  list(z = pmax(pmin(z, clip), -clip), ok = ok)
}

## 把特征折叠成储存损伤模块分数：模块内稳健 z 的均值。
.module_fold <- function(features, x, space, modules, min_features = 2L) {
  rz <- .robust_row_z(x)
  Z <- rz$z
  ok <- rz$ok
  f <- if (identical(space, "metabolite")) tolower(trimws(features)) else toupper(trimws(features))
  scores <- list()
  cover <- list()
  for (key in names(modules)) {
    m <- modules[[key]]
    want <- if (identical(space, "metabolite")) m$metabolites else m$genes
    want <- unlist(want, use.names = FALSE)
    if (is.null(want) || !length(want)) next
    want <- if (identical(space, "metabolite")) tolower(trimws(want)) else toupper(trimws(want))
    idx <- which(f %in% want & ok)
    if (length(idx) < min_features) next
    scores[[key]] <- apply(Z[idx, , drop = FALSE], 2, .nanmean0)
    cover[[key]] <- length(idx)
  }
  list(scores = scores, cover = cover)
}

#' Storage clock: how many days a unit has drifted from fresh
#'
#' Projects each sample onto a pre-trained storage clock and returns predicted
#' days. When the annotation carries real storage days, those days are used to
#' **validate** the predictions (MAE, bias, rho) rather than to fit them.
#'
#' Two paths exist. miRNA and metabolomics input use a same-layer model. Gene
#' and protein input have no day labels anywhere in the reference, so they are
#' bridged through storage-lesion module scores learned on metabolomics: that is
#' a cross-layer extrapolation and is reported as such, with a resolution check
#' rather than a confident number.
#'
#' @param x Numeric matrix or data.frame, features in rows and samples in
#'   columns.
#' @param days Optional numeric vector of true storage days, named by sample.
#'   When given, predictions are compared against it.
#' @param space `"auto"` (default) detects miRNA / metabolite from the feature
#'   names; otherwise one of `"mirna"`, `"metabolite"`, `"gene"`.
#' @param path Which route to take when both are available: `"auto"` (default)
#'   uses the same-layer reference model whenever the space has one and falls
#'   back to the cross-layer module bridge otherwise; `"reference"` forces the
#'   reference route, `"module"` forces the bridge. Forcing is mainly useful for
#'   auditing the extrapolation path, which is an extrapolation either way.
#'
#' @return A data.frame with `sample`, `day`, `lo`, `hi`, `coverage`,
#'   `raw`, `clipped`, `true_day`. Attributes carry the mode, matched features,
#'   the model's own MAE and the caveats.
#'
#' @examples
#' \dontrun{
#' ery_clock(read.delim("metabolites.tsv", row.names = 1))
#' }
#' @export
ery_clock <- function(x, days = NULL, space = "auto",
                      path = c("auto", "reference", "module")) {
  path <- match.arg(path)
  cl <- .clean_matrix(x)
  if (identical(space, "auto")) {
    id <- ery_id_space(cl$features)
    space <- switch(id$verdict,
                    mirna = "mirna",
                    metabolite = "metabolite",
                    gene = "gene",
                    protein = "gene",
                    "gene")
  }
  sp <- ery_to_space(cl$x, cl$features, space)
  ## 与引擎 _prepare 一致：**不做** as_log2（时钟模型训练时用的是原尺度），
  ## 缺失用行中位数填补，并单独记下哪些格子原本有值——否则"整行缺失"的样本
  ## 会被当成"测到了、值等于中位数"，给出看似正常却没有依据的天数。
  xf <- sp$x
  finite <- is.finite(xf)
  for (i in seq_len(nrow(xf))) {
    row <- xf[i, ]
    if (any(!is.finite(row))) {
      fill <- stats::median(row, na.rm = TRUE)
      row[!is.finite(row)] <- if (is.finite(fill)) fill else 0
      xf[i, ] <- row
    }
  }
  prep <- list(names = sp$features, x = xf, finite = finite, space = space)

  use_ref <- switch(path,
                    reference = TRUE,
                    module = FALSE,
                    space %in% c("mirna", "metabolite"))
  out <- if (use_ref) .clock_reference(prep) else .clock_module(prep)
  ## 判断是否加验证要看模式，不看 out$ok：out 是 data.frame，没有 ok 列，
  ## 用 isTRUE(out$ok) 会恒为 FALSE，验证块永远不会挂上。
  if (!is.null(days) && attr(out, "mode") %in% c("reference", "extrapolated")) {
    out <- .clock_attach_validation(out, days)
  }
  out
}

.clock_reference <- function(prep) {
  refs <- .ewb_read_rds("toolkit_refs.rds")$clock
  models <- Filter(function(v) identical(v$space, prep$space), refs)
  if (!length(models)) {
    stop("No in-house storage clock for space '", prep$space, "'. Available: ",
         paste(sort(unique(vapply(refs, function(v) v$space, ""))), collapse = ", "),
         call. = FALSE)
  }
  ## 选重合特征最多的模型
  pos <- seq_along(prep$names)
  names(pos) <- prep$names
  n_hit <- vapply(models, function(m) {
    sum(unlist(m$features, use.names = FALSE) %in% prep$names)
  }, 0)
  key <- names(models)[which.max(n_hit)]
  model <- models[[key]]
  want <- unlist(model$features, use.names = FALSE)
  hit <- which(want %in% prep$names)
  if (length(hit) < 10) {
    stop("Only ", length(hit), " of ", length(want),
         " reference features are present; a projection needs at least 10.",
         call. = FALSE)
  }

  ridx <- pos[want[hit]]
  v <- prep$x[ridx, , drop = FALSE]
  covered <- colMeans(prep$finite[ridx, , drop = FALSE])
  mu <- .num_vec(model$mu)[hit]
  sdv <- .num_vec(model$sd)[hit]
  sdv[!is.finite(sdv) | sdv <= 0] <- 1
  w <- .num_vec(model$w)[hit]
  Z <- sweep(sweep(v, 1, mu, "-"), 1, sdv, "/")

  ## 逐样本平台平移：中位 z 超过 1 才校，且各自校各自的
  off <- apply(Z, 2, stats::median, na.rm = TRUE)
  off[!is.finite(off) | abs(off) <= 1] <- 0
  zz <- as.numeric(crossprod(w, sweep(Z, 2, off, "-"))) + as.numeric(model$intercept)
  zz[!is.finite(zz)] <- 0
  pred <- expm1(pmax(pmin(zz, 10), -8))
  lo_d <- as.numeric(model$day_min)
  hi_d <- as.numeric(model$day_max)
  clipped <- pmax(pmin(pred, hi_d), lo_d)
  m <- model$metrics
  q10 <- if (is.null(m$resid_q10)) -as.numeric(m$mae_days) else as.numeric(m$resid_q10)
  q90 <- if (is.null(m$resid_q90)) as.numeric(m$mae_days) else as.numeric(m$resid_q90)

  n <- ncol(prep$x)
  rows <- lapply(seq_len(n), function(j) {
    if (covered[j] < 0.5) {
      return(data.frame(sample = colnames(prep$x)[j], day = NA_real_,
                        lo = NA_real_, hi = NA_real_, coverage = round(covered[j], 2),
                        raw = NA_real_, clipped = NA, true_day = NA_real_,
                        note = "insufficient_features", stringsAsFactors = FALSE))
    }
    data.frame(sample = colnames(prep$x)[j], day = round(clipped[j], 1),
               lo = round(max(0, clipped[j] + q10), 1),
               hi = round(clipped[j] + q90, 1),
               coverage = round(covered[j], 2), raw = round(pred[j], 1),
               clipped = abs(pred[j] - clipped[j]) > 0.5, true_day = NA_real_,
               note = NA_character_, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  res <- res[order(is.na(res$day), ifelse(is.na(res$day), 0, res$day)), , drop = FALSE]

  attr(res, "mode") <- "reference"
  attr(res, "space") <- prep$space
  attr(res, "model") <- key
  attr(res, "dataset") <- model$dataset
  attr(res, "matched_features") <- length(hit)
  attr(res, "total_features") <- length(want)
  attr(res, "day_range") <- c(lo_d, hi_d)
  attr(res, "n_offset") <- sum(abs(off) > 0)
  attr(res, "metrics") <- m[c("mae_days", "mae_days_nested", "rho", "baseline_mae_days")]
  attr(res, "caveats") <- list(
    zh = c("\u5929\u6570\u662f**\u9884\u8bad\u7ec3\u6a21\u578b**\u6295\u5f71\u51fa\u6765\u7684\uff0c\u4e0d\u662f\u5728\u8fd9\u4efd\u6570\u636e\u4e0a\u62df\u5408\u7684\uff1b\u51c6\u786e\u5ea6\u4ee5\u8be5\u6a21\u578b\u7684 MAE \u4e3a\u51c6\u3002",
           sprintf("\u9884\u6d4b\u88ab\u622a\u65ad\u5728\u8bad\u7ec3\u96c6\u8986\u76d6\u7684 %s\u2013%s \u5929\uff1a\u8d85\u51fa\u8fd9\u4e2a\u8303\u56f4\u662f\u5916\u63a8\uff0c\u6a21\u578b\u6ca1\u6709\u4f9d\u636e\u3002", lo_d, hi_d),
           sprintf("\u547d\u4e2d %d/%d \u4e2a\u53c2\u8003\u7279\u5f81\u3002", length(hit), length(want)),
           "\u53c2\u8003\u6a21\u578b\u53ea\u5728\u5355\u4e00\u5e73\u53f0/\u5355\u4e00\u6279\u6b21\u4e0a\u8bad\u7ec3\u8fc7\uff0c\u6362\u5e73\u53f0\u3001\u6362\u6297\u51dd\u5242\u3001\u6362\u4fdd\u5b58\u6db2\u90fd\u4f1a\u5e26\u6765\u7cfb\u7edf\u504f\u5dee\u3002"),
    en = c("The days are projected by a **pre-trained** clock, not fitted on this data; its MAE is the accuracy you get.",
           sprintf("Predictions are clipped to the training range %s-%s days; beyond it the model extrapolates without support.", lo_d, hi_d),
           sprintf("%d/%d reference features matched.", length(hit), length(want)),
           "The reference clock was trained on one platform and one batch; a different platform, anticoagulant or storage solution adds systematic bias.")
  )
  res
}

.clock_module <- function(prep) {
  mod <- .ewb_read_rds("toolkit_refs.rds")$clock$module
  if (is.null(mod)) stop("No module clock in the reference bundle.", call. = FALSE)
  if (!prep$space %in% c("gene", "metabolite")) {
    stop("Cross-layer bridging only applies to gene or metabolite input.", call. = FALSE)
  }
  mf <- .module_fold(prep$names, prep$x, prep$space, mod$modules, 2L)
  order_k <- unlist(mod$module_order, use.names = FALSE)
  order_k <- order_k[order_k %in% names(mf$scores)]
  if (length(order_k) < 3) {
    stop("Only ", length(order_k), " storage-lesion modules could be mapped; ",
         "at least 3 are needed to extrapolate.", call. = FALSE)
  }
  M <- do.call(cbind, mf$scores[order_k])
  if (is.null(dim(M))) M <- matrix(M, nrow = 1)
  full_order <- unlist(mod$module_order, use.names = FALSE)
  pick <- match(order_k, full_order)
  mu <- .num_vec(mod$mu)[pick]
  sdv <- .num_vec(mod$sd)[pick]
  sdv[!is.finite(sdv) | sdv <= 0] <- 1
  w <- .num_vec(mod$w)[pick]
  lo_d <- as.numeric(mod$day_min)
  hi_d <- as.numeric(mod$day_max)
  cal <- mod$calibration
  cal_x <- .num_vec(cal$x)
  cal_y <- .num_vec(cal$y)

  to_days <- function(Mx) {
    Z <- sweep(sweep(Mx, 2, mu, "-"), 2, sdv, "/")
    lin <- as.numeric(Z %*% w) + as.numeric(mod$intercept)
    lin[!is.finite(lin)] <- 0
    p <- expm1(pmax(pmin(lin, 10), -8))
    if (length(cal_x) >= 2) {
      stats::approx(cal_x, cal_y, xout = p, rule = 2)$y   # 等渗校准
    } else {
      p
    }
  }

  center <- identical(as.character(mod$centering), "cross_layer")
  Mcentered <- if (center) sweep(M, 2, colMeans(M, na.rm = TRUE), "-") else M
  adj <- to_days(Mcentered)
  clipped <- pmax(pmin(adj, hi_d), lo_d)
  m <- mod$metrics
  mae <- if (!is.null(m$nested_mae_days)) as.numeric(m$nested_mae_days) else
    if (!is.null(m$mae_days)) as.numeric(m$mae_days) else 10

  free <- adj[abs(adj - clipped) <= 0.5 & is.finite(adj)]
  frac_clipped <- mean(abs(adj - clipped) > 0.5)
  saturated <- frac_clipped >= 0.25
  spread_free <- if (length(free) >= 3) diff(range(free)) else 0
  resolvable <- spread_free > mae && !saturated

  res <- data.frame(
    sample = colnames(prep$x),
    day = round(clipped, 1),
    lo = round(pmax(0, clipped - mae), 1),
    hi = round(pmin(hi_d, clipped + mae), 1),
    raw = round(adj, 1),
    clipped = abs(adj - clipped) > 0.5,
    coverage = NA_real_, true_day = NA_real_, note = NA_character_,
    stringsAsFactors = FALSE
  )
  res <- res[order(res$day), , drop = FALSE]

  attr(res, "mode") <- "extrapolated"
  attr(res, "space") <- prep$space
  attr(res, "model") <- "module"
  attr(res, "matched_modules") <- length(order_k)
  attr(res, "module_coverage") <- unlist(mf$cover[order_k], use.names = TRUE)
  attr(res, "day_range") <- c(lo_d, hi_d)
  attr(res, "metrics") <- list(mae_days = mae, nested_mae_days = m$nested_mae_days,
                               bias_days = m$bias_days, rho = m$rho,
                               calibrated = length(cal_x) >= 2, centered = center)
  attr(res, "resolvable") <- resolvable
  attr(res, "saturated") <- saturated
  attr(res, "frac_clipped") <- round(frac_clipped, 3)
  attr(res, "caveats") <- list(
    zh = c("\u8fd9\u662f**\u8de8\u5c42\u5916\u63a8**\uff1a\u6a21\u5757\u5230\u5929\u6570\u7684\u7cfb\u6570\u5728\u4ee3\u8c22\u7ec4\u4e0a\u5b66\u5230\uff0c\u5957\u5230\u522b\u7684\u7ec4\u5b66\u5c42\uff1b\u6a21\u5757\u662f\u751f\u7269\u5b66\u5355\u5143\uff0c\u6240\u4ee5\u80fd\u8de8\u5c42\uff0c\u8fd9\u4e5f\u662f\u5b83\u552f\u4e00\u7684\u4f9d\u636e\u3002",
           sprintf("\u547d\u4e2d %d \u4e2a\u50a8\u5b58\u635f\u4f24\u6a21\u5757\uff1b\u6a21\u578b\u8bef\u5dee %.1f \u5929\u3002", length(order_k), mae),
           if (center) "\u672c\u6279\u6a21\u5757\u5747\u503c\u5df2\u5bf9\u9f50\u5230\u8bad\u7ec3\u4e2d\u5fc3\uff0c\u53ea\u4fdd\u7559\u6837\u672c\u4e4b\u95f4\u7684\u76f8\u5bf9\u5dee\u5f02\u2014\u2014\u5c42\u95f4\u6574\u4f53\u504f\u79fb\u4e0d\u662f\u5929\u6570\u3002" else "\u672a\u505a\u8de8\u5c42\u4e2d\u5fc3\u5316\u3002",
           if (!resolvable) "\u672c\u6279\u6837\u672c\u7684\u6a21\u5757\u5dee\u5f02\u5c0f\u4e8e\u6a21\u578b\u8bef\u5dee\uff0c\u5929\u6570\u5dee\u5f02\u4e0d\u53ef\u5206\u8fa8\uff0c\u8bf7\u5f53\u4f5c\u540c\u4e00\u6c34\u5e73\u3002" else "\u672c\u6279\u5dee\u5f02\u5927\u4e8e\u6a21\u578b\u8bef\u5dee\uff0c\u53ef\u5206\u8fa8\u3002",
           if (saturated) "\u8d85\u8fc7\u56db\u5206\u4e4b\u4e00\u7684\u6837\u672c\u88ab\u622a\u65ad\u5728\u8bad\u7ec3\u533a\u95f4\u8fb9\u7f18\uff0c\u8bf4\u660e\u5b83\u4eec\u8fdc\u5728\u6a21\u578b\u8303\u56f4\u5185\uff0c\u5929\u6570\u53ea\u80fd\u7b97\u5916\u63a8\u3002" else ""),
    en = c("This is a **cross-layer extrapolation**: the module-to-day coefficients were learned on metabolomics and applied to another layer.",
           sprintf("%d storage-lesion modules matched; model error %.1f days.", length(order_k), mae),
           if (center) "Batch module means were aligned to the training centre; only relative differences are kept - a layer shift is not extra storage time." else "No cross-layer centring applied.",
           if (!resolvable) "Module differences in this batch are smaller than the model error, so the day differences are not resolvable." else "Batch differences exceed the model error.",
           if (saturated) "More than a quarter of samples are clipped at the training edge, so they sit well outside the model's range." else "")
  )
  res
}

.clock_attach_validation <- function(res, days) {
  truth <- as.numeric(days[res$sample])
  res$true_day <- truth
  ok <- is.finite(res$day) & is.finite(truth)
  if (sum(ok) < 4) {
    attr(res, "validation") <- list(ok = FALSE, n = sum(ok),
                                    note = "Too few samples with true days to report an error.")
    return(res)
  }
  err <- res$day[ok] - truth[ok]
  attr(res, "validation") <- list(
    ok = TRUE, n = sum(ok),
    mae_days = round(mean(abs(err)), 2),
    rmse_days = round(sqrt(mean(err^2)), 2),
    bias_days = round(mean(err), 2),
    rho = round(stats::cor(res$day[ok], truth[ok], method = "spearman"), 3),
    true_range = round(range(truth[ok]), 1),
    note = paste("These days are pre-trained predictions, not fitted values;",
                 "the MAE above is their external error on this data.")
  )
  attr(res, "mode") <- "reference_verified"
  res
}
