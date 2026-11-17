## utils.R — 与 Python 引擎逐位对齐的计算原语。
##
## 这个文件里的每个函数都对应 app/backend/engines/toolkit_core.py 或
## sample_types.py 里的一处实现。**不要"改进"它们**：一旦与引擎分叉，
## 包和网站就会给出不同的数字，而这种差异不会报错、只会静默出现。
## 每处都注明了 Python 侧对应实现与容易踩的差异（尤其是 std 的自由度）。
##
## 零依赖：只用 base R。

## ---------------------------------------------------------------------------
## 内部：包级缓存环境（参照数据只读一次）
## ---------------------------------------------------------------------------
.ewb_cache <- new.env(parent = emptyenv())

## 读取随包分发的 extdata；找不到时给出可操作的报错，而不是返回 NULL。
.ewb_extdata <- function(name) {
  p <- system.file("extdata", name, package = "ErythroWorkbench")
  if (!nzchar(p)) {
    p <- file.path("inst", "extdata", name)   # 源码树里直接 source 时的退路
  }
  if (!file.exists(p)) {
    stop("Reference bundle missing: ", name,
         "\nReinstall the package, or regenerate it with data-raw/export_refs.R.",
         call. = FALSE)
  }
  p
}

.ewb_read_rds <- function(name) {
  key <- paste0("rds:", name)
  if (!is.null(.ewb_cache[[key]])) return(.ewb_cache[[key]])
  obj <- readRDS(.ewb_extdata(name))
  .ewb_cache[[key]] <- obj
  obj
}

.ewb_read_tsv <- function(name) {
  key <- paste0("tsv:", name)
  if (!is.null(.ewb_cache[[key]])) return(.ewb_cache[[key]])
  obj <- utils::read.delim(.ewb_extdata(name), stringsAsFactors = FALSE,
                           check.names = FALSE, quote = "")
  .ewb_cache[[key]] <- obj
  obj
}

## 把 jsonlite 读出的 list 转成数值向量，**NULL 必须保留成 NA 并占位**。
## 直接 unlist() 会把 NULL 元素整个丢掉，后面的值全体左移，末尾再补 NA——
## 缺失总数看起来一样，位置却错了。图谱里 373/381 个样本有缺失 marker，
## 这个 bug 会让相似度算在错位的向量上，而且不报任何错。
.num_vec <- function(x) {
  vapply(x, function(z) {
    if (is.null(z) || length(z) == 0L || !is.finite(suppressWarnings(as.numeric(z[1])))) {
      NA_real_
    } else {
      as.numeric(z[1])
    }
  }, 0, USE.NAMES = FALSE)
}

## 按 ASCII key 找层参照。**不要**改成按中文层名索引 list：
## 在非 UTF-8 locale 的 Windows 上，源码里的中文字面量与该索引比较会让 R
## 段错误（不是抛错，是进程直接挂掉）。key 由 data-raw/export_refs.R 写入，
## 中文名放在 $key_label 里只作显示。
.find_layer <- function(refs, key) {
  for (nm in names(refs$layers)) {
    lyr <- refs$layers[[nm]]
    if (identical(lyr$key, key)) return(lyr)
  }
  NULL
}

## ---------------------------------------------------------------------------
## 参照数据访问器
## ---------------------------------------------------------------------------

#' Reference bundle provenance
#'
#' Returns the fingerprints recorded when the reference bundle was exported
#' from the web-server engine: which engine version produced it, when, and the
#' md5 of every source file. Use this to state exactly which reference a
#' numeric result came from.
#'
#' @return A list with `exported_at`, `engine_version`, `r_version` and
#'   `sources` (named md5 vector).
#' @export
ery_provenance <- function() .ewb_read_rds("provenance.rds")

#' Marker panels and per-layer calibration
#'
#' Erythroid (24) and leukocyte (23) marker panels, the per-layer rank
#' consensus over in-house red-cell samples, and the calibration quantities
#' (`self_rho_p50`, `self_rho_p95`) that set the prior and the ceiling.
#'
#' @return The parsed reference list.
#' @export
ery_layer_reference <- function() .ewb_read_rds("layer_references.rds")

#' Fixed reference atlas
#'
#' 381 in-house samples profiled on the 47 lineage markers, with the per-sample
#' rank vector used for label transfer, plus the out-of-distribution gate.
#'
#' @return A list with `meta` (data.frame), `rank` (matrix samples x markers),
#'   `markers`, `ood_gate` and build metadata.
#' @export
ery_atlas <- function() .ewb_read_rds("reference_atlas.rds")

## ---------------------------------------------------------------------------
## 计算原语
## ---------------------------------------------------------------------------

## np.nanstd / np.nanmean 用 ddof = 0（总体标准差），R 的 sd() 用 ddof = 1。
## 直接换成 sd() 会让每个 z 值差一个 sqrt(n/(n-1)) 的因子，n 小时肉眼可见，
## 且不会报任何错。所以这里一律手写 ddof = 0。
.nanmean0 <- function(x) {
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  mean(x[ok])
}

.nanvar0 <- function(x) {
  ok <- is.finite(x)
  if (sum(ok) < 1L) return(NA_real_)
  v <- x[ok]
  mean((v - mean(v))^2)
}

.nansd0 <- function(x) {
  v <- .nanvar0(x)
  if (is.na(v)) NA_real_ else sqrt(v)
}

#' Column-wise z-score (within sample)
#'
#' Standardises each **sample** (column) across features, so that "high" means
#' high within that sample rather than measured deeper overall. Uses the
#' population standard deviation (ddof = 0) to match the engine.
#'
#' @param x A numeric matrix, features in rows and samples in columns.
#' @return A matrix of the same shape; zero-variance columns become NA.
#' @keywords internal
.within_sample_z <- function(x) {
  mu <- apply(x, 2, .nanmean0)
  sdv <- apply(x, 2, .nansd0)
  sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_
  sweep(sweep(x, 2, mu, "-"), 2, sdv, "/")
}

#' Row-wise z-score (across samples)
#'
#' The cohort view: standardises each feature across samples, answering
#' "high relative to this batch". Population SD (ddof = 0), zero-variance rows
#' become NA.
#'
#' @param x A numeric matrix, features in rows and samples in columns.
#' @return A matrix of the same shape.
#' @keywords internal
.zscore_rows <- function(x) {
  mu <- apply(x, 1, .nanmean0)
  sdv <- apply(x, 1, .nansd0)
  sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_
  sweep(sweep(x, 1, mu, "-"), 1, sdv, "/")
}

#' Rank a vector, keeping NA
#'
#' Same as `rank(x, ties.method = "average")` with NAs preserved; returns NA
#' when fewer than three finite values are present (as the engine does).
#'
#' @param v Numeric vector.
#' @return Numeric vector of 0-based ranks, or all NA.
#' @keywords internal
.rank_profile <- function(v) {
  out <- rep(NA_real_, length(v))
  ok <- is.finite(v)
  if (sum(ok) < 3L) return(out)
  out[ok] <- rank(v[ok], ties.method = "average") - 1
  out
}

#' Scale ranks to 0-1
#'
#' `.rank_profile()` divided by `n_finite - 1`, matching `_rank01()`.
#' @keywords internal
.rank01 <- function(v) {
  ok <- is.finite(v)
  n <- sum(ok)
  out <- rep(NA_real_, length(v))
  if (n < 1L) return(out)
  out[ok] <- .rank_profile(v)[ok] / max(1, n - 1)
  out
}

#' Pearson correlation on the finite pairs
#'
#' Same contract as the engine's `_pearson()`: fewer than three finite pairs
#' gives NaN, and a zero denominator gives NaN rather than an error. Returning
#' NaN (not 0) matters downstream: 0 would be a real similarity, NaN means "no
#' information" and the callers test for it explicitly.
#' @keywords internal
.pearson <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NaN)
  x <- a[ok] - mean(a[ok])
  y <- b[ok] - mean(b[ok])
  den <- sqrt(sum(x * x) * sum(y * y))
  if (!is.finite(den) || den <= 0) return(NaN)
  sum(x * y) / den
}

#' Spearman correlation, computed as Pearson on average ranks
#'
#' Deliberately not `stats::cor(method = "spearman")`: that route goes through
#' a different code path for ties and the last digit can drift from the engine,
#' which is fatal when the whole point is bit-for-bit parity.
#' @keywords internal
.spearman <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NaN)
  ra <- .rank_profile(a[ok])
  rb <- .rank_profile(b[ok])
  ra <- ra - mean(ra)
  rb <- rb - mean(rb)
  den <- sqrt(sum(ra * ra) * sum(rb * rb))
  if (!is.finite(den) || den <= 0) return(NaN)
  sum(ra * rb) / den
}

#' Decide whether a matrix is on the linear scale, then log2 it
#'
#' In-house matrices mix linear intensities/counts with already-log2 data.
#' Applying `log2(x + 1)` to log2 data silently compresses everything into
#' 0-5 and every threshold loses meaning, so the scale is decided from the
#' 99th percentile of finite values: above 50 means linear.
#'
#' @param x Numeric matrix.
#' @return A list with `x` (transformed matrix) and `scale`
#'   (`"linear"` when log2 was applied, `"log2"` when data were already log2).
#' @export
ery_as_log2 <- function(x) {
  v <- x[is.finite(x)]
  if (!length(v)) return(list(x = x, scale = "log2"))
  p99 <- as.numeric(stats::quantile(v, 0.99, type = 7, names = FALSE))
  if (p99 > 50) {
    return(list(x = log2(pmax(x, 0) + 1), scale = "linear"))
  }
  list(x = x, scale = "log2")
}

#' Metabolite name normalisation
#'
#' Lower-case, strip whitespace, hyphens and underscores, drop stereo prefixes.
#' Deliberately conservative: no synonym merging, because merging needs a
#' mapping table, not a guess.
#'
#' @param name Character vector of metabolite names.
#' @return Normalised character vector (`NA` for empty input).
#' @export
ery_norm_metabolite <- function(name) {
  s <- tolower(trimws(as.character(name)))
  s <- gsub("[[:space:]]+", "", s)
  s <- gsub("[-_]", "", s)
  s[!nzchar(s)] <- NA_character_
  s
}

## ---------------------------------------------------------------------------
## 输入校验
## ---------------------------------------------------------------------------

## 代谢物判定的缩写表：这些全大写短名是辅因子 / 核苷酸，不是基因
.met_acronym <- c("ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "NAD", "NADH",
                  "NADP", "NADPH", "FAD", "FADH2", "COA", "GSH", "GSSG",
                  "CAMP", "CGMP", "SAM", "SAH", "PEP", "DHAP", "GAP", "UDP",
                  "UDPG", "TCA")

.ewb_metabolite_lexicon <- function() {
  if (!is.null(.ewb_cache$met_lex)) return(.ewb_cache$met_lex)
  ## 由 data-raw/export_lexicon.py 从引擎的 SQLite 导出：base R 读不了 SQLite，
  ## 而"零依赖"是这个包的硬承诺，所以这份 183 个名字的名单只能冻结成静态资源
  p <- .ewb_extdata("metabolite_lexicon.txt")
  lex <- if (file.exists(p)) trimws(readLines(p, warn = FALSE, encoding = "UTF-8")) else character(0)
  lex <- unique(lex[nzchar(lex)])
  .ewb_cache$met_lex <- lex
  lex
}

.is_metabolite_name <- function(x, lex) {
  if (toupper(x) %in% .met_acronym) return(TRUE)
  nm <- ery_norm_metabolite(x)
  if (!is.na(nm) && nm %in% lex) return(TRUE)
  if (grepl("^(?:C[0-9]{5}|HMDB[0-9]+|LM[A-Z]{2}[0-9]+)", x,
            ignore.case = TRUE, perl = TRUE)) return(TRUE)
  if (grepl("^[0-9]+[-_][A-Za-z]", x, perl = TRUE)) return(TRUE)
  ## 带空格的长名（"Glucose 6 phosphate"）
  if (grepl("^[A-Za-z0-9][A-Za-z0-9\\-\\(\\)\\[\\]',. ]*\\s+[A-Za-z0-9-]+",
            x, perl = TRUE)) return(TRUE)
  ## 立体前缀（L-alanine）、脂质缩写（PC:34:1 / PE(36:2)）
  if (grepl("^[LDld]-[A-Za-z]", x, perl = TRUE)) return(TRUE)
  if (grepl("^[A-Z]{2,4}[:(]", x, perl = TRUE)) return(TRUE)
  ## 后缀规则只在"一个词"上生效：蛋白组的长描述（...chain OS=... SV=2）里
  ## 也含 -ine/-ate，按后缀判会把整份蛋白组误判成代谢组
  if (grepl(" ", x, fixed = TRUE) || grepl("=", x, fixed = TRUE)) return(FALSE)
  grepl("(?:ate|ine|osine|oside|ose|ol|acid|amine|amide|carnitine|phosphate|pyruvate|lactate)$",
        x, ignore.case = TRUE, perl = TRUE)
}

#' Diagnose the identifier space of a feature vector
#'
#' Decides whether the identifiers can be used at all and, when they cannot,
#' which kind they are: gene symbol, microRNA, metabolite, probe or unknown.
#' The verdict is deliberately **not** based on canonicalisation hit rates —
#' the metabolite branch of `ery_canonical()` accepts any non-probe string, so
#' its hit rate is always 100% and carries no information. This function looks
#' at the written form instead, and returns all four ratios so the caller can
#' explain what to do.
#'
#' The order matters and is the same as the engine's: probe, microRNA,
#' protein/gene identifier, metabolite, then bare symbol. Metabolites are tested
#' before symbols because `L-alanine` and `PC:34:1` also look like symbol-ish
#' strings, and because a proteomics long description can contain `GN=` *and*
#' an `-ine` suffix.
#'
#' @param features Character vector of feature identifiers.
#' @return A list with `verdict` (one of `gene`, `mirna`, `metabolite`, `probe`,
#'   `unknown`), `n`, `symbol_ratio` and the four per-class `ratios`.
#' @examples
#' ery_id_space(c("HBB", "HBA1", "SLC4A1"))$verdict
#' ery_id_space(c("L-alanine", "L-arginine"))$verdict
#' @export
ery_id_space <- function(features) {
  f <- trimws(as.character(features))
  f <- f[!is.na(f) & nzchar(f)]
  n <- length(f)
  if (!n) {
    return(list(verdict = "unknown", n = 0L, symbol_ratio = 0,
                ratios = list(probe = 0, mirna = 0, gene = 0, met = 0)))
  }
  lex <- .ewb_metabolite_lexicon()
  cnt <- c(probe = 0L, mirna = 0L, gene = 0L, met = 0L)
  for (x in f) {
    k <- NA_character_
    if (grepl("(?:_[a-z]{0,2}_?at$|^AFFX-|^A_\\d+_P\\d+|^ILMN_\\d+|^CUST_|^anti-sense|^chr[0-9vxy]|_(?:Exon|Intron)\\d*$|^\\d+$|^\\d+_\\d+$)",
              x, ignore.case = TRUE, perl = TRUE)) {
      k <- "probe"
    } else if (grepl("^(?:(?:hsa|mmu|rno|ath|cel|dre)-)?-?mir", x,
                     ignore.case = TRUE, perl = TRUE)) {
      k <- "mirna"
    } else if (grepl("^(?:sp|tr)\\|[A-Z0-9]+\\||\\bGN=[A-Za-z0-9_.-]+", x,
                     ignore.case = TRUE, perl = TRUE) ||
               grepl("^ENS[A-Z]*[0-9]{6,}", x, ignore.case = TRUE, perl = TRUE)) {
      k <- "gene"
    } else if (.is_metabolite_name(x, lex)) {
      k <- "met"
    } else if (grepl("^[A-Z][A-Z0-9.-]{1,19}$", x, perl = TRUE) ||
               grepl("^[A-Z][a-z][a-z0-9]{0,13}$", x, perl = TRUE)) {
      k <- "gene"
    }
    if (!is.na(k)) cnt[k] <- cnt[k] + 1L
  }
  ratios <- as.list(round(cnt / n, 3))
  verdict <- if (ratios$mirna >= 0.5) "mirna"
             else if (ratios$probe >= 0.5) "probe"
             else if (ratios$gene >= 0.6) "gene"
             else if (ratios$met >= 0.5) "metabolite"
             else if (ratios$gene >= 0.3) "gene"
             else "unknown"
  list(verdict = verdict, n = n, symbol_ratio = ratios$gene, ratios = ratios)
}

## ---------------------------------------------------------------------------
## 特征名标准化：canonical() + to_space()
##
## 与 backend/engines/disease_match.canonical 和 toolkit_core.to_space 逐条对齐。
## 不做这一层，别名输入（UniProt、带版本号的 Ensembl、旧符号）会被当成不同
## 特征，与网站算出的东西不是一回事。
## ---------------------------------------------------------------------------

.ewb_display_map <- function() {
  if (!is.null(.ewb_cache[["display_map"]])) return(.ewb_cache[["display_map"]])
  m <- tryCatch(readRDS(.ewb_extdata("gene_display_map.rds")),
                error = function(e) character(0))
  .ewb_cache[["display_map"]] <- m
  m
}

## 探针/序列 ID：这些不是基因符号，归不到就丢弃
.PROBE_RE <- "^(AFFX|ILMN_|[0-9]+_s?_at$|[0-9]+_at$|NM_|NR_|XM_|XR_|ENST|AMBI_)"
.SYMBOL_RE <- "^[A-Z][A-Z0-9.-]{0,14}$"
.MIRNA_RE <- "^hsa-(mir|let)-"
.UNI_RE <- "^(sp|tr)\\|([A-Z0-9]+)(-[0-9]+)?\\|"

#' Canonicalise feature identifiers for a given space
#'
#' Mirrors `disease_match.canonical()`: miRNA names are lower-cased and
#' validated, metabolite names pass through unless they look like probes, and
#' gene identifiers are resolved through UniProt `GN=`, `sp|ACC|` wrappers,
#' version suffixes and the shipped alias table before being checked against the
#' symbol pattern. Anything unresolvable becomes `NA` and is dropped.
#'
#' @param fids Character vector of identifiers.
#' @param space One of `"gene"`, `"mirna"`, `"metabolite"`.
#' @return Character vector, `NA` where unresolvable.
#' @export
ery_canonical <- function(fids, space = "gene") {
  s <- trimws(as.character(fids))
  s <- gsub('^"+|"+$', "", s)
  low <- tolower(s)
  bad <- !nzchar(s) | low %in% c("nan", "none", "unk")
  out <- rep(NA_character_, length(s))

  if (identical(space, "mirna")) {
    is_mir <- grepl(.MIRNA_RE, s, ignore.case = TRUE) & !bad
    out[is_mir] <- tolower(s[is_mir])
    return(out)
  }
  if (identical(space, "metabolite")) {
    is_probe <- grepl(.PROBE_RE, s, ignore.case = TRUE)
    out[!bad & !is_probe] <- s[!bad & !is_probe]
    return(out)
  }

  d <- .ewb_display_map()
  for (i in seq_along(s)) {
    if (bad[i]) next
    si <- s[i]
    ## 1) GN=SYMBOL（UniProt 头）
    m <- regmatches(si, regexpr("GN=[A-Za-z0-9_.-]+", si))
    if (length(m)) {
      g <- toupper(sub("^GN=", "", m))
      hit <- if (length(d)) d[g] else NA_character_
      out[i] <- if (!is.na(hit)) unname(hit) else
        if (grepl(.SYMBOL_RE, g)) g else NA_character_
      next
    }
    ## 2) sp|ACC|NAME / tr|ACC|NAME 外壳
    m2 <- regmatches(si, regexpr(.UNI_RE, si, ignore.case = TRUE))
    if (length(m2)) si <- toupper(sub(.UNI_RE, "\\2", m2, ignore.case = TRUE))
    u <- toupper(si)
    ## 3) 别名表（先原名，再去版本号）
    hit <- if (length(d)) d[u] else NA_character_
    if (!is.na(hit)) { out[i] <- unname(hit); next }
    base <- sub("\\.\\d+$", "", u)
    hit <- if (length(d)) d[base] else NA_character_
    if (!is.na(hit)) { out[i] <- unname(hit); next }
    ## 4) 探针丢弃，其余按符号正则校验
    if (grepl(.PROBE_RE, u, ignore.case = TRUE)) next
    if (grepl(.SYMBOL_RE, u)) out[i] <- u
  }
  out
}

#' Collapse a matrix onto a canonical feature space
#'
#' Same rule as the engine: when several rows map to the same canonical name,
#' keep the row with the **largest spread** (max - min). Never average —
#' averaging multi-probe rows cancels signals that point in opposite directions.
#'
#' @param x Numeric matrix, features in rows.
#' @param features Character vector of raw identifiers (defaults to row names).
#' @param space One of `"gene"`, `"mirna"`, `"metabolite"`.
#' @return A list with `x` (collapsed matrix) and `features` (canonical names).
#' @export
ery_to_space <- function(x, features = NULL, space = "gene") {
  if (is.data.frame(x)) x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(features)) features <- rownames(x)
  if (is.null(features)) stop("`x` needs row names.", call. = FALSE)
  nm <- ery_canonical(features, space)
  keep <- which(!is.na(nm))
  if (!length(keep)) {
    return(list(x = x[0, , drop = FALSE], features = character(0)))
  }
  nm <- nm[keep]
  xk <- x[keep, , drop = FALSE]
  ## 离散度：有限值的 max - min；全缺失记 -1
  spread <- apply(xk, 1, function(v) {
    v <- v[is.finite(v)]
    if (!length(v)) -1 else max(v) - min(v)
  })
  ord <- order(spread, decreasing = TRUE)
  first <- !duplicated(nm[ord])
  sel <- ord[first]
  out <- xk[sel, , drop = FALSE]
  rownames(out) <- nm[sel]
  ## 输出顺序按首次出现，便于与输入对照
  o <- order(sel)
  list(x = out[o, , drop = FALSE], features = nm[sel][o])
}

## 把上传矩阵整理成"特征 × 样本"的数值矩阵。
## 只做类型整理，不做名空间折叠——折叠交给 ery_to_space()，因为要先知道
## 数据属于哪个空间（基因 / miRNA / 代谢物）才能归一。
.clean_matrix <- function(x, features = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("`x` must be a matrix or data.frame.", call. = FALSE)
  storage.mode(x) <- "double"
  if (is.null(features)) {
    if (is.null(rownames(x))) stop("`x` needs rownames (feature identifiers).", call. = FALSE)
    features <- rownames(x)
  }
  list(x = x, features = as.character(features))
}
