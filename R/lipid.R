## ---------------------------------------------------------------------------
## 红细胞膜脂质成熟轨迹：网织红（0）→ 衰老红细胞（4）
##
## 移植自 app/backend/engines/age_axis.py。四步，每一步都有它非如此不可的理由：
##
##   1. 解析脂质名。库里同时存在 `PC:34:1 (16:0-18:1w9)`、`PC(34:1)`、
##      `SM(t42:3)`、`TG(18:2_17:1_18:2)` 四种写法；字符串等值匹配一条都对不上，
##      而"对不上"在页面会显示成"你的数据里没有这些脂质"——那是误导。
##   2. 按 **A 级（类别+总碳+总双键+醚键+碱基全同）/ B 级（只对类别+碳+双键）**
##      把输入行对到参考的 45 个键上。
##   3. **在样本实际测到的那个子集上重建曲线**，不是拿全 45 个的均值/标准差硬套。
##      占比依赖分母：样本只测到 17 个分子种时，这 17 个占它自己的 100%，而参考里
##      只占约 40%，两者的 z 不在同一个坐标系。
##   4. 对 81 个网格点一次算完相关系数，取峰值位置；再做置换零分布、逐脂质位置
##      一致率与二项检验。
## ---------------------------------------------------------------------------

.ery_grid <- seq(0, 4, length.out = 81)

## 命中脂质少于这个数就不投影：12 个分子种已低于"一类三四个"，谱形相关不可信
.ery_min_lipids <- 12L
## 谱形相关门槛：参考自投影 1.00，留一阶段 0.90+；跨平台实测只有 0.22–0.29，
## 所以 0.45 能把自己人放进去、把对不上的挡在外面
.ery_min_rho <- 0.45
## 判定"区间"时取 r 不低于峰值 90% 的位置范围
.ery_band_frac <- 0.90

.ery_lipid_alias <- c(
  SPHINGOMYELIN = "SM", CERAMIDE = "CER", TRIACYLGLYCEROL = "TG",
  DIACYLGLYCEROL = "DG", CHOLESTEROL = "CHOL",
  PHOSPHATIDYLCHOLINE = "PC", PHOSPHATIDYLETHANOLAMINE = "PE",
  PHOSPHATIDYLSERINE = "PS", PHOSPHATIDYLINOSITOL = "PI",
  LYSOPHOSPHATIDYLCHOLINE = "LPC", LYSOPHOSPHATIDYLETHANOLAMINE = "LPE"
)

.ery_zones <- list(
  list(lo = 0.0, hi = 0.5, key = "retic_young", zh = "\u7f51\u7ec7\u7ea2\uff08\u5e74\u8f7b\uff09", en = "reticulocyte, young"),
  list(lo = 0.5, hi = 1.5, key = "retic_mature", zh = "\u7f51\u7ec7\u7ea2\uff08\u6210\u719f\uff09", en = "reticulocyte, mature"),
  list(lo = 1.5, hi = 2.5, key = "ery_young", zh = "\u7ea2\u7ec6\u80de\uff08\u5e74\u8f7b\uff09", en = "erythrocyte, young"),
  list(lo = 2.5, hi = 3.5, key = "ery_mature", zh = "\u7ea2\u7ec6\u80de\uff08\u6210\u719f\uff09", en = "erythrocyte, mature"),
  list(lo = 3.5, hi = 4.01, key = "ery_old", zh = "\u7ea2\u7ec6\u80de\uff08\u8870\u8001\uff09", en = "erythrocyte, old")
)

.ery_zone <- function(pos) {
  for (z in .ery_zones) {
    if (pos >= z$lo && pos < z$hi) return(z[c("key", "zh", "en")])
  }
  z <- .ery_zones[[length(.ery_zones)]]
  z[c("key", "zh", "en")]
}

#' Parse a lipid name into (class, total carbons, total double bonds, ether, base)
#'
#' Handles the four naming conventions that coexist in public lipidomics data:
#' `PC:34:1 (16:0-18:1w9)`, `PC(34:1)`, `SM(t42:3)` and `TG(18:2_17:1_18:2)`.
#' Exact string matching fails on all four, and the failure is silent — the
#' sample simply looks like it has no lipids.
#'
#' @param name Character vector of lipid names.
#' @return A data.frame with columns `name`, `cls`, `c`, `db`, `ether`, `base`;
#'   rows that cannot be parsed are dropped.
#' @examples
#' parse_lipid(c("PC:34:1 (16:0-18:1w9)", "SM(d42:1)", "TG(18:2_17:1_18:2)"))
#' @export
parse_lipid <- function(name) {
  one <- function(nm) {
    s <- trimws(as.character(nm)[1])
    if (!nzchar(s) || is.na(s)) return(NULL)
    mh <- regexpr("^[A-Za-z-]+", s, perl = TRUE)
    if (mh[1] != 1L) return(NULL)
    head <- regmatches(s, mh)
    ## 注意别用 [[ ]]：命名向量取不存在的名字会抛错（不是返回 NULL）
    cls <- unname(.ery_lipid_alias[toupper(head)])
    if (is.na(cls)) cls <- toupper(head)
    body <- trimws(substring(s, attr(mh, "match.length") + 1L))
    if (!nzchar(body)) return(NULL)

    ether <- FALSE
    base <- ""
    total <- NULL

    ## 写法一：34:1 / (34:1) / :34:1 / (d34:1) / (34:1e)
    m1 <- regexec(paste0("^[(\\s:-]*(?:([dtm])-?)?\\s*(?:O-|P-)?\\s*",
                         "(\\d+)\\s*:\\s*(\\d+)\\s*([ep])?\\s*\\)?"),
                  body, ignore.case = TRUE, perl = TRUE)
    g <- regmatches(body, m1)[[1]]
    if (length(g) >= 4L) {
      base <- tolower(g[2])
      total <- c(as.numeric(g[3]), as.numeric(g[4]))
      if (length(g) >= 5L && nzchar(g[5])) ether <- TRUE
    }

    ## 括号里的链式写法：(16:0-18:1) / (16:0_18:1) / (d18:1-16:0)
    ms <- regexpr("\\(([^()]*)\\)", body, perl = TRUE)
    spec <- if (ms[1] != -1L) trimws(gsub("^\\(|\\)$", "", regmatches(body, ms))) else ""
    if (nzchar(spec)) {
      mb <- regmatches(spec, regexec("^([dtm])(?=[0-9])", spec,
                                     ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(mb) >= 2L) base <- tolower(mb[2])
      ## 只在**没有给出总式**时才求和，判据是括号之前有没有「碳:双键」：
      ##   PC:34:1 (16:0-18:1w9) → 括号前的 :34:1 就是总式，必须尊重
      ##   PC:38:5 (18:0-20:4w6) → 数据写 38:5、链注释只有 38:4，求和等于改写原始数据
      ##   TG(18:2_17:1_18:2)    → 括号前没有总式，求和才是这个分子种的总式
      segs <- regmatches(spec, gregexpr("(?:[dtm]-?)?(?:O-|P-)?(\\d+)\\s*:\\s*(\\d+)",
                                       spec, ignore.case = TRUE, perl = TRUE))[[1]]
      pre_curly <- sub("\\(.*$", "", body)
      has_total <- grepl("\\d+\\s*:\\s*\\d+", pre_curly, perl = TRUE)
      if (length(segs) >= 2L && !has_total) {
        nums <- lapply(segs, function(g) {
          regmatches(g, regexec("(\\d+)\\s*:\\s*(\\d+)", g, perl = TRUE))[[1]]
        })
        nums <- Filter(function(v) length(v) >= 3L, nums)
        if (length(nums)) {
          total <- c(sum(vapply(nums, function(v) as.numeric(v[2]), 0)),
                     sum(vapply(nums, function(v) as.numeric(v[3]), 0)))
        }
      }
    }

    if (is.null(total)) return(NULL)
    ## 醚键三种写法：O-34:1 / P-36:2 前缀，以及 PC(34:1)e 这种**后缀**
    ## （后缀里 e 在括号之后，原来的 [ep]\s*\) 匹配不到，plasmanyl 被当成二酰基）
    if (grepl("\\(\\s*[OP]-", body, ignore.case = TRUE, perl = TRUE) ||
        grepl("[ep]\\s*\\)?\\s*$", body, ignore.case = TRUE, perl = TRUE)) ether <- TRUE
    if (total[1] <= 0) return(NULL)
    list(cls = cls, c = total[1], db = total[2], ether = ether, base = base)
  }
  out <- lapply(name, one)
  keep <- !vapply(out, is.null, TRUE)
  if (!any(keep)) {
    return(data.frame(name = character(0), cls = character(0), c = numeric(0),
                      db = numeric(0), ether = logical(0), base = character(0),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    name = as.character(name)[keep],
    cls = vapply(out[keep], `[[`, "", "cls"),
    c = vapply(out[keep], `[[`, 0, "c"),
    db = vapply(out[keep], `[[`, 0, "db"),
    ether = vapply(out[keep], `[[`, TRUE, "ether"),
    base = vapply(out[keep], `[[`, "", "base"),
    stringsAsFactors = FALSE
  )
}

.ery_lipid_key_strict <- function(cls, c, db, ether, base) {
  paste(cls, c, db, ifelse(ether, "T", "F"), base, sep = "\r")
}

.ery_lipid_key_loose <- function(cls, c, db) paste(cls, c, db, sep = "\r")

#' Fit the reference maturation trajectory from the in-house lipid panel
#'
#' @return A list with `keys`, `names`, `stages`, `raw` (anchor x stage means),
#'   `pct` (closed to 100%), `mu`, `sd`, `z`, `ok`, `curve` (anchor x 81 grid)
#'   and `stage_meta`.
#' @keywords internal
ery_lipid_model <- function() {
  if (!is.null(.ewb_cache$lipid_model)) return(.ewb_cache$lipid_model)
  panel <- .ewb_read_rds("lipid_panel.rds")
  raw_mat <- panel$raw
  stages <- panel$stages

  ## 按 A 级键聚合成锚点：同一键的多行取均值（SM:d42:1 有两种链式写法）
  keys_chr <- character(0)
  key_meta <- list()
  bucket <- list()
  name_of <- character(0)
  for (i in seq_len(nrow(raw_mat))) {
    nm <- rownames(raw_mat)[i]
    p <- parse_lipid(nm)
    if (!nrow(p)) next
    k <- .ery_lipid_key_strict(p$cls, p$c, p$db, p$ether, p$base)
    if (is.null(bucket[[k]])) {
      bucket[[k]] <- numeric(0)
      keys_chr <- c(keys_chr, k)
      key_meta[[length(key_meta) + 1L]] <- list(
        cls = p$cls, c = p$c, db = p$db, ether = p$ether, base = p$base)
      name_of[k] <- nm
    }
    bucket[[k]] <- c(bucket[[k]], i)
  }
  A <- t(vapply(keys_chr, function(k) {
    colMeans(raw_mat[bucket[[k]], , drop = FALSE])
  }, numeric(length(stages))))
  rownames(A) <- keys_chr

  ## 参考数据是组成型（每列总和归到 300），五个阶段分母几乎恒等，但只取 45 个
  ## 分子种时总和会飘（104.37 → 100.73）。先闭合到 100% 再标准化，用户数据也
  ## 一样闭合，"占比"这一层才是同一种量。
  colsum <- colSums(A)
  P <- sweep(A, 2L, ifelse(colsum > 0, colsum, NA_real_), "/") * 100
  mu <- rowMeans(P)
  sd <- apply(P, 1L, stats::sd)
  sd[!is.finite(sd) | sd < 1e-9] <- NA_real_
  Z <- sweep(sweep(P, 1L, mu, "-"), 1L, sd, "/")
  ok <- apply(Z, 1L, function(v) all(is.finite(v)))

  S <- seq_along(stages) - 1
  curve <- matrix(NA_real_, nrow = length(keys_chr), ncol = length(.ery_grid))
  for (i in seq_along(keys_chr)) {
    curve[i, ] <- if (all(is.finite(Z[i, ]))) {
      stats::approx(S, Z[i, ], xout = .ery_grid, rule = 2)$y
    } else NA_real_
  }

  ## 不要用 strsplit 反解键串：base 为空时末尾的空字段会被丢掉（只剩 4 段），
  ## rbind 出来列数就错位，后面 keys$c/keys$db 全是错的值。这里直接记录字段。
  key_df <- data.frame(
    cls = vapply(key_meta, `[[`, "", "cls"),
    c = vapply(key_meta, `[[`, 0, "c"),
    db = vapply(key_meta, `[[`, 0, "db"),
    ether = vapply(key_meta, `[[`, TRUE, "ether"),
    base = vapply(key_meta, `[[`, "", "base"),
    stringsAsFactors = FALSE)

  out <- list(
    keys = key_df, key_chr = keys_chr,
    names = unname(name_of[keys_chr]),
    stages = stages, raw = A, pct = P, mu = mu, sd = sd, z = Z,
    ok = ok, curve = curve, grid = .ery_grid,
    stage_meta = panel$stage_meta
  )
  .ewb_cache$lipid_model <- out
  out
}

#' Leave-one-stage-out reconstruction error
#'
#' The only honest accuracy figure for this module: self-projection is 1.00 by
#' construction (the curve is built from those five points) and must not be
#' quoted as accuracy. Here each stage is dropped, the curve rebuilt from the
#' other four, and the held-out stage projected back.
#'
#' @return A list with `loso_errors`, `loso_mae`, `loso_max`, `loso_mae_pct`
#'   and `n_stages`.
#' @keywords internal
ery_lipid_validation <- function() {
  m <- ery_lipid_model()
  ok <- m$ok
  n_stage <- length(m$stages)
  S <- seq_len(n_stage) - 1
  errs <- numeric(n_stage)
  for (h in seq_len(n_stage)) {
    cols <- setdiff(seq_len(n_stage), h)
    P <- m$pct[ok, cols, drop = FALSE]
    cs <- colSums(P)
    P <- sweep(P, 2L, ifelse(cs > 0, cs, NA_real_), "/") * 100
    mu <- rowMeans(P)
    sd <- apply(P, 1L, stats::sd)
    sd[!is.finite(sd) | sd < 1e-9] <- NA_real_
    Zh <- sweep(sweep(P, 1L, mu, "-"), 1L, sd, "/")
    cur <- matrix(NA_real_, nrow(Zh), length(.ery_grid))
    for (i in seq_len(nrow(Zh))) {
      cur[i, ] <- if (all(is.finite(Zh[i, ]))) {
        stats::approx(S[cols], Zh[i, ], xout = .ery_grid, rule = 2)$y
      } else NA_real_
    }
    prof <- m$z[ok, h]
    r <- vapply(seq_along(.ery_grid), function(g) .pearson(prof, cur[, g]), 0)
    r[!is.finite(r)] <- -1
    errs[h] <- .ery_grid[which.max(r)] - S[h]
  }
  ae <- abs(errs)
  list(loso_errors = round(errs, 3), loso_mae = round(mean(ae), 3),
       loso_max = round(max(ae), 3),
       loso_mae_pct = round(mean(ae) / 4 * 100, 1), n_stages = n_stage)
}

## 一个样本对 81 个网格点的相关系数，一次算完。逐个网格点调 cor() 在置换检验
## 里会被放大约 200 倍，慢到不可用。
.corr_to_grid <- function(prof, C) {
  p <- prof - mean(prof)
  ## 按**列**中心化：每个网格点是一列，减的是这一列 45 个锚点的均值。
  ## 三处都容易写错，任一处错都会让 |rho| > 1：
  ##   sweep 的 MARGIN 必须是 2（列），STATS 用 colMeans；
  ##   den 用 colSums（每列平方和），不是 rowSums。
  Cm <- sweep(C, 2L, colMeans(C), "-")
  num <- as.numeric(p %*% Cm)
  den <- sqrt(sum(p * p) * colSums(Cm * Cm))
  out <- rep(NaN, ncol(C))
  okd <- is.finite(den) & den > 0
  out[okd] <- num[okd] / den[okd]
  out
}

## 把锚点限制到样本真正测到的那些脂质上，重新闭合 → 标准化 → 插值
.lipid_subset_model <- function(idx) {
  m <- ery_lipid_model()
  P <- m$raw[idx, , drop = FALSE]
  cs <- colSums(P)
  P <- sweep(P, 2L, ifelse(cs > 0, cs, NA_real_), "/") * 100
  mu <- rowMeans(P)
  sd <- apply(P, 1L, stats::sd)
  ok <- is.finite(sd) & sd > 1e-9
  Z <- matrix(NaN, nrow(P), ncol(P))
  if (any(ok)) {
    Z[ok, ] <- sweep(sweep(P[ok, , drop = FALSE], 1L, mu[ok], "-"), 1L, sd[ok], "/")
  }
  S <- seq_len(ncol(P)) - 1
  C <- matrix(NaN, nrow(P), length(.ery_grid))
  for (i in seq_len(nrow(P))) {
    if (ok[i]) C[i, ] <- stats::approx(S, Z[i, ], xout = .ery_grid, rule = 2)$y
  }
  list(mu = mu, sd = sd, ok = ok, C = C)
}

## 把输入行名对到参考的键上：A 级严格，B 级忽略醚键 / 碱基
.lipid_annotate_rows <- function(features) {
  m <- ery_lipid_model()
  strict <- stats::setNames(seq_along(m$key_chr), m$key_chr)
  loose <- list()
  for (i in seq_along(m$key_chr)) {
    kl <- .ery_lipid_key_loose(m$keys$cls[i], m$keys$c[i], m$keys$db[i])
    if (is.null(loose[[kl]])) loose[[kl]] <- i
  }
  out <- list()
  for (nm in features) {
    p <- parse_lipid(nm)
    if (!nrow(p)) next
    ks <- .ery_lipid_key_strict(p$cls, p$c, p$db, p$ether, p$base)
    kl <- .ery_lipid_key_loose(p$cls, p$c, p$db)
    if (!is.na(strict[ks]) && m$ok[strict[[ks]]]) {
      i <- strict[[ks]]
      out[[length(out) + 1L]] <- list(row = i, name = as.character(nm),
                                      anchor = m$names[i], tier = "A")
    } else if (!is.null(loose[[kl]]) && m$ok[loose[[kl]]]) {
      i <- loose[[kl]]
      out[[length(out) + 1L]] <- list(row = i, name = as.character(nm),
                                      anchor = m$names[i], tier = "B")
    }
  }
  out
}

## 谱形相关的置换零分布：打乱样本在各脂质上的取值，看最大相关能到多少。
## 没有这个数，"ρ=0.52" 在 12 个脂质上到底是强信号还是噪声就没法说。
.perm_p_lipid <- function(prof, C, observed, n_perm = 200L, seed = 20261112L,
                          permutations = NULL) {
  hits <- 0L
  if (!is.null(permutations)) {
    for (t in seq_len(nrow(permutations))) {
      r <- .corr_to_grid(prof[permutations[t, ]], C)
      if (is.finite(max(r, na.rm = TRUE)) && max(r, na.rm = TRUE) >= observed) hits <- hits + 1L
    }
    return((hits + 1) / (nrow(permutations) + 1))
  }
  ## 保存并恢复 .Random.seed：置换检验不该改变调用方后续的随机流
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else NULL
  set.seed(seed)
  n <- length(prof)
  for (t in seq_len(n_perm)) {
    r <- .corr_to_grid(prof[sample.int(n)], C)
    if (is.finite(max(r, na.rm = TRUE)) && max(r, na.rm = TRUE) >= observed) hits <- hits + 1L
  }
  if (is.null(old_seed)) {
    rm(".Random.seed", envir = globalenv())
  } else {
    assign(".Random.seed", old_seed, envir = globalenv())
  }
  (hits + 1) / (n_perm + 1)
}

## P(X >= k)，X ~ Bin(n, p)。n <= 200 时直接累加，别为这点事引依赖
.binom_sf <- function(k, n, p) {
  if (n <= 0) return(1)
  sum(vapply(k:n, function(i) choose(n, i) * p^i * (1 - p)^(n - i), 0))
}

## 把一个样本投到轨迹上，返回 0–4 的阶段坐标
.lipid_project_sample <- function(mat, col, rows, permutations = NULL) {
  m <- ery_lipid_model()
  values <- mat[, col]
  by_row <- list(); name_of <- character(0); tier_of <- character(0)
  for (r in rows) {
    v <- values[match(r$name, rownames(mat))]
    if (!length(v) || !is.finite(v)) next
    k <- as.character(r$row)
    by_row[[k]] <- c(by_row[[k]], v)
    if (!nzchar(name_of[k]) || is.na(name_of[k])) name_of[k] <- r$name
    if (!nzchar(tier_of[k]) || is.na(tier_of[k])) tier_of[k] <- r$tier
  }
  idx <- sort(as.integer(names(by_row)))
  base <- list(sample = col, n_matched = length(idx), n_anchor = sum(m$ok))

  if (length(idx) < .ery_min_lipids) {
    return(c(base, list(ok = FALSE, reason = "insufficient_overlap",
                        detail = list(threshold = .ery_min_lipids))))
  }
  x <- vapply(idx, function(i) mean(by_row[[as.character(i)]]), 0)
  if (!all(is.finite(x))) {
    return(c(base, list(ok = FALSE, reason = "non_numeric", detail = list())))
  }
  if (any(x <= 0)) {
    ## 组成型闭合要求正值；log2 强度、log 比值都会跌破 0，硬闭合会造出一个
    ## 没有意义的比例，宁可退回去让用户换一份线性尺度（强度 / 百分比）
    return(c(base, list(ok = FALSE, reason = "requires_positive_values",
                        detail = list(min_value = round(min(x), 3)))))
  }

  sm <- .lipid_subset_model(idx)
  keep <- sm$ok & is.finite(x)
  if (sum(keep) < .ery_min_lipids) {
    return(c(base, list(ok = FALSE, reason = "insufficient_varying_lipids",
                        detail = list(threshold = .ery_min_lipids,
                                      n_usable = sum(keep)))))
  }
  pct <- x / sum(x) * 100
  prof <- as.numeric(sweep(sweep(matrix(pct, ncol = 1L), 1L, sm$mu, "-"),
                           1L, sm$sd, "/"))
  rs <- .corr_to_grid(prof[keep], sm$C[keep, , drop = FALSE])
  rs[!is.finite(rs)] <- -1
  best <- which.max(rs)
  pos <- .ery_grid[best]
  rho <- rs[best]
  p_perm <- .perm_p_lipid(prof[keep], sm$C[keep, , drop = FALSE], rho,
                          permutations = permutations)

  detail <- list(rho = round(rho, 3), p_perm = round(p_perm, 4),
                 threshold_rho = .ery_min_rho, n_used = sum(keep))
  if (rho < .ery_min_rho || p_perm > 0.01) {
    return(c(base, list(ok = FALSE, reason = "off_trajectory", detail = detail)))
  }

  band <- .ery_grid[rs >= .ery_band_frac * rho]
  lo <- if (length(band)) min(band) else pos
  hi <- if (length(band)) max(band) else pos

  ## 逐脂质位置：这个分子种的取值落在它自己轨迹的哪一段
  zall <- prof[keep]
  Ck <- sm$C[keep, , drop = FALSE]
  per <- vapply(seq_len(nrow(Ck)), function(k) {
    .ery_grid[which.min(abs(Ck[k, ] - zall[k]))]
  }, 0)
  within <- mean(abs(per - pos) <= 0.5)
  p_binom <- .binom_sf(as.integer(round(within * nrow(Ck))), nrow(Ck), 0.25)

  zone <- .ery_zone(pos)
  anchors <- m$stage_meta
  nearest <- anchors[which.min(abs(anchors$position - pos)), ]
  bracket <- NULL
  for (i in seq_len(nrow(anchors) - 1L)) {
    if (anchors$position[i] <= pos && pos <= anchors$position[i + 1L]) {
      bracket <- list(from = as.list(anchors[i, ]), to = as.list(anchors[i + 1L, ]))
      break
    }
  }

  idx_keep <- idx[keep]
  lip_rows <- lapply(seq_along(idx_keep), function(k) {
    i <- idx_keep[k]
    list(lipid = name_of[as.character(i)], anchor = m$names[i],
         tier = tier_of[as.character(i)],
         pct = round(pct[match(i, idx)], 3), z = round(zall[k], 2),
         position = round(per[k], 2))
  })
  if (length(lip_rows)) {
    ord <- order(-abs(vapply(lip_rows, `[[`, 0, "z")))
    lip_rows <- lip_rows[ord]
  }

  c(base, list(
    ok = TRUE, position = round(pos, 2), rho = round(rho, 3),
    p_perm = round(p_perm, 4), band = c(round(lo, 2), round(hi, 2)),
    zone = zone,
    nearest = list(key = nearest$key, position = nearest$position),
    nearest_label = list(zh = nearest$zh, en = nearest$en),
    bracket = bracket,
    n_within = round(within, 3), p_value = round(p_binom, 4),
    lipids = lip_rows,
    tiers = list(A = sum(tier_of[as.character(idx_keep)] == "A", na.rm = TRUE),
                 B = sum(tier_of[as.character(idx_keep)] == "B", na.rm = TRUE))
  ))
}

## 整份矩阵：逐样本投影 + 汇总 + 弃权原因统计
.lipid_project_matrix <- function(mat, permutations = NULL) {
  m <- ery_lipid_model()
  rows <- .lipid_annotate_rows(rownames(mat))
  if (!length(rows)) {
    return(list(ok = FALSE, reason = "no_lipid_names",
                n_anchor = sum(m$ok),
                examples_anchor = m$names[seq_len(min(5L, length(m$names)))],
                examples_input = utils::head(rownames(mat), 5)))
  }
  samples <- lapply(colnames(mat), function(cl) {
    .lipid_project_sample(mat, cl, rows, permutations)
  })
  good <- Filter(function(s) isTRUE(s$ok), samples)
  summary <- list(n_samples = length(samples), n_ok = length(good),
                  n_refused = length(samples) - length(good))
  if (length(good)) {
    pos <- vapply(good, `[[`, 0, "position")
    med <- stats::median(pos)
    summary <- c(summary, list(
      median_position = round(med, 2), min_position = round(min(pos), 2),
      max_position = round(max(pos), 2), median_zone = .ery_zone(med),
      n_young = sum(pos < 1.5), n_old = sum(pos >= 2.5),
      n_between = sum(pos >= 1.5 & pos < 2.5)))
  }
  reasons <- table(vapply(Filter(function(s) !isTRUE(s$ok), samples),
                          function(s) as.character(s$reason), ""))
  all_matched <- unique(vapply(rows, `[[`, 0, "row"))
  list(
    ok = length(good) > 0,
    model = list(dataset = "PMC11868425_lipidomics", n_lipids = sum(m$ok),
                 stages = m$stage_meta, min_lipids = .ery_min_lipids,
                 min_rho = .ery_min_rho),
    coverage = list(n_matched = length(all_matched), n_anchor = sum(m$ok),
                    tier_A = sum(vapply(rows, `[[`, "", "tier") == "A"),
                    tier_B = sum(vapply(rows, `[[`, "", "tier") == "B"),
                    classes = sort(unique(m$keys$cls[all_matched]))),
    validation = ery_lipid_validation(),
    samples = samples, summary = summary,
    refusal_reasons = as.list(reasons)
  )
}

## 便宜版的可用性检查：只算命中数与最好的谱形相关，不做置换检验
ery_lipid_probe <- function(mat) {
  rows <- .lipid_annotate_rows(rownames(mat))
  if (!length(rows)) return(list(n_matched = 0, max_rho = NULL, classes = character(0)))
  best <- NULL
  for (cl in colnames(mat)) {
    values <- mat[, cl]
    by_row <- list()
    for (r in rows) {
      v <- values[match(r$name, rownames(mat))]
      if (length(v) && is.finite(v)) by_row[[as.character(r$row)]] <- c(by_row[[as.character(r$row)]], v)
    }
    idx <- sort(as.integer(names(by_row)))
    if (length(idx) < .ery_min_lipids) next
    x <- vapply(idx, function(i) mean(by_row[[as.character(i)]]), 0)
    if (any(x <= 0)) next
    sm <- .lipid_subset_model(idx)
    keep <- sm$ok & is.finite(x)
    if (sum(keep) < .ery_min_lipids) next
    pct <- x / sum(x) * 100
    prof <- as.numeric(sweep(sweep(matrix(pct, ncol = 1L), 1L, sm$mu, "-"),
                             1L, sm$sd, "/"))
    rs <- .corr_to_grid(prof[keep], sm$C[keep, , drop = FALSE])
    rs[!is.finite(rs)] <- -1
    r <- max(rs)
    best <- if (is.null(best)) r else max(best, r)
  }
  mi <- unique(vapply(rows, `[[`, 0, "row"))
  m <- ery_lipid_model()
  list(n_matched = length(mi),
       max_rho = if (is.null(best)) NULL else round(best, 3),
       classes = sort(unique(m$keys$cls[mi])))
}
