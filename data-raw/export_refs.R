## export_refs.R — 把网站引擎的参照数据冻结成 RDS，随包分发。
##
## 为什么不用 Python 直接写 RDS：RDS 是 R 的序列化格式，手写要复刻
## 类型标签与属性，出错时不会报错只会静默读成错的东西。这里用 R 自己
## 读 JSON 再写 RDS，序列化交给 R。
##
## 为什么 simplifyVector = FALSE：toolkit_refs.json 的 purity / attribution
## 是深嵌套结构（templates、cell_prior、centroids），jsonlite 默认会把
## 同构数组简化成 data.frame，NULL 变成 NA，取字段时路径全变。保持嵌套
## list 后，R 侧访问路径与 Python 侧完全一致（$purity$templates$genes），
## 对拍时不会因为结构差异产生假阳性。
##
## 用法：在包根目录执行
##   Rscript data-raw/export_refs.R
## 源目录默认 ../../../app/backend/data（相对本脚本），可用环境变量
## ERYWB_APP_DATA 覆盖。

suppressPackageStartupMessages({
  library(jsonlite)
})

## 路径一律用相对路径：R 在 Windows 上遇到含非 ASCII 的绝对路径时，
## normalizePath() / dir.exists() 会抛 "file name conversion problem"，
## 而本仓库的上级目录名带中文。相对路径不经过那次转换，能正常工作。
## 因此本脚本必须在包根目录（含 DESCRIPTION 的那一层）执行。
stopifnot(file.exists("DESCRIPTION"))

## 计算原语（.num_vec / .rank_profile）直接复用包内实现，避免两套定义分叉。
source("R/utils.R")

pkg <- "."
src <- Sys.getenv("ERYWB_APP_DATA", unset = file.path("..", "..", "app", "backend", "data"))
out <- file.path("inst", "extdata")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(src)) {
  stop("找不到源数据目录（相对路径 ", src, "）。请在包根目录执行，",
       "或用 ERYWB_APP_DATA 指定 app/backend/data 的位置。", call. = FALSE)
}
message("source: ", src)
message("output: ", out)

sha256 <- function(p) {
  if (!file.exists(p)) return(NA_character_)
  unname(tools::md5sum(p))
}

json_ref <- function(name) {
  p <- file.path(src, name)
  if (!file.exists(p)) {
    warning("missing: ", name, call. = FALSE)
    return(NULL)
  }
  fromJSON(p, simplifyVector = FALSE)
}

## ---- 1. 层级校准：marker 名单、每层秩共识、p50/p95、端点切点 ----------
layers <- json_ref("layer_references.json")
if (!is.null(layers)) {
  ## 给每个层加 ASCII key。R 侧不能用中文层名做 list 索引：在非 UTF-8
  ## locale 的 Windows 上，源码里的中文字面量与数据里的 UTF-8 名字比较
  ## 会让 R 直接段错误（不是报错，是崩溃），这一层映射把它挡在数据侧。
  ## 顺序由 JSON 保持，同时把原名记进 key_label 供显示。
  nms <- names(layers$layers)
  keys <- c("transcriptome", "proteome")
  if (length(nms) != length(keys)) {
    warning("unexpected number of omics layers: ", length(nms),
            " (expected ", length(keys), "); keys may be misassigned", call. = FALSE)
  }
  for (i in seq_along(nms)) {
    if (i <= length(keys)) {
      layers$layers[[i]]$key <- keys[i]
      layers$layers[[i]]$key_label <- nms[i]
    }
  }
  saveRDS(layers, file.path(out, "layer_references.rds"), compress = "xz")
  message(sprintf("  layer_references.rds   layers=%d (%s) markers=%d",
                  length(layers$layers), paste(keys, collapse = "/"),
                  length(layers$markers$all)))
}

## ---- 2. 固定参考图谱：381 样本 × 47 marker 的秩/值向量 + OOD 闸门 -----
atlas <- json_ref("reference_atlas.json")
if (!is.null(atlas)) {
  ## 381 × 47 的长表更好用，但秩向量按 marker 顺序对齐才安全，所以保留
  ## 原结构，另外附一个矩阵视图供 kNN 用。
  ## 必须逐元素取：unlist() 会把 null 整个丢掉，后面的值左移、末尾补 NA，
  ## 缺失总数看起来一样但位置错了（实测 RBC1 本来是第 3 位缺失，unlist 后
  ## 变成第 24 位起连续缺失）。
  mk <- unlist(atlas$markers, use.names = FALSE)
  rank_mat <- t(vapply(atlas$samples, function(s) .num_vec(s$rank), numeric(length(mk))))
  colnames(rank_mat) <- mk
  meta <- data.frame(
    dataset_id = vapply(atlas$samples, function(s) as.character(s$dataset_id), ""),
    sample     = vapply(atlas$samples, function(s) as.character(s$sample), ""),
    layer      = vapply(atlas$samples, function(s) as.character(s$layer), ""),
    condition  = vapply(atlas$samples, function(s) as.character(s$condition), ""),
    compartment= vapply(atlas$samples, function(s) as.character(s$compartment), ""),
    labeled    = vapply(atlas$samples, function(s) isTRUE(s$labeled), logical(1)),
    stringsAsFactors = FALSE
  )
  stopifnot(nrow(meta) == nrow(rank_mat))

  ## rank 空间的投票 margin 门槛：留一数据集时跨数据集投票 margin 的 p10。
  ## 这个量与图谱固定绑定，不必每次调用重算（引擎里是 lru_cache 缓存的同一
  ## 个量），导出来更省事，也让「门槛从哪来」可被审稿人直接核对。
  ##
  ## 相似度必须**逐对在共同非缺失子集上重排秩**再求相关，不能直接用
  ## cor(..., use="pairwise.complete.obs")。图谱里 373/381 个样本有缺失
  ## marker（最多缺 31 个），子集重排会改变相对秩；省掉这一步算出的门槛是
  ## 0.154，而引擎是 0.433，差三倍——投票因此会在该弃权时给标签。
  sim_pair <- function(a, b) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 8L) return(NA_real_)
    x <- .rank_profile(a[ok])
    y <- .rank_profile(b[ok])
    sdx <- sqrt(mean((x - mean(x))^2))
    sdy <- sqrt(mean((y - mean(y))^2))
    if (sdx <= 0 || sdy <= 0) return(NA_real_)
    mean((x - mean(x)) * (y - mean(y))) / (sdx * sdy)
  }
  K <- 7L
  n_ref <- nrow(rank_mat)
  S <- matrix(NA_real_, n_ref, n_ref)
  rows <- lapply(seq_len(n_ref), function(i) rank_mat[i, ])
  for (i in seq_len(n_ref - 1L)) {
    for (j in (i + 1L):n_ref) {
      v <- sim_pair(rows[[i]], rows[[j]])
      S[i, j] <- v
      S[j, i] <- v
    }
  }
  ds <- meta$dataset_id
  lab <- meta$labeled
  margins <- numeric(0)
  for (i in which(lab)) {
    others <- which(lab & ds != ds[i])
    if (length(others) < K) next
    v <- S[i, others]
    keep <- others[is.finite(v)]
    if (!length(keep)) next
    top <- keep[order(-S[i, keep])][seq_len(min(K, length(keep)))]
    w <- pmax(S[i, top], 0)
    if (sum(w) <= 0) next
    tal <- tapply(w, meta$compartment[top], sum)
    rk <- sort(as.numeric(tal), decreasing = TRUE)
    margins <- c(margins, (rk[1] - if (length(rk) > 1) rk[2] else 0) / sum(w))
  }
  margin_gate <- if (length(margins)) {
    as.numeric(stats::quantile(margins, 0.10, type = 7, names = FALSE))
  } else 0

  saveRDS(list(meta = meta, rank = rank_mat,
               markers = mk,
               ood_gate = atlas$ood_gate,
               margin_gate_rank = margin_gate,
               k = K,
               built_at = atlas$built_at,
               scheme = atlas$scheme,
               n_labeled = atlas$n_labeled,
               n_samples = atlas$n_samples,
               n_datasets = atlas$n_datasets,
               n_layers = atlas$n_layers,
               min_markers = atlas$min_markers,
               skipped = atlas$skipped),
          file.path(out, "reference_atlas.rds"), compress = "xz")
  message(sprintf("  reference_atlas.rds    %d samples x %d markers, labeled=%d, margin_gate=%.4f",
                  nrow(rank_mat), ncol(rank_mat), sum(meta$labeled), margin_gate))
}

## ---- 3. 引擎参照：clock / purity / attribution ------------------------
refs <- json_ref("toolkit_refs.json")
if (!is.null(refs)) {
  saveRDS(refs, file.path(out, "toolkit_refs.rds"), compress = "xz")
  message(sprintf("  toolkit_refs.rds       blocks=%s", paste(names(refs), collapse = "/")))
}

## ---- 3b. 投影的留一数据集验证结果（引擎缓存的同一份数字）---------------
val <- json_ref("reference_projection_validation.json")
if (!is.null(val)) {
  saveRDS(val, file.path(out, "reference_projection_validation.rds"), compress = "xz")
  message(sprintf("  reference_projection_validation.rds  spaces=%s",
                  paste(names(val$spaces), collapse = "/")))
}

## ---- 4. 知识库：通路与膜脂先验（原样 TSV，R 侧 base::read.delim 读）----
## 注意两处不在同一目录：参照 JSON 在 app/backend/data/，知识库 TSV 在
## app/data/kb/（与引擎的 DATA_DIR 一致）。
kb_src <- Sys.getenv("ERYWB_KB", unset = file.path("..", "..", "app", "data", "kb"))
kb_files <- c("pathways.tsv", "pathway_genes.tsv", "pathway_metabolites.tsv",
              "lipid_regulators.tsv", "lipid_regulator_sources.tsv")
for (f in kb_files) {
  p <- file.path(kb_src, f)
  if (file.exists(p)) {
    file.copy(p, file.path(out, f), overwrite = TRUE)
    message(sprintf("  %-26s %d rows", f, length(readLines(p)) - 1L))
  } else {
    warning("missing kb file: ", f, call. = FALSE)
  }
}

## ---- 4b. 特征名标准化表：canonical() 的 display map -------------------
## 不带这份表，R 侧只能做「去版本号 + 符号正则校验」，别名（UniProt、
## Ensembl、旧符号）就归不到标准名，与网站的 to_space() 会产生差异。
dm <- file.path("..", "..", "app", "data", "gene_display_map.tsv")
if (file.exists(dm)) {
  tab <- utils::read.delim(dm, stringsAsFactors = FALSE, quote = "",
                           col.names = c("feature_id", "display_name", "molecule_type", "source"))
  keep <- !is.na(tab$display_name) & nzchar(tab$display_name)
  map <- stats::setNames(toupper(tab$display_name[keep]), toupper(tab$feature_id[keep]))
  map <- map[!duplicated(names(map))]
  saveRDS(map, file.path(out, "gene_display_map.rds"), compress = "xz")
  message(sprintf("  gene_display_map.rds   %d aliases", length(map)))
} else {
  warning("gene_display_map.tsv not found; aliases will not be resolved", call. = FALSE)
}

## ---- 4c. 红系衰老：膜脂参考面板 + 阶段标签 ----------------------------
## 只导出**原始**面板与阶段名，闭合 / 标准化 / 插值留在 R 侧重建：审稿人可以
## 直接核对 45 个键是怎么从 151 行聚出来的，而不是相信一个打包好的曲线。
lip_mat <- file.path("..", "..", "app", "data", "matrices",
                     "PMC11868425_lipidomics", "matrix.tsv")
lip_json <- file.path("..", "..", "app", "data", "disease_focus",
                      "lipid_maturation.json")
if (file.exists(lip_mat) && file.exists(lip_json)) {
  lm_dat <- utils::read.delim(lip_mat, row.names = 1, check.names = FALSE,
                              quote = "", stringsAsFactors = FALSE)
  lm_dat[] <- lapply(lm_dat, function(v) suppressWarnings(as.numeric(v)))
  ## 引擎只用在五个阶段都测到的分子种（complete.cases），未测到的行不参与
  ## 平均值，也不进键表——否则「某阶段缺失」会被当成 0 参与闭合。
  comp <- lm_dat[stats::complete.cases(lm_dat), , drop = FALSE]
  stages <- colnames(comp)

  ## 阶段标签从疾病聚焦页那份 JSON 取，保证两处措辞一致。零依赖解析：
  ## 只认 {"key": "...", "label": {"zh": "...", "en": "..."}} 这一种写法。
  jtxt <- paste(readLines(lip_json, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  pat <- "\"key\"\\s*:\\s*\"([A-Za-z0-9_]+)\"\\s*,\\s*\"label\"\\s*:\\s*\\{\\s*\"zh\"\\s*:\\s*\"([^\"]*)\"\\s*,\\s*\"en\"\\s*:\\s*\"([^\"]*)\""
  hits <- regmatches(jtxt, gregexpr(pat, jtxt, perl = TRUE))[[1]]
  sm <- do.call(rbind, lapply(hits, function(h) {
    g <- regmatches(h, regexec(pat, h, perl = TRUE))[[1]]
    data.frame(key = g[2], zh = g[3], en = g[4], stringsAsFactors = FALSE)
  }))
  if (is.null(sm) || nrow(sm) != length(stages)) {
    ## 解析不出来就退回列名：位置还在，只是少了人话标签
    warning("stage labels not parsed from lipid_maturation.json; using column names",
            call. = FALSE)
    sm <- data.frame(key = stages, zh = stages, en = stages, stringsAsFactors = FALSE)
  }
  sm$position <- seq_len(nrow(sm)) - 1

  panel <- list(raw = as.matrix(comp), stages = stages, stage_meta = sm,
                dataset = "PMC11868425_lipidomics")
  saveRDS(panel, file.path(out, "lipid_panel.rds"), compress = "xz")
  message(sprintf("  lipid_panel.rds        %d complete rows x %d stages",
                  nrow(comp), ncol(comp)))
} else {
  warning("lipid reference panel or stage json not found; ery_age() will not work",
          call. = FALSE)
}

## ---- 5. 溯源清单：源文件指纹 + 引擎版本 ------------------------------
prov <- list(
  exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  engine_version = tryCatch(readLines(file.path("..", "..", "VERSION.txt"), n = 1),
                            error = function(e) NA_character_),
  r_version = R.version.string,
  sources = {
    paths <- c(file.path(src, c("layer_references.json", "reference_atlas.json",
                                "toolkit_refs.json")),
               file.path(kb_src, kb_files),
               file.path("..", "..", "app", "data", "matrices",
                         "PMC11868425_lipidomics", "matrix.tsv"),
               file.path("..", "..", "app", "data", "disease_focus",
                         "lipid_maturation.json"))
    paths <- paths[file.exists(paths)]
    setNames(vapply(paths, sha256, ""), basename(paths))
  }
)
saveRDS(prov, file.path(out, "provenance.rds"), compress = "xz")

## 人读版，便于审稿人直接核对
writeLines(c(
  "# ErythroWorkbench reference bundle",
  sprintf("exported_at: %s", prov$exported_at),
  sprintf("engine_version: %s", prov$engine_version),
  sprintf("r_version: %s", prov$r_version),
  "",
  "## source file fingerprints (md5)",
  sprintf("  %-32s %s", names(prov$sources), unname(prov$sources))
), file.path(out, "PROVENANCE.txt"))

message("\n== inst/extdata ==")
fs <- file.info(list.files(out, full.names = TRUE))
print(data.frame(file = basename(rownames(fs)), KB = round(fs$size / 1024, 1), row.names = NULL))
message(sprintf("total: %.1f KB", sum(fs$size) / 1024))
