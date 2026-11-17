## ---------------------------------------------------------------------------
## 样本身份判定
##
## 移植自 app/backend/engines/sample_identity.py。用途有二：
##   1. 条件谱匹配时判断「这个样本到底是不是红细胞」——对白细胞制品做疾病归因
##      是错位提问，必须先把这句话说出来；
##   2. 红系衰老的身份闸门——年龄轴只在红细胞膜脂轨迹上定义，PBMC / 血浆投上来
##      的「位置」没有生物学解释。
##
## 规则本身不复杂，难的是**判据必须与引擎逐条一致**：正则里少一个词，某些数据集
## 就会被分到别的区室，闸门放行/拦截名单随之改变，而这种差异不会报错。
## 因此下面每条正则都保留 Python 版的原样顺序与权重，中文关键词用 \u 转义写
## （直接写中文会让 Encoding 变成 unknown，PCRE 按字节处理，匹配静默失败）。
## ---------------------------------------------------------------------------

## 细胞系名单绝对优先：培养条件已经把谱系改造过，不能与「原代红细胞」混为一谈
.si_rx_cell_line_nonery <- paste0(
  "hela|hek.?293|\\bjurkat\\b|mcf.?7|\\bthp.?1\\b|\\bu937\\b|\\bdt40\\b|\\bcho\\b|",
  "nih.?3t3|\\ba549\\b|hct.?116|\\bhela\\b|lymphoblastoid|ebv.?transformed|",
  "\\bimr.?90\\b|\\bshsy5y\\b|\\bpc.?3\\b|\\blncap\\b|\\bhuh.?7\\b|\\bhepg2\\b|",
  "ebv.?immortal")
.si_rx_cell_line_ery <- paste0(
  "hudep|k562|\\bhel\\b|ut.?7|erythroleukemi|mel.?cell|\\bhela\\b(?!.*non)")

.si_rx_ery_strong <- paste0(
  "erythrocyt|reticulocyt|\\brbc\\d*\\b|red blood cell|red cell|normoblast|",
  "proerythroblast|erythroblast|megaloblast|\\bbfu.?e\\b|\\bcfu.?e\\b|",
  "glycophorin|\\bgypa\\b|band.?3\\b|slc4a1|ter119|\\bcd71\\b|",
  "\u7ea2\u7ec6\u80de|\u7f51\u7ec7\u7ea2|\u7ea2\u7cfb\u7956|",
  "\u539f\u7ea2\u7ec6\u80de|\u5e7c\u7ea2\u7ec6\u80de")
.si_rx_ery_weak <- paste0(
  "erythroid|h[ae]moglobin|\\bglobin\\b|\\bhba\\b|\\bhbb\\b|",
  "\u7ea2\u7cfb|\u7f51\u7ec7|\u8840\u7ea2\u7d20|",
  "\u8840\u7ea2\u86cb\u767d|\u8f6c\u94c1\u86cb\u767d\u53d7\u4f53")
.si_rx_leu_strong <- paste0(
  "pbmc|peripheral blood mononuclear|mononuclear|leukocyt|\\bwbc\\b|",
  "lymphocyt|\\bt.?cell\\b|\\bb.?cell\\b|\\bnk.?cell\\b|monocyt|neutrophil|",
  "eosinophil|basophil|platelet|megakaryocyt|macrophage|dendritic|",
  "mast.?cell|myeloid|granulocyt|",
  "\u5355\u4e2a\u6838|\u767d\u7ec6\u80de|\u6dcb\u5df4\u7ec6\u80de|\u5355\u6838\u7ec6\u80de|",
  "\u4e2d\u6027\u7c92|\u55dc\u9178|\u55dc\u78b1|\u8840\u5c0f\u677f|",
  "\u5de8\u566c|\u6811\u7a81")
.si_rx_plasma <- paste0(
  "\\bplasma\\b|\\bserum\\b|extracellular vesicle|exosom|small ev|\\bev\\b|",
  "microparticle|\u8840\u6d46|\u8840\u6e05|\u5916\u6ccc\u4f53|\u56ca\u6ce1")
.si_rx_other_tissue <- paste0(
  "kidney|renal|hepat|liver|brain|cortex|hippocamp|\\blung\\b|airway|",
  "skin|dermis|epiderm|muscle|cardiac|heart|adipos|intestin|colon|",
  "\\bgut\\b|synovi|cartilag|osteoblast|retina|cornea|pancrea|endotheli|",
  "fibroblast|epitheli|tumou?r|carcinoma|xenograft|placenta|",
  "\u809d|\u80be|\u8111|\u80ba|\u76ae\u80a4|\u808c\u8089|\u5fc3\u810f|",
  "\u8102\u80aa|\u80a0|\u7ed3\u80a0|\u6ed1\u819c|\u8f6f\u9aa8|\u9aa8|",
  "\u89c6\u7f51\u819c|\u89d2\u819c|\u80f0\u817a|\u5185\u76ae|",
  "\u6210\u7ea4\u7ef4|\u4e0a\u76ae|\u80bf\u7624|\u764c")
.si_rx_whole_blood <- paste0(
  "whole blood|\u5168\u8840|peripheral blood(?! mononuclear)$|",
  "\\bblood\\d*\\b|venous blood|cord blood|\u8110\u8840|umbilical")

## 区室与放行规则：unknown 放行但标低置信（不知道不等于排除）
.si_scope <- c(ery_mature = TRUE, ery_progenitor = TRUE, ery_cell_line = TRUE,
               mixed_blood = TRUE, leukocyte = FALSE, plasma = FALSE,
               other_tissue = FALSE, unknown = TRUE)
.si_block <- c("leukocyte", "plasma", "other_tissue", "ery_cell_line")

## 把样本名与该行所有注释拼成一段文本（去重、限长 600 字符）
.si_text_of <- function(name, row = NULL) {
  parts <- as.character(name)
  if (!is.null(row) && length(row)) {
    vals <- vapply(row, function(v) {
      s <- trimws(as.character(v)[1])
      if (is.na(s) || !nzchar(s) || tolower(s) %in% c("nan", "none", "unknown", "na")) "" else s
    }, "")
    parts <- c(parts, vals[nzchar(vals)])
  }
  parts <- unique(parts)
  substr(paste(parts, collapse = " | "), 1, 600)
}

.si_hit <- function(pat, txt) {
  if (!nzchar(txt)) return(character(0))
  m <- regexpr(pat, txt, perl = TRUE, ignore.case = TRUE)
  if (m[1] == -1L) return(character(0))
  regmatches(txt, m)
}

.si_out <- function(name, comp, hit, txt, strength) {
  list(sample = as.character(name), compartment = comp,
       in_scope = isTRUE(.si_scope[[comp]]),
       evidence = utils::head(hit, 6),
       confidence = if (strength >= 3) "high" else if (strength == 2) "medium" else "low",
       text = substr(txt, 1, 200))
}

## 判定单个样本。row 可以是 list / 单行 data.frame / NULL。
.si_classify <- function(name, row = NULL) {
  txt <- .si_text_of(name, row)
  hit <- character(0)
  score <- 0

  ## 1) 细胞系名单绝对优先
  m_non_line <- .si_hit(.si_rx_cell_line_nonery, txt)
  m_ery_line <- .si_hit(.si_rx_cell_line_ery, txt)
  if (length(m_non_line)) {
    comp <- "other_tissue"
    hit <- c(hit, paste0("cell_line:", m_non_line[1]))
    ## EBV 淋巴母细胞归 leukocyte 更准确（血液来源、淋巴谱系）
    if (grepl("lymphoblastoid|ebv", txt, ignore.case = TRUE, perl = TRUE)) comp <- "leukocyte"
    return(.si_out(name, comp, hit, txt, 3))
  }
  if (length(m_ery_line)) {
    hit <- c(hit, paste0("ery_cell_line:", m_ery_line[1]))
    return(.si_out(name, "ery_cell_line", hit, txt, 3))
  }

  ## 2) 加权证据：红系 +3/+2，白细胞/血浆/其他组织各 -3
  flags <- list()
  for (spec in list(list(.si_rx_ery_strong, 3, "ery_strong"),
                    list(.si_rx_ery_weak, 2, "ery_weak"),
                    list(.si_rx_leu_strong, -3, "leukocyte"),
                    list(.si_rx_plasma, -3, "plasma"),
                    list(.si_rx_other_tissue, -3, "other_tissue"))) {
    ms <- .si_hit(spec[[1]], txt)
    flags[[spec[[3]]]] <- length(ms) > 0
    if (length(ms)) {
      score <- score + spec[[2]]
      hit <- c(hit, paste0(spec[[3]], ":", ms[1]))
    }
  }
  whole <- length(.si_hit(.si_rx_whole_blood, txt)) > 0
  if (whole) hit <- c(hit, "whole_blood")

  prog_rx <- paste0("progenitor|precursor|\\bbfu|\\bcfu|cd34|erythroblast|",
                    "\u7ea2\u7cfb\u7956|\u7956\u7ec6\u80de|\u524d\u4f53")

  if (score >= 3) {
    comp <- if (grepl(prog_rx, txt, ignore.case = TRUE, perl = TRUE)) {
      "ery_progenitor"
    } else "ery_mature"
  } else if (score >= 2 && grepl("\u7ea2\u7cfb|erythroid", txt, ignore.case = TRUE, perl = TRUE)) {
    comp <- if (grepl("progenitor|precursor|bfu|cfu|cd34|\u7ea2\u7cfb\u7956|\u7956\u7ec6\u80de",
                      txt, ignore.case = TRUE, perl = TRUE)) "ery_progenitor" else "ery_mature"
  } else if (score <= -3) {
    ## 血浆优先于白细胞：血浆样本里也常写「去除白细胞」
    comp <- if (isTRUE(flags$plasma)) "plasma"
            else if (isTRUE(flags$leukocyte)) "leukocyte" else "other_tissue"
  } else if (score == 0 && (isTRUE(flags$ery_strong) || isTRUE(flags$ery_weak)) &&
             (isTRUE(flags$plasma) || isTRUE(flags$leukocyte))) {
    ## 红系与血浆/白细胞证据等强（「红细胞+血浆」这类配对研究）：是同一份血的
    ## 不同组分，不是污染——单列一档，既不丢也不当成纯红细胞。
    comp <- "mixed_blood"
  } else if (whole && score < 3) {
    ## 全血/脐血：含红系但未分选。只有在**没有**明确红系强证据时才降级，
    ## 否则「脐血 + cell_type=成熟红细胞」会被错降成混合血，纯度整体偏低。
    comp <- "mixed_blood"
  } else {
    comp <- "unknown"
  }
  .si_out(name, comp, hit, txt, abs(score))
}

## ---- 注释对齐：键 → 文本 → 位置 + 分块检验 --------------------------------

.si_key_align <- function(cols, ann) {
  keys <- c("sample_id", "matrix_column", "sample", "sampleid", "gsm", "id")
  n <- length(cols)
  for (key in keys) {
    idx <- which(tolower(trimws(names(ann))) == key)
    if (!length(idx)) next
    col <- names(ann)[idx[1]]
    vals <- trimws(as.character(ann[[col]]))
    m <- !duplicated(vals)
    lookup <- stats::setNames(seq_len(nrow(ann))[m], vals[m])
    hit <- lookup[cols]
    if (sum(!is.na(hit)) >= max(1, as.integer(0.5 * n))) {
      return(list(rows = hit, how = key))
    }
  }
  for (col in names(ann)) {
    vals <- trimws(as.character(ann[[col]]))
    m <- !duplicated(vals)
    lookup <- stats::setNames(seq_len(nrow(ann))[m], vals[m])
    hit <- lookup[cols]
    if (sum(!is.na(hit)) >= max(1, as.integer(0.8 * n))) {
      return(list(rows = hit, how = col))
    }
  }
  list(rows = NULL, how = NULL)
}

## 注释的文字列（notes / title / summary）里常常写着矩阵列名。
## GSE63703 的 notes 是 "GEO sample GSM1830782 (Long RNA-seq: RBC_1)"——
## 列名就在里面，比位置对齐可靠得多。
.si_text_align <- function(cols, ann) {
  want <- unique(cols)
  if (!length(want)) return(list(rows = NULL, how = NULL))
  best <- NULL
  best_hit <- 0L
  for (col in names(ann)) {
    key <- tolower(trimws(col))
    if (key %in% c("sample_id", "matrix_column", "dataset_id", "gsm", "id")) next
    m <- integer(0)
    names(m) <- character(0)
    hit <- 0L
    for (i in seq_len(nrow(ann))) {
      blob <- as.character(ann[[col]][i])[1]
      if (is.na(blob) || !nzchar(blob) || tolower(blob) %in% c("nan", "none")) next
      toks <- regmatches(blob, gregexpr("[A-Za-z0-9_.-]+", blob, perl = TRUE))[[1]]
      for (tk in toks) {
        if (tk %in% want && !(tk %in% names(m))) {
          m[tk] <- i
          hit <- hit + 1L
          break
        }
      }
    }
    if (hit > best_hit) {
      best <- m
      best_hit <- hit
    }
  }
  if (!is.null(best) && best_hit >= max(1, as.integer(0.8 * length(want)))) {
    return(list(rows = unname(best[cols]), how = paste0("text:", best_hit)))
  }
  list(rows = NULL, how = NULL)
}

.si_runs <- function(seq) {
  if (length(seq) < 2) return(length(seq))
  1L + sum(seq[-1] != seq[-length(seq)])
}

## 行数 == 列数时按位置对齐，并做**分块一致性检验**：把列名剥掉尾部编号得到
## 「命名块」（NF54_plusCHO_rep1..3 → NF54_plusCHO），块数与注释任一条件列的
## 取值块数相同，才认为顺序一致。顺序错了身份判定不受影响（cell_type 是数据集
## 级常量），但分组会全错，所以这里只把「验证过」的标出来。
.si_position_align <- function(cols, ann) {
  if (nrow(ann) != length(cols)) return(list(rows = NULL, how = NULL))
  heads <- sub("[_-]?(rep|r|s)?[_-]?[0-9]+.*$", "", cols, ignore.case = TRUE)
  hr <- .si_runs(heads)
  for (col in names(ann)) {
    if (tolower(trimws(col)) %in% c("sample_id", "dataset_id", "gsm", "id",
                                    "matrix_column", "platform", "batch")) next
    vals <- trimws(as.character(ann[[col]]))
    nu <- length(unique(vals))
    if (nu < 2 || nu > 12) next
    if (.si_runs(vals) == hr) {
      return(list(rows = seq_len(nrow(ann)),
                  how = paste0("position:verified(", col, ")")))
    }
  }
  list(rows = seq_len(nrow(ann)), how = "position:unverified")
}

.si_row_for <- function(cols, ann) {
  if (is.null(ann) || !nrow(ann)) return(list(rows = NULL, how = "no_annotation"))
  cols <- trimws(as.character(cols))
  for (f in list(.si_key_align, .si_text_align, .si_position_align)) {
    r <- f(cols, ann)
    if (!is.null(r$rows)) {
      how <- r$how
      if (f == .si_key_align) how <- paste0("key:", how)
      return(list(rows = r$rows, how = how))
    }
  }
  list(rows = NULL, how = "failed")
}

## 注释缺行时，用同前缀样本的身份补齐（D8_3 缺注释 → 继承 D8_1）。
## 只在同数据集内、同前缀、且已有判定不是 unknown 时继承。
.si_prefix <- function(name) {
  s <- as.character(name)
  m <- regmatches(s, regexec("^(.*?)[_\\-.]?(rep|r|s)?[_\\-.]?[0-9]+[a-z]?[_\\-.]?$",
                             s, ignore.case = TRUE, perl = TRUE))[[1]]
  out <- if (length(m) >= 2 && nzchar(m[2])) m[2] else s
  tolower(out)
}

.si_inherit <- function(samples) {
  known <- list()
  for (s in samples) {
    if (s$compartment != "unknown") {
      p <- .si_prefix(s$sample)
      if (is.null(known[[p]])) known[[p]] <- s$compartment
    }
  }
  n <- 0L
  for (i in seq_along(samples)) {
    s <- samples[[i]]
    if (s$compartment != "unknown") next
    p <- .si_prefix(s$sample)
    if (!is.null(known[[p]])) {
      s$compartment <- known[[p]]
      s$in_scope <- isTRUE(.si_scope[[known[[p]]]])
      s$confidence <- "inherited"
      s$evidence <- c(s$evidence, paste0("prefix:", p))
      samples[[i]] <- s
      n <- n + 1L
    }
  }
  attr(samples, "n_inherited") <- n
  samples
}

#' Judge the sample identity of every column
#'
#' Reads the annotation (if any) and the sample names, and assigns each sample to
#' a compartment: mature erythrocyte, erythroid progenitor, erythroid cell line,
#' mixed blood, leukocyte, plasma, other tissue, or unknown. Only the first four
#' count as erythroid for downstream modules.
#'
#' The rules are ported one-to-one from the web server, so a sample excluded on
#' the site is excluded here too. Identity is decided by weighted keyword
#' evidence (erythroid +3 / weak erythroid +2 / leukocyte, plasma, other tissue
#' -3 each), with cell-line names taking absolute precedence.
#'
#' @param x Numeric matrix with features in rows and samples in columns. Only the
#'   column names are used.
#' @param annotation Optional data.frame of sample annotation. It is aligned to
#'   the matrix columns by key, then by text, then by position with a block
#'   consistency check; each alignment mode is reported in `attr(, "aligned_by")`.
#' @return A data.frame with one row per sample: `sample`, `compartment`,
#'   `in_scope`, `confidence`, `evidence`.
#' @examples
#' m <- matrix(1:6, 3, dimnames = list(c("A", "B", "C"), c("RBC_1", "PBMC_1")))
#' ery_identify(m)
#' @export
ery_identify <- function(x, annotation = NULL) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop("`x` must be a matrix or data.frame.", call. = FALSE)
  cols <- colnames(x)
  if (is.null(cols)) stop("`x` needs column names (sample identifiers).", call. = FALSE)
  if (!is.null(annotation) && !is.data.frame(annotation)) {
    annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  }
  al <- .si_row_for(cols, annotation)
  samples <- vector("list", length(cols))
  for (i in seq_along(cols)) {
    row <- NULL
    if (!is.null(al$rows)) {
      j <- al$rows[i]
      if (!is.na(j)) row <- as.list(annotation[j, , drop = FALSE])
    }
    samples[[i]] <- .si_classify(cols[i], row)
  }
  ## 继承只做一次：注释缺行时用同前缀样本补齐（D8_3 缺行 → 继承 D8_1）
  samples <- .si_inherit(samples)
  n_inh <- attr(samples, "n_inherited")
  out <- data.frame(
    sample = vapply(samples, `[[`, "", "sample"),
    compartment = vapply(samples, `[[`, "", "compartment"),
    in_scope = vapply(samples, `[[`, TRUE, "in_scope"),
    confidence = vapply(samples, `[[`, "", "confidence"),
    evidence = vapply(samples, function(s) paste(s$evidence, collapse = ";"), ""),
    stringsAsFactors = FALSE
  )
  attr(out, "aligned_by") <- al$how
  attr(out, "n_inherited") <- as.integer(n_inh %||% 0L)
  out
}
