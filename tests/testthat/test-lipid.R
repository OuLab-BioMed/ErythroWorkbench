## 脂质名解析与参考模型：逐名对拍 Python 的 parse_lipid() 与 model()。
##
## 解析器错一个字符，某些分子种就永远对不上参考轨迹，而症状只是"命中数偏低"——
## 看不出是解析错了。所以这里对 180 个名字逐个比。

fx <- function(...) testthat::test_path("fixtures", ...)

test_that("lipid name parsing matches the engine name by name", {
  ## colClasses = "character"：默认的 type.convert 会把 "FALSE" 变成 logical、
  ## 把空 base 变成 NA，比到后面全是类型不符的假失败
  exp <- utils::read.csv(fx("lipid_parse_expected.csv"), colClasses = "character",
                         check.names = FALSE)
  expect_equal(nrow(exp), 180L)
  n_parsed <- 0L
  for (i in seq_len(nrow(exp))) {
    nm <- exp$name[i]
    got <- parse_lipid(nm)
    if (!nzchar(exp$cls[i])) {
      expect_equal(nrow(got), 0L, info = sprintf("should not parse: %s", nm))
      next
    }
    n_parsed <- n_parsed + 1L
    expect_equal(nrow(got), 1L, info = sprintf("should parse: %s", nm))
    expect_equal(got$cls[1], exp$cls[i], info = nm)
    expect_equal(got$c[1], as.numeric(exp$c[i]), info = nm)
    expect_equal(got$db[1], as.numeric(exp$db[i]), info = nm)
    expect_equal(as.character(got$ether[1]), exp$ether[i], info = nm)
    expect_equal(got$base[1], exp$base[i], info = nm)
  }
  expect_equal(n_parsed, sum(nzchar(exp$cls)))
})

test_that("parsing handles the four conventions that coexist in public data", {
  ## 四种写法必须给出同一个键，否则"同一批数据换了个导出工具"就对不上了
  ref <- parse_lipid("PC:34:1 (16:0-18:1w9)")
  expect_equal(ref$cls, "PC"); expect_equal(ref$c, 34); expect_equal(ref$db, 1)
  expect_equal(parse_lipid("PC(34:1)")$c, 34)
  expect_equal(parse_lipid("PC(34:1)")$db, 1)

  ## 醚键与碱基必须能区分，否则 A 级会退化成 B 级
  expect_true(parse_lipid("PC(O-34:1)")$ether)
  expect_true(parse_lipid("PC(34:1)e")$ether)
  expect_equal(parse_lipid("SM(d42:1)")$base, "d")
  expect_equal(parse_lipid("SM(t42:3)")$base, "t")

  ## 链式写法要和总式写法算出同一个总数
  expect_equal(parse_lipid("TG(18:2_17:1_18:2)")$c, 53)
  expect_equal(parse_lipid("TG(18:2_17:1_18:2)")$db, 5)
  expect_equal(parse_lipid("DG(16:0-18:1)")$c, 34)

  ## 全称别名
  expect_equal(parse_lipid("SPHINGOMYELIN 34:1")$cls, "SM")
  expect_equal(parse_lipid("PHOSPHATIDYLCHOLINE(36:2)")$cls, "PC")
})

test_that("the reference trajectory has 45 anchors across five stages", {
  m <- ery_lipid_model()
  expect_equal(nrow(m$keys), 45L)
  expect_equal(m$stages, c("RY", "RM", "EY", "EM", "EO"))
  expect_equal(sum(m$ok), 45L)
  expect_equal(dim(m$raw), c(45L, 5L))
  expect_equal(dim(m$curve), c(45L, 81L))
  expect_equal(m$grid[1], 0); expect_equal(m$grid[81], 4)

  ## 阶段标签来自疾病聚焦页那份 JSON，两处措辞必须一致
  expect_equal(m$stage_meta$key, m$stages)
  expect_true(all(nzchar(m$stage_meta$zh)))

  ## 组成型：每个阶段闭合到 100%，所以每列之和都是 100
  expect_true(all(abs(colSums(m$pct) - 100) < 1e-8))

  ## 45 个键里 SM 的 d 系列必须保留 base，否则 A 级会全部退化成 B 级
  expect_true(any(m$keys$base == "d"))
})

test_that("leave-one-stage-out error is the only accuracy figure quoted", {
  v <- ery_lipid_validation()
  expect_equal(v$n_stages, 5L)
  expect_equal(v$loso_errors, c(0.0, 0.0, 0.1, 0.1, -1.0))
  expect_equal(v$loso_mae, 0.24)
  expect_equal(v$loso_max, 1.0)
  expect_equal(v$loso_mae_pct, 6.0)
})

test_that("identity gate blocks leukocyte and plasma samples", {
  ## 脂质层没有 RNA，闸门走注释/样本名身份；白细胞、血浆、其他组织必须被拦
  m <- matrix(c(1, 2, 3, 4), 2,
              dimnames = list(c("PC:34:1 (16:0-18:1)", "PE:36:2 (18:1-18:1)"),
                              c("RBC_1", "PBMC_1")))
  g <- ErythroWorkbench:::.age_gate(m)
  expect_equal(g$level, "identity")
  expect_equal(g$passed, "RBC_1")
  expect_equal(g$blocked, "PBMC_1")

  ## 一个都认不出身份时不能假装核验过，要明确 unverifiable
  m2 <- matrix(1, nrow = 1, ncol = 1, dimnames = list("PC:34:1", "S1"))
  g2 <- ErythroWorkbench:::.age_gate(m2)
  expect_equal(g2$level, "unverifiable")
  expect_true(grepl("unverifiable|UNVERIFIED", paste(g2$note, g2$note_en)))
})
