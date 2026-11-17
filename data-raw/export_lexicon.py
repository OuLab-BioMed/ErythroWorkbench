#!/usr/bin/env python
"""把引擎的代谢物词典导出成包内资源（inst/extdata/metabolite_lexicon.txt）。

为什么这一步必须在 Python 侧做
------------------------------
词典来自 SQLite 的 pathway_metabolites 与 metabolite_annotations 两张表。
R 包的全部承诺是**零第三方依赖**，而 base R 读不了 SQLite（需要 RSQLite）。
所以这份名单只能在这里生成一次、作为静态资源随包冻结，并在 provenance 里
记下指纹。它是 183 个归一化名字，4 KB，不存在"过期"风险——库里新增代谢物时
重跑这个脚本即可。

归一化口径必须与 toolkit_core.norm_metabolite() 一致（小写、去空格、去连字符
与下划线），否则 R 侧查不到。
"""

from __future__ import annotations

import pathlib
import sqlite3
import sys

REPO = pathlib.Path(__file__).resolve().parents[3]
DB = REPO / "app" / "data" / "eiod.db"
OUT = pathlib.Path(__file__).resolve().parents[1] / "inst" / "extdata"


def main() -> int:
    if not DB.exists():
        print("database not found: %s" % DB, file=sys.stderr)
        return 1
    names: set[str] = set()
    con = sqlite3.connect(DB)
    try:
        for (n,) in con.execute(
                "select distinct metabolite_name from pathway_metabolites"):
            k = norm(n)
            if k:
                names.add(k)
        for a, b in con.execute(
                "select original_feature_name, standardized_name from metabolite_annotations"):
            for n in (a, b):
                k = norm(n)
                if k:
                    names.add(k)
    finally:
        con.close()

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "metabolite_lexicon.txt"
    path.write_text("\n".join(sorted(names)) + "\n", encoding="utf-8",
                    newline="\n")
    print("wrote %d names to %s (%.1f KB)"
          % (len(names), path, path.stat().st_size / 1024))
    return 0


def norm(name) -> str | None:
    """与 toolkit_core.norm_metabolite() 同一口径：小写、去空格、去 - 与 _。"""
    s = str(name or "").strip().lower()
    if not s:
        return None
    s = "".join(s.split())
    s = s.replace("-", "").replace("_", "")
    return s or None


if __name__ == "__main__":
    sys.exit(main())
