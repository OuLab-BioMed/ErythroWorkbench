#!/usr/bin/env python
"""把 R 源码**字符串字面量**里的非 ASCII 字符转成 \\uXXXX 转义。

为什么必须这么做
----------------
R 在 Windows 上默认按 native locale（GBK）解析源码字节。源码是 UTF-8 时，
字符串会被标记成 `Encoding == "unknown"`，内容仍是 UTF-8 字节。后果有两层：

  1. `print()` 在 GBK 控制台显示成 `<U+7EA2>` 转义，中文用户看到的是乱码；
  2. 更要命的是**正则**：`grepl("红细胞", x, perl = TRUE)` 里 pattern 的编码
     未标记，PCRE 会按字节处理，非 ASCII 的匹配静默失败——不报错，只是永远
     匹配不上。身份判定那种依赖中文关键词的地方会**静默给出错误分类**。

`"\\u7ea2\\u7ec6\\u80de"` 则总是生成带 UTF-8 标记的字符串，与 locale 无关。
所以：字符串一律转义，注释保留中文（注释不参与执行）。

幂等：已经是 \\uXXXX 的部分不会被再次处理。

用法
----
    python data-raw/escape_cjk.py R/*.R      # 就地改写
    python data-raw/escape_cjk.py --check R  # 只报告，不写
"""

from __future__ import annotations

import pathlib
import re
import sys

# 这些位置的 "#" 是颜色值（#2563eb）或注释，不能当字符串处理
_HEXCOLOR = re.compile(r"#[0-9A-Fa-f]{3,8}\b")


def escape_literals(line: str) -> tuple[str, int]:
    """转换一行里的双引号字符串。返回 (新行, 替换字符数)。"""
    out: list[str] = []
    i, n = 0, len(line)
    changed = 0
    in_str = False
    in_comment = False
    while i < n:
        ch = line[i]
        if in_comment:
            out.append(line[i:])
            break
        if not in_str:
            if ch == "#":
                in_comment = True
                out.append(ch)
                i += 1
                continue
            if ch == '"':
                in_str = True
            out.append(ch)
            i += 1
            continue
        # 已在字符串内
        if ch == "\\":
            # 原样保留转义序列（含 \\u、\\\\、\\d 等），并跳过被转义的字符
            out.append(line[i:i + 2])
            i += 2
            continue
        if ch == '"':
            in_str = False
            out.append(ch)
            i += 1
            continue
        if ord(ch) > 127:
            out.append("\\u%04x" % ord(ch))
            changed += 1
            i += 1
            continue
        out.append(ch)
        i += 1
    return "".join(out), changed


def main(argv: list[str]) -> int:
    check = "--check" in argv
    args = [a for a in argv if not a.startswith("--")]
    files: list[pathlib.Path] = []
    for a in args:
        p = pathlib.Path(a)
        files.extend(sorted(p.glob("*.R")) if p.is_dir() else [p])
    if not files:
        print("usage: escape_cjk.py [--check] R/*.R", file=sys.stderr)
        return 2

    total = 0
    for f in files:
        raw = f.read_text(encoding="utf-8")
        lines = raw.split("\n")
        new_lines, n_here = [], 0
        for line in lines:
            new, k = escape_literals(line)
            new_lines.append(new)
            n_here += k
        if n_here:
            total += n_here
            print("%-24s %4d chars" % (f.name, n_here))
            if not check:
                f.write_text("\n".join(new_lines), encoding="utf-8", newline="\n")
    print("%s %d characters in %d file(s)"
          % ("would escape" if check else "escaped", total, len(files)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
