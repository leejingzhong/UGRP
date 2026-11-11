"""把 src/*.md 转成 notebooks/*.ipynb。

约定：```r 围栏内是代码格，其余全部是 markdown 格。
这样写课程内容比手写 JSON 快一个数量级，也不容易出错。

用法：
    python tools/md2nb.py            # 转换 src/ 下全部 .md
    python tools/md2nb.py S01        # 只转换名字含 S01 的
"""
import json
import sys
import uuid
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
BASE = HERE.parent
SRC = BASE / "src"
OUT = BASE / "notebooks"

KERNEL = {
    "display_name": "R (4.6.1)",
    "language": "R",
    "name": "ir",
}


def parse(md_text):
    """把 markdown 文本切成 (kind, source) 序列。"""
    cells = []
    buf = []
    lines = md_text.split("\n")
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.lower() in ("```r", "```{r}"):
            if buf:
                cells.append(("md", "\n".join(buf)))
                buf = []
            i += 1
            code = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code.append(lines[i])
                i += 1
            i += 1  # 跳过收尾围栏
            cells.append(("code", "\n".join(code)))
        elif stripped.startswith("```"):
            # 非 R 的围栏（例如输出示例），原样保留为 markdown
            block = [line]
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                block.append(lines[i])
                i += 1
            if i < len(lines):
                block.append(lines[i])
                i += 1
            buf.extend(block)
        else:
            buf.append(line)
            i += 1
    if buf:
        cells.append(("md", "\n".join(buf)))

    # 去掉首尾空行，丢弃空单元格
    cleaned = []
    for kind, src in cells:
        src = src.strip("\n")
        if src.strip():
            cleaned.append((kind, src))
    return cleaned


def to_lines(src):
    """nbformat 要求每行（最后一行除外）自带换行符，否则所有行会被拼成一行。"""
    lines = src.split("\n")
    return [ln + "\n" for ln in lines[:-1]] + [lines[-1]]


def build(cells):
    nb_cells = []
    for kind, src in cells:
        cid = uuid.uuid4().hex[:8]
        if kind == "md":
            nb_cells.append({
                "cell_type": "markdown",
                "id": cid,
                "metadata": {},
                "source": to_lines(src),
            })
        else:
            nb_cells.append({
                "cell_type": "code",
                "id": cid,
                "metadata": {},
                "execution_count": None,
                "outputs": [],
                "source": to_lines(src),
            })
    return {
        "cells": nb_cells,
        "metadata": {
            "kernelspec": KERNEL,
            "language_info": {"name": "R"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main():
    pattern = sys.argv[1] if len(sys.argv) > 1 else ""
    OUT.mkdir(parents=True, exist_ok=True)
    files = sorted(SRC.glob("*.md"))
    if not files:
        print("src/ 下没有 .md 文件")
        return
    for f in files:
        if pattern and pattern not in f.name:
            continue
        nb = build(parse(f.read_text(encoding="utf-8")))
        dst = OUT / (f.stem + ".ipynb")
        dst.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"{f.name:38s} -> {dst.name}  ({len(nb['cells'])} cells)")


if __name__ == "__main__":
    main()
