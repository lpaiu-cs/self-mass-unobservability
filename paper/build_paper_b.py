#!/usr/bin/env python3
"""Build the Paper B LaTeX draft from the markdown source.

Adapted from paper/build_paper_a.py (branch main) with three extensions:
fenced code blocks -> verbatim, wide tables -> proportional p-columns, and
the **Note:** metadata line included in the front-matter block; dash bullets
are accepted alongside asterisk bullets.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_SOURCE = ROOT / "paper-B-dynamic-sep-limit.md"
DEFAULT_OUTPUT = ROOT / "paper-B-dynamic-sep-limit.tex"

SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
}

UNICODE_REPLACEMENTS = {
    "–": "--",
    "—": "---",
    "‘": "`",
    "’": "'",
    "“": "``",
    "”": "''",
    " ": " ",
    "É": "E",
    "é": "e",
    "è": "e",
    "ö": "o",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build the Paper B LaTeX draft.")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def parse_document(text: str) -> tuple[str, dict[str, str], str, list[str]]:
    lines = text.splitlines()
    if not lines or not lines[0].startswith("# "):
        raise ValueError("Paper B source must start with a level-1 title")

    title = lines[0][2:].strip()
    metadata: dict[str, str] = {}
    i = 1

    while i < len(lines) and not lines[i].strip():
        i += 1

    while i < len(lines):
        stripped = lines[i].strip()
        if not stripped:
            i += 1
            continue
        if stripped == "---":
            i += 1
            break
        match = re.match(r"^\*\*([^*]+):\*\*\s*(.*)$", stripped)
        if not match:
            break
        metadata[match.group(1).strip().lower()] = match.group(2).strip()
        i += 1

    while i < len(lines) and not lines[i].strip():
        i += 1

    if i >= len(lines) or lines[i].strip() != "## Abstract":
        raise ValueError("Paper B source must contain a '## Abstract' section")
    i += 1

    while i < len(lines) and not lines[i].strip():
        i += 1

    abstract_lines: list[str] = []
    in_fence = False
    while i < len(lines):
        if lines[i].strip().startswith("```"):
            in_fence = not in_fence
        if lines[i].strip() == "---" and not in_fence:
            i += 1
            break
        abstract_lines.append(lines[i])
        i += 1

    while i < len(lines) and not lines[i].strip():
        i += 1

    return title, metadata, "\n".join(abstract_lines).strip(), lines[i:]


def strip_wrapping_backticks(value: str) -> str:
    value = value.strip()
    if value.startswith("`") and value.endswith("`"):
        return value[1:-1]
    return value


def latex_escape(text: str) -> str:
    for old, new in UNICODE_REPLACEMENTS.items():
        text = text.replace(old, new)
    return "".join(SPECIALS.get(ch, ch) for ch in text)


def convert_code_and_math(text: str) -> str:
    parts = re.split(r"(`[^`]+`|\\\(.+?\\\)|\$[^$]+\$)", text)
    converted: list[str] = []
    for part in parts:
        if not part:
            continue
        if part.startswith("`") and part.endswith("`"):
            converted.append(r"\texttt{" + latex_escape(part[1:-1]) + "}")
        elif part.startswith(r"\(") and part.endswith(r"\)"):
            converted.append("$" + part[2:-2] + "$")
        elif part.startswith("$") and part.endswith("$"):
            converted.append(part)
        else:
            converted.append(latex_escape(part))
    return "".join(converted)


def inline_latex(text: str) -> str:
    pattern = re.compile(r"(\*\*.*?\*\*|\*.*?\*)")
    parts = pattern.split(text)
    converted: list[str] = []
    for part in parts:
        if not part:
            continue
        if part.startswith("**") and part.endswith("**"):
            converted.append(r"\textbf{" + inline_latex(part[2:-2]) + "}")
        elif part.startswith("*") and part.endswith("*"):
            converted.append(r"\emph{" + inline_latex(part[1:-1]) + "}")
        else:
            converted.append(convert_code_and_math(part))
    return "".join(converted)


def heading_latex(line: str) -> str | None:
    match = re.match(r"^(#{2,4})\s+(.*)$", line)
    if not match:
        return None
    level = len(match.group(1))
    title = re.sub(r"^\d+(?:\.\d+)*\.?\s+", "", match.group(2).strip())
    title_tex = inline_latex(title)
    command = {2: "section", 3: "subsection", 4: "subsubsection"}[level]
    return rf"\{command}{{{title_tex}}}"


def table_latex(lines: list[str]) -> str:
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in lines]
    if len(rows) < 2:
        return ""
    header = rows[0]
    body = rows[2:]
    maxlens = [
        max(len(r[k]) if k < len(r) else 0 for r in [header] + body)
        for k in range(len(header))
    ]
    if max(maxlens) > 40:
        total = sum(maxlens)
        widths = [max(0.07, 0.92 * ml / total) for ml in maxlens]
        budget = 0.92 - 0.02 * len(maxlens)   # leave room for column padding
        scale = budget / sum(widths)
        widths = [w * scale for w in widths]
        cols = "".join(r"p{%.2f\textwidth}" % w for w in widths)
    else:
        cols = "l" * len(header)

    def row(cells: list[str]) -> str:
        padded = cells + [""] * (len(header) - len(cells))
        return " & ".join(inline_latex(cell) for cell in padded[: len(header)]) + r" \\"

    out = [r"\begin{center}", rf"\begin{{tabular}}{{{cols}}}", r"\hline", row(header), r"\hline"]
    out.extend(row(cells) for cells in body)
    out.extend([r"\hline", r"\end{tabular}", r"\end{center}"])
    return "\n".join(out)


def list_latex(lines: list[str], ordered: bool) -> str:
    env = "enumerate" if ordered else "itemize"
    items: list[str] = []
    current: list[str] = []
    marker = re.compile(r"^\s*(?:\d+\.|\*|-)\s+(.*)$")
    for line in lines:
        match = marker.match(line)
        if match:
            if current:
                items.append(" ".join(part.strip() for part in current))
            current = [match.group(1)]
        elif line.strip():
            current.append(line.strip())
    if current:
        items.append(" ".join(part.strip() for part in current))

    out = [rf"\begin{{{env}}}"]
    out.extend(r"\item " + inline_latex(item) for item in items)
    out.append(rf"\end{{{env}}}")
    return "\n".join(out)


def paragraph_latex(lines: list[str]) -> str:
    return inline_latex(" ".join(line.strip() for line in lines))


def body_latex(lines: list[str]) -> str:
    out: list[str] = []
    para: list[str] = []
    i = 0

    def flush_para() -> None:
        nonlocal para
        if para:
            out.append(paragraph_latex(para))
            para = []

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        if not stripped or stripped == "---":
            flush_para()
            i += 1
            continue

        heading = heading_latex(line)
        if heading:
            flush_para()
            out.append(heading)
            i += 1
            continue

        if stripped.startswith("```"):
            flush_para()
            code_lines: list[str] = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("```"):
                code_lines.append(lines[i])
                i += 1
            if i >= len(lines):
                raise ValueError("Unclosed fenced code block")
            i += 1
            out.append(
                "\\begin{verbatim}\n" + "\n".join(code_lines) + "\n\\end{verbatim}"
            )
            continue

        if stripped == r"\[":
            flush_para()
            math_lines: list[str] = []
            i += 1
            while i < len(lines) and lines[i].strip() != r"\]":
                math_lines.append(lines[i])
                i += 1
            if i >= len(lines):
                raise ValueError("Unclosed display math block")
            out.append(r"\[" + "\n" + "\n".join(math_lines) + "\n" + r"\]")
            i += 1
            continue

        if stripped.startswith("|"):
            flush_para()
            table_lines: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i])
                i += 1
            out.append(table_latex(table_lines))
            continue

        if re.match(r"^\s*(?:\*|-)\s+", line):
            flush_para()
            list_lines: list[str] = []
            while i < len(lines):
                candidate = lines[i]
                if (
                    re.match(r"^\s*(?:\*|-)\s+", candidate)
                    or (candidate.startswith("  ") and candidate.strip())
                    or not candidate.strip()
                ):
                    list_lines.append(candidate)
                    i += 1
                    continue
                break
            out.append(list_latex(list_lines, ordered=False))
            continue

        if re.match(r"^\s*\d+\.\s+", line):
            flush_para()
            list_lines = []
            while i < len(lines):
                candidate = lines[i]
                if (
                    re.match(r"^\s*\d+\.\s+", candidate)
                    or (candidate.startswith("  ") and candidate.strip())
                    or not candidate.strip()
                ):
                    list_lines.append(candidate)
                    i += 1
                    continue
                break
            out.append(list_latex(list_lines, ordered=True))
            continue

        para.append(line)
        i += 1

    flush_para()
    return "\n\n".join(out).strip()


def build_tex(title: str, metadata: dict[str, str], abstract: str, body: list[str], source: Path) -> str:
    status = metadata.get("status", "")
    repo = metadata.get("repository", "")
    author = strip_wrapping_backticks(metadata.get("author", ""))
    if "replace the GitHub alias" in author:
        author = ""

    note_lines = []
    if status:
        note_lines.append(inline_latex(status))
    if repo:
        note_lines.append("Repository: " + inline_latex(repo))
    extra_note = metadata.get("note", "")
    if extra_note:
        note_lines.append(inline_latex(extra_note))
    note_block = ""
    if note_lines:
        note_block = (
            "\\begin{center}\n"
            "\\small\n"
            + "\\\\\n".join(note_lines)
            + "\n\\end{center}\n\n"
        )

    return f"""% Auto-generated by paper/build_paper_b.py from {source.name}.
% Edit the markdown source and rerun the build script.
\\documentclass[11pt]{{article}}

\\usepackage[T1]{{fontenc}}
\\usepackage[utf8]{{inputenc}}
\\usepackage{{lmodern}}
\\usepackage{{microtype}}
\\usepackage[margin=1in]{{geometry}}
\\usepackage{{amsmath,amssymb,bm}}
\\usepackage{{array}}
\\usepackage{{hyperref}}
\\usepackage{{xurl}}

\\hypersetup{{
  colorlinks=true,
  linkcolor=blue,
  urlcolor=blue
}}

\\title{{{inline_latex(title)}}}
\\author{{{inline_latex(author)}}}
\\date{{2026-07-12}}

\\begin{{document}}
\\maketitle

{note_block}\\begin{{abstract}}
{body_latex(abstract.splitlines())}
\\end{{abstract}}

{body_latex(body)}

\\end{{document}}
"""


def main() -> int:
    args = parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    title, metadata, abstract, body = parse_document(source.read_text(encoding="utf-8"))
    output.write_text(build_tex(title, metadata, abstract, body, source), encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
