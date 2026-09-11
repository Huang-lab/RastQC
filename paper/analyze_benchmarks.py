#!/usr/bin/env python3
"""
Generate SVG figures and markdown tables from RastQC benchmark results.

Reads:   benchmark/results/benchmark_results.csv (from benchmark/run_benchmark.sh)
Outputs: paper/figures/fig1_short_read_speed.svg
         paper/figures/fig2_long_read_speed.svg
         paper/figures/fig3_memory_comparison.svg
         markdown tables on stdout, for pasting into the paper

The tool ids below must match the ones benchmark/run_benchmark.sh writes to the
CSV. They are the contract between the two scripts: a series whose id is not
produced is simply absent from the figures, so a rename on either side loses
data silently rather than erroring.
"""

import csv
import os
import sys

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(PAPER_DIR)
CSV_PATH = os.environ.get(
    "BENCHMARK_CSV",
    os.path.join(PROJECT_DIR, "benchmark", "results", "benchmark_results.csv"),
)
FIG_DIR = os.environ.get("FIGURE_DIR", os.path.join(PAPER_DIR, "figures"))

# (csv tool id, display name, colour)
FASTQC = ("fastqc", "FastQC", "#e74c3c")
FALCO = ("falco", "Falco", "#f39c12")
RASTQC_T1 = ("rastqc_t1", "RastQC -t 1", "#5dade2")
RASTQC_TN = ("rastqc_tN", "RastQC -t N", "#2980b9")
RASTQC_LR = ("rastqc_lr", "RastQC --long-read", "#27ae60")


def load_results():
    with open(CSV_PATH) as f:
        return list(csv.DictReader(f))


def index(rows):
    """{filename: {tool_id: row}}, plus per-file metadata, in CSV order."""
    by_file, meta, order = {}, {}, []
    for r in rows:
        f = r["file"]
        if f not in by_file:
            by_file[f] = {}
            order.append(f)
            meta[f] = {
                "type": r["type"],
                "size_mb": int(r["size_mb"] or 0),
                "reads": int(r["reads"] or 0),
                "mean_len": int(r["mean_read_len"] or 0),
            }
        by_file[f][r["tool"]] = r
    return by_file, meta, order


def short_label(fname):
    return fname.replace(".fastq.gz", "").replace(".fastq", "")


# ─── One grouped-bar renderer, used by all three figures ────────────────────

def grouped_bar_svg(path, title, y_label, groups, series, unit, decimals=1):
    """groups: [{"label", "sublabel", "values": {series_id: float}}]
       series: [(series_id, display_name, colour)] — drawn left to right."""
    series = [s for s in series if any(s[0] in g["values"] for g in groups)]
    if not groups or not series:
        print(f"  (skipped {os.path.basename(path)} — no matching data)")
        return

    vals = [v for g in groups for v in g["values"].values()]
    y_max = max(vals) * 1.18 or 1

    # Round the axis up to a round number so gridlines land on readable values.
    step = 10 ** (len(str(int(y_max))) - 1)
    y_max = (int(y_max / step) + 1) * step
    ticks = [y_max * i / 5 for i in range(6)]

    n_groups, n_series = len(groups), len(series)
    w = max(640, 130 * n_groups + 190)
    h = 420
    ml, mr, mt, mb = 78, 160, 46, 74
    pw, ph = w - ml - mr, h - mt - mb
    gw = pw / n_groups
    bw = min(26, (gw * 0.8) / n_series)

    s = []
    s.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
             f'font-family="Helvetica, Arial, sans-serif">')
    s.append(f'<rect width="{w}" height="{h}" fill="white"/>')
    s.append(f'<text x="{w/2:.0f}" y="24" text-anchor="middle" font-size="14" '
             f'font-weight="bold">{title}</text>')

    for t in ticks:
        y = mt + ph - (t / y_max * ph)
        s.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+pw:.0f}" y2="{y:.1f}" '
                 f'stroke="#e3e3e3" stroke-width="0.5"/>')
        lab = f"{t:.0f}" if t >= 10 or t == 0 else f"{t:.1f}"
        s.append(f'<text x="{ml-8}" y="{y+4:.1f}" text-anchor="end" font-size="10.5" '
                 f'fill="#555">{lab}</text>')

    for gi, g in enumerate(groups):
        block = bw * n_series
        x0 = ml + gw * gi + (gw - block) / 2
        yb = mt + ph
        for si, (sid, _, colour) in enumerate(series):
            if sid not in g["values"]:
                continue
            v = g["values"][sid]
            bh = (v / y_max) * ph
            x = x0 + si * bw
            s.append(f'<rect x="{x:.1f}" y="{yb-bh:.1f}" width="{bw-2:.1f}" '
                     f'height="{bh:.1f}" fill="{colour}" rx="1.5"/>')
            s.append(f'<text x="{x+(bw-2)/2:.1f}" y="{yb-bh-3:.1f}" '
                     f'text-anchor="middle" font-size="7.5" fill="{colour}">'
                     f'{v:.{decimals}f}</text>')
        xc = ml + gw * gi + gw / 2
        s.append(f'<text x="{xc:.0f}" y="{yb+16}" text-anchor="middle" '
                 f'font-size="10" fill="#333">{g["label"]}</text>')
        if g.get("sublabel"):
            s.append(f'<text x="{xc:.0f}" y="{yb+29}" text-anchor="middle" '
                     f'font-size="8.5" fill="#888">{g["sublabel"]}</text>')

    s.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+ph:.0f}" stroke="#333" stroke-width="1.4"/>')
    s.append(f'<line x1="{ml}" y1="{mt+ph:.0f}" x2="{ml+pw:.0f}" y2="{mt+ph:.0f}" stroke="#333" stroke-width="1.4"/>')
    s.append(f'<text x="17" y="{mt+ph/2:.0f}" text-anchor="middle" font-size="11.5" '
             f'fill="#333" transform="rotate(-90 17 {mt+ph/2:.0f})">{y_label}</text>')

    lx, ly = w - mr + 14, mt + 14
    for si, (_, name, colour) in enumerate(series):
        y = ly + si * 20
        s.append(f'<rect x="{lx}" y="{y}" width="12" height="12" fill="{colour}" rx="2"/>')
        s.append(f'<text x="{lx+18}" y="{y+11}" font-size="10.5" fill="#333">{name}</text>')
    s.append(f'<text x="{lx}" y="{ly + n_series*20 + 14}" font-size="9" fill="#888">({unit})</text>')

    s.append("</svg>")
    with open(path, "w") as f:
        f.write("\n".join(s))
    print(f"  saved {os.path.relpath(path, PROJECT_DIR)}")


def build_groups(by_file, meta, order, ftype, value_key, cast=float, skip_agg=True):
    groups = []
    for fname in order:
        m = meta[fname]
        if m["type"] != ftype:
            continue
        if skip_agg and fname.startswith("ALL"):
            continue
        values = {tool: cast(r[value_key]) for tool, r in by_file[fname].items()}
        groups.append({
            "label": short_label(fname),
            "sublabel": f"{m['size_mb']} MB",
            "values": values,
        })
    return sorted(groups, key=lambda g: int(g["sublabel"].split()[0]))


# ─── Markdown tables ─────────────────────────────────────────────────────────

def md_tables(by_file, meta, order):
    def cell(row, key, fmt="{:.2f}"):
        return fmt.format(float(row[key])) if row else "—"

    for ftype, heading in (("short", "Short-read"), ("long", "Long-read")):
        files = [f for f in order if meta[f]["type"] == ftype]
        if not files:
            continue
        print(f"\n### {heading} performance\n")
        cols = ["File", "Size", "Reads", "Mean len", "FastQC", "Falco",
                "RastQC -t 1", "RastQC -t N", "vs FastQC", "vs Falco"]
        if ftype == "long":
            cols.insert(8, "RastQC --long-read")
        print("| " + " | ".join(cols) + " |")
        print("|" + "|".join(["---"] * len(cols)) + "|")

        for fname in files:
            t, m = by_file[fname], meta[fname]
            fq, fa = t.get("fastqc"), t.get("falco")
            r1, rn, rl = t.get("rastqc_t1"), t.get("rastqc_tN"), t.get("rastqc_lr")
            best = float(rn["real_sec"]) if rn else None
            vs_fq = f"{float(fq['real_sec'])/best:.1f}x" if fq and best else "—"
            vs_fa = f"{float(fa['real_sec'])/best:.1f}x" if fa and best else "—"
            row = [
                short_label(fname),
                f"{m['size_mb']} MB",
                f"{m['reads']:,}" if m["reads"] else "—",
                f"{m['mean_len']} bp" if m["mean_len"] else "—",
                cell(fq, "real_sec") + " s" if fq else "—",
                cell(fa, "real_sec") + " s" if fa else "—",
                cell(r1, "real_sec") + " s" if r1 else "—",
                "**" + cell(rn, "real_sec") + " s**" if rn else "—",
                vs_fq, vs_fa,
            ]
            if ftype == "long":
                row.insert(8, cell(rl, "real_sec") + " s" if rl else "—")
            print("| " + " | ".join(row) + " |")

        print(f"\n**Peak resident memory ({heading.lower()})**\n")
        print("| File | FastQC | Falco | RastQC -t 1 | RastQC -t N |")
        print("|---|---|---|---|---|")
        for fname in files:
            t = by_file[fname]
            def mem(tool):
                return f"{t[tool]['max_rss_mb']} MB" if tool in t else "—"
            print(f"| {short_label(fname)} | {mem('fastqc')} | {mem('falco')} "
                  f"| {mem('rastqc_t1')} | {mem('rastqc_tN')} |")


# ─── Main ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: Benchmark results not found at {CSV_PATH}")
        print("Run benchmark/run_benchmark.sh first.")
        sys.exit(1)

    os.makedirs(FIG_DIR, exist_ok=True)
    rows = load_results()
    by_file, meta, order = index(rows)

    print("=" * 68)
    print("RastQC benchmark figures")
    print(f"Source: {os.path.relpath(CSV_PATH, PROJECT_DIR)}")
    print("=" * 68)

    grouped_bar_svg(
        os.path.join(FIG_DIR, "fig1_short_read_speed.svg"),
        "Short-read performance: RastQC vs FastQC and Falco",
        "Wall-clock time (seconds)",
        build_groups(by_file, meta, order, "short", "real_sec"),
        [FASTQC, FALCO, RASTQC_T1, RASTQC_TN],
        "lower is better",
    )
    grouped_bar_svg(
        os.path.join(FIG_DIR, "fig2_long_read_speed.svg"),
        "Long-read performance: RastQC vs FastQC and Falco",
        "Wall-clock time (seconds)",
        build_groups(by_file, meta, order, "long", "real_sec"),
        [FASTQC, FALCO, RASTQC_TN, RASTQC_LR],
        "lower is better",
    )

    mem_groups = (build_groups(by_file, meta, order, "short", "max_rss_mb", int)
                  + build_groups(by_file, meta, order, "long", "max_rss_mb", int))
    grouped_bar_svg(
        os.path.join(FIG_DIR, "fig3_memory_comparison.svg"),
        "Peak memory: RastQC vs FastQC and Falco",
        "Peak resident set size (MB)",
        mem_groups,
        [FASTQC, FALCO, RASTQC_T1, RASTQC_TN],
        "lower is better",
        decimals=0,
    )

    md_tables(by_file, meta, order)
