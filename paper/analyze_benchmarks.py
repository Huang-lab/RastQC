#!/usr/bin/env python3
"""
Generate SVG figures and markdown tables from RastQC benchmark results.

Reads: benchmark/results/benchmark_results.csv (produced by benchmark/run_benchmark.sh)
Outputs: paper/figures/fig1_short_read_speed.svg
         paper/figures/fig2_long_read_speed.svg
         paper/figures/fig3_memory_comparison.svg
"""

import csv
import os
import sys

PAPER_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(PAPER_DIR)
CSV_PATH = os.path.join(PROJECT_DIR, "benchmark", "results", "benchmark_results.csv")
FIG_DIR = os.path.join(PAPER_DIR, "figures")


def load_results():
    rows = []
    with open(CSV_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


# ─── Figure 1: Short-Read Speed Comparison ──────────────────────────────────

def fig1_short_read_speed(rows):
    short = [r for r in rows if r["type"] == "short" and r["file"] != "ALL_SHORT"]
    files = []
    seen = set()
    for r in short:
        if r["file"] not in seen:
            files.append(r["file"])
            seen.add(r["file"])

    labels = []
    fastqc_times = []
    rastqc_times = []
    for fname in files:
        fqc = [r for r in short if r["file"] == fname and r["tool"] == "fastqc"]
        rqc = [r for r in short if r["file"] == fname and r["tool"] == "rastqc"]
        if fqc and rqc:
            label = fname.replace(".fastq.gz", "").replace("_", " ")
            if len(label) > 15:
                label = label[:12] + "..."
            labels.append(label)
            fastqc_times.append(float(fqc[0]["real_sec"]))
            rastqc_times.append(float(rqc[0]["real_sec"]))

    max_time = max(max(fastqc_times), max(rastqc_times)) * 1.15

    w, h = 700, 400
    ml, mr, mt, mb = 80, 120, 40, 80
    pw = w - ml - mr
    ph = h - mt - mb
    n = len(labels)
    gw = pw / n
    bw = gw * 0.35

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')
    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Short-Read Performance: RastQC vs FastQC</text>')

    # Y axis
    y_max = int(max_time / 10 + 1) * 10
    y_ticks = list(range(0, y_max + 1, max(y_max // 5, 1)))
    for t in y_ticks:
        y = mt + ph - (t / y_max * ph)
        svg.append(f'<line x1="{ml}" y1="{y}" x2="{ml + pw}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{ml - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}s</text>')

    cr, cf = "#2980b9", "#e74c3c"
    for i in range(n):
        xc = ml + gw * i + gw / 2
        xr = xc - bw - 2
        xf = xc + 2
        yb = mt + ph

        hr = (rastqc_times[i] / y_max) * ph
        hf = (fastqc_times[i] / y_max) * ph

        svg.append(f'<rect x="{xr}" y="{yb - hr}" width="{bw}" height="{hr}" fill="{cr}" rx="2"/>')
        svg.append(f'<rect x="{xf}" y="{yb - hf}" width="{bw}" height="{hf}" fill="{cf}" rx="2"/>')

        svg.append(f'<text x="{xr + bw/2}" y="{yb - hr - 4}" text-anchor="middle" font-size="9" fill="{cr}">{rastqc_times[i]:.1f}s</text>')
        svg.append(f'<text x="{xf + bw/2}" y="{yb - hf - 4}" text-anchor="middle" font-size="9" fill="{cf}">{fastqc_times[i]:.1f}s</text>')

        speedup = fastqc_times[i] / rastqc_times[i]
        svg.append(f'<text x="{xc}" y="{yb + 18}" text-anchor="middle" font-size="10" fill="#333">{labels[i]}</text>')
        svg.append(f'<text x="{xc}" y="{yb + 32}" text-anchor="middle" font-size="9" fill="#888">{speedup:.1f}x</text>')

    # Axes
    svg.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<text x="18" y="{mt + ph/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {mt + ph/2})">Wall-Clock Time (seconds)</text>')

    # Legend
    lx = w - mr + 10
    ly = mt + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{cr}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="11" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{cf}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="11" fill="#333">FastQC</text>')

    svg.append('</svg>')

    path = os.path.join(FIG_DIR, "fig1_short_read_speed.svg")
    with open(path, "w") as f:
        f.write("\n".join(svg))
    print(f"  Figure 1 saved: {path}")


# ─── Figure 2: Long-Read Speed Comparison ───────────────────────────────────

def fig2_long_read_speed(rows):
    long_rows = [r for r in rows if r["type"] == "long"]
    files = []
    seen = set()
    for r in long_rows:
        if r["file"] not in seen:
            files.append(r["file"])
            seen.add(r["file"])

    labels = []
    fastqc_t = []
    rastqc_t = []
    rastqc_lr_t = []
    for fname in files:
        fqc = [r for r in long_rows if r["file"] == fname and r["tool"] == "fastqc"]
        rqc = [r for r in long_rows if r["file"] == fname and r["tool"] == "rastqc"]
        rlr = [r for r in long_rows if r["file"] == fname and r["tool"] == "rastqc_lr"]
        if fqc and rqc:
            label = "ONT MinION" if "ont" in fname else "PacBio Revio"
            labels.append(label)
            fastqc_t.append(float(fqc[0]["real_sec"]))
            rastqc_t.append(float(rqc[0]["real_sec"]))
            rastqc_lr_t.append(float(rlr[0]["real_sec"]) if rlr else 0)

    max_time = max(max(fastqc_t), max(rastqc_t)) * 1.15

    w, h = 500, 400
    ml, mr, mt, mb = 80, 140, 40, 80
    pw = w - ml - mr
    ph = h - mt - mb
    n = len(labels)
    gw = pw / n
    bw = gw * 0.25

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')
    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Long-Read Performance: RastQC vs FastQC</text>')

    y_max = int(max_time / 5 + 1) * 5
    y_ticks = list(range(0, y_max + 1, max(y_max // 4, 1)))
    for t in y_ticks:
        y = mt + ph - (t / y_max * ph)
        svg.append(f'<line x1="{ml}" y1="{y}" x2="{ml + pw}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{ml - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}s</text>')

    cr, cf, clr = "#2980b9", "#e74c3c", "#27ae60"
    for i in range(n):
        xc = ml + gw * i + gw / 2
        xr = xc - bw * 1.5 - 2
        xf = xc - bw * 0.5
        xlr = xc + bw * 0.5 + 2
        yb = mt + ph

        hr = (rastqc_t[i] / y_max) * ph
        hf = (fastqc_t[i] / y_max) * ph
        hlr = (rastqc_lr_t[i] / y_max) * ph if rastqc_lr_t[i] > 0 else 0

        svg.append(f'<rect x="{xr}" y="{yb - hr}" width="{bw}" height="{hr}" fill="{cr}" rx="2"/>')
        svg.append(f'<rect x="{xf}" y="{yb - hf}" width="{bw}" height="{hf}" fill="{cf}" rx="2"/>')
        if hlr > 0:
            svg.append(f'<rect x="{xlr}" y="{yb - hlr}" width="{bw}" height="{hlr}" fill="{clr}" rx="2"/>')

        svg.append(f'<text x="{xr + bw/2}" y="{yb - hr - 4}" text-anchor="middle" font-size="8" fill="{cr}">{rastqc_t[i]:.1f}s</text>')
        svg.append(f'<text x="{xf + bw/2}" y="{yb - hf - 4}" text-anchor="middle" font-size="8" fill="{cf}">{fastqc_t[i]:.1f}s</text>')
        if rastqc_lr_t[i] > 0:
            svg.append(f'<text x="{xlr + bw/2}" y="{yb - hlr - 4}" text-anchor="middle" font-size="8" fill="{clr}">{rastqc_lr_t[i]:.1f}s</text>')

        speedup = fastqc_t[i] / rastqc_t[i]
        svg.append(f'<text x="{xc}" y="{yb + 18}" text-anchor="middle" font-size="11" fill="#333">{labels[i]}</text>')
        svg.append(f'<text x="{xc}" y="{yb + 32}" text-anchor="middle" font-size="9" fill="#888">{speedup:.1f}x faster</text>')

    svg.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<text x="18" y="{mt + ph/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {mt + ph/2})">Wall-Clock Time (seconds)</text>')

    lx = w - mr + 10
    ly = mt + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{cr}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="10" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{cf}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="10" fill="#333">FastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 44}" width="14" height="14" fill="{clr}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 56}" font-size="10" fill="#333">RastQC +LR</text>')

    svg.append('</svg>')

    path = os.path.join(FIG_DIR, "fig2_long_read_speed.svg")
    with open(path, "w") as f:
        f.write("\n".join(svg))
    print(f"  Figure 2 saved: {path}")


# ─── Figure 3: Memory Comparison ────────────────────────────────────────────

def fig3_memory(rows):
    # All per-file results (skip ALL_SHORT)
    files_data = []
    seen = set()
    for r in rows:
        if r["file"].startswith("ALL"):
            continue
        if r["file"] not in seen:
            seen.add(r["file"])
            fqc = [x for x in rows if x["file"] == r["file"] and x["tool"] == "fastqc"]
            rqc = [x for x in rows if x["file"] == r["file"] and x["tool"] == "rastqc"]
            if fqc and rqc:
                label = r["file"].replace(".fastq.gz", "").replace("_", " ")
                if len(label) > 18:
                    label = label[:15] + "..."
                files_data.append({
                    "label": label,
                    "type": r["type"],
                    "fqc_rss": int(fqc[0]["max_rss_mb"]),
                    "rqc_rss": int(rqc[0]["max_rss_mb"]),
                })

    n = len(files_data)
    max_mem = max(max(d["fqc_rss"] for d in files_data), max(d["rqc_rss"] for d in files_data)) * 1.15

    w, h = 800, 400
    ml, mr, mt, mb = 80, 120, 40, 90
    pw = w - ml - mr
    ph = h - mt - mb
    gw = pw / n
    bw = gw * 0.35

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" font-family="Arial, sans-serif">')
    svg.append(f'<rect width="{w}" height="{h}" fill="white"/>')
    svg.append(f'<text x="{w/2}" y="22" text-anchor="middle" font-size="14" font-weight="bold">Peak Memory (RSS): RastQC vs FastQC</text>')

    y_max = int(max_mem / 200 + 1) * 200
    for t in range(0, y_max + 1, max(y_max // 5, 1)):
        y = mt + ph - (t / y_max * ph)
        svg.append(f'<line x1="{ml}" y1="{y}" x2="{ml + pw}" y2="{y}" stroke="#e0e0e0" stroke-width="0.5"/>')
        svg.append(f'<text x="{ml - 8}" y="{y + 4}" text-anchor="end" font-size="11" fill="#555">{t}</text>')

    cr, cf = "#2980b9", "#e74c3c"
    for i, d in enumerate(files_data):
        xc = ml + gw * i + gw / 2
        xr = xc - bw - 2
        xf = xc + 2
        yb = mt + ph

        hr = (d["rqc_rss"] / y_max) * ph
        hf = (d["fqc_rss"] / y_max) * ph

        svg.append(f'<rect x="{xr}" y="{yb - hr}" width="{bw}" height="{hr}" fill="{cr}" rx="2"/>')
        svg.append(f'<rect x="{xf}" y="{yb - hf}" width="{bw}" height="{hf}" fill="{cf}" rx="2"/>')

        svg.append(f'<text x="{xr + bw/2}" y="{yb - hr - 4}" text-anchor="middle" font-size="8" fill="{cr}">{d["rqc_rss"]}</text>')
        svg.append(f'<text x="{xf + bw/2}" y="{yb - hf - 4}" text-anchor="middle" font-size="8" fill="{cf}">{d["fqc_rss"]}</text>')

        svg.append(f'<text x="{xc}" y="{yb + 16}" text-anchor="middle" font-size="9" fill="#333">{d["label"]}</text>')
        tag = "(LR)" if d["type"] == "long" else ""
        if tag:
            svg.append(f'<text x="{xc}" y="{yb + 28}" text-anchor="middle" font-size="8" fill="#888">{tag}</text>')

    svg.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<line x1="{ml}" y1="{mt + ph}" x2="{ml + pw}" y2="{mt + ph}" stroke="#333" stroke-width="1.5"/>')
    svg.append(f'<text x="18" y="{mt + ph/2}" text-anchor="middle" font-size="12" fill="#333" transform="rotate(-90 18 {mt + ph/2})">Memory (MB)</text>')

    lx = w - mr + 10
    ly = mt + 20
    svg.append(f'<rect x="{lx}" y="{ly}" width="14" height="14" fill="{cr}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 12}" font-size="11" fill="#333">RastQC</text>')
    svg.append(f'<rect x="{lx}" y="{ly + 22}" width="14" height="14" fill="{cf}" rx="2"/>')
    svg.append(f'<text x="{lx + 20}" y="{ly + 34}" font-size="11" fill="#333">FastQC</text>')

    svg.append('</svg>')

    path = os.path.join(FIG_DIR, "fig3_memory_comparison.svg")
    with open(path, "w") as f:
        f.write("\n".join(svg))
    print(f"  Figure 3 saved: {path}")


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if not os.path.exists(CSV_PATH):
        print(f"ERROR: Benchmark results not found at {CSV_PATH}")
        print("Run benchmark/run_benchmark.sh first.")
        sys.exit(1)

    os.makedirs(FIG_DIR, exist_ok=True)
    rows = load_results()

    print("=" * 60)
    print("RastQC Benchmark Figure Generation")
    print(f"Source: {CSV_PATH}")
    print("=" * 60)

    # Print summary tables
    print("\n━━━ SHORT-READ RESULTS ━━━")
    print(f"{'File':<30} {'Size':>6} {'FastQC':>9} {'RastQC':>9} {'Speedup':>8} {'FQC MB':>7} {'RQC MB':>7}")
    for r in rows:
        if r["type"] == "short" and r["tool"] == "fastqc":
            fqc = r
            rqc = [x for x in rows if x["file"] == r["file"] and x["tool"] == "rastqc"]
            if rqc:
                rqc = rqc[0]
                ft = float(fqc["real_sec"])
                rt = float(rqc["real_sec"])
                sp = ft / rt if rt > 0 else 0
                print(f"{r['file']:<30} {r['size_mb']:>5}M {ft:>8.2f}s {rt:>8.2f}s {sp:>7.1f}x {fqc['max_rss_mb']:>6}M {rqc['max_rss_mb']:>6}M")

    long_rows = [r for r in rows if r["type"] == "long"]
    if long_rows:
        print("\n━━━ LONG-READ RESULTS ━━━")
        print(f"{'File':<35} {'FastQC':>9} {'RastQC':>9} {'RQC+LR':>9} {'Speedup':>8}")
        seen = set()
        for r in long_rows:
            if r["file"] in seen:
                continue
            seen.add(r["file"])
            fqc = [x for x in long_rows if x["file"] == r["file"] and x["tool"] == "fastqc"]
            rqc = [x for x in long_rows if x["file"] == r["file"] and x["tool"] == "rastqc"]
            rlr = [x for x in long_rows if x["file"] == r["file"] and x["tool"] == "rastqc_lr"]
            if fqc and rqc:
                ft = float(fqc[0]["real_sec"])
                rt = float(rqc[0]["real_sec"])
                lt = float(rlr[0]["real_sec"]) if rlr else 0
                sp = ft / rt if rt > 0 else 0
                lt_str = f"{lt:>8.2f}s" if lt > 0 else "       -"
                print(f"{r['file']:<35} {ft:>8.2f}s {rt:>8.2f}s {lt_str} {sp:>7.1f}x")

    print("\nGenerating figures...")
    fig1_short_read_speed(rows)
    fig2_long_read_speed(rows)
    fig3_memory(rows)
    print("\nDone.")
