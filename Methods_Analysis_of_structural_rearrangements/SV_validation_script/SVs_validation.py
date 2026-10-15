#!/usr/bin/env python3
"""
make_sv_cov_panels.py: coverage panels and a coverage-shift test for SVs.

For each SV, the script builds fixed-bin windows around the breakpoints, computes
per-window mean depth with mosdepth, normalizes each sample by its genome-wide mean
depth, plots all samples, and (for real intra-chromosomal SVs only) tests whether
the target-vs-background log2 coverage ratio inside the event differs from the
flanks, using an exact circular-shift permutation test.

See README.md for full documentation.
"""

import argparse
import csv
import gzip
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

__version__ = "2.0.0"

# ------------------------------- Defaults -------------------------------------
DEFAULT_INTER_FLANK_BP = 1000   # flank around each inter-chr breakend
DEFAULT_INTRA_FLANK_PCT = 0.20  # intra-chr flank as fraction of SV length
DEFAULT_INTRA_MIN_FLANK_BP = 400
DEFAULT_BINS_PER_BLOCK = 50
DEFAULT_THREADS = 4
EPS = 1e-9                      # pseudocount for log ratios
TIE_TOL = 1e-12                 # null values within this of |obs| count as ties (conservative)

Block = tuple[str, int, int]    # (chrom, start, end), 0-based half-open


# ------------------------------- Messages -------------------------------------
def warn(msg: str) -> None:
    print(f"[warn] {msg}", file=sys.stderr)


def info(msg: str) -> None:
    print(f"[info] {msg}", file=sys.stderr)


class SVSkip(ValueError):
    """Raised for an SV that cannot be processed; the SV is skipped with a warning"""


# ------------------------------- Small utilities ------------------------------
def expand_abs_path(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def sanitize_sv_id(sv_id: str) -> str:
    return re.sub(r"[|:/\\ ]", "_", sv_id)


def is_like_sv(sv_id: str) -> bool:
    """SVs whose ID contains 'like' (case-insensitive) are plotted but never tested"""
    return "like" in str(sv_id).lower()


# ------------------------------- SV ID parsing --------------------------------
# Two-breakend form:  TYPE|chrA:posA-chrB:posB   (chrA may equal chrB)
# Span form:          TYPE|chr:start-end
_ID_TWO = re.compile(r"^[^|]+\|([^:|]+):(\d+)-([^:|]+):(\d+)")
_ID_ONE = re.compile(r"^[^|]+\|([^:|]+):(\d+)-(\d+)")


def parse_sv_coords(sv_id: str) -> tuple[str, int, str, int]:
    m = _ID_TWO.match(sv_id)
    if m:
        return m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
    m = _ID_ONE.match(sv_id)
    if m:
        return m.group(1), int(m.group(2)), m.group(1), int(m.group(3))
    raise SVSkip(f"cannot parse INFO_variantID '{sv_id}'")


def compute_flank(c1: str, p1: int, c2: str, p2: int, intra_pct: float,
                  intra_min_bp: int, intra_max_bp: int | None,
                  inter_bp: int) -> tuple[int, bool]:
    """Return (flank_bp, is_inter_chrom)"""
    if c1 != c2:
        return inter_bp, True
    f = max(int(round(abs(p2 - p1) * intra_pct)), intra_min_bp)
    if intra_max_bp is not None:
        f = min(f, intra_max_bp)
    return max(1, f), False


def clamp_block(chrom: str, start: int, end: int, chrom_sizes: dict[str, int]) -> Block:
    if chrom not in chrom_sizes:
        raise SVSkip(f"contig '{chrom}' not found in the BAM header")
    s = max(0, int(start))
    e = min(int(end), chrom_sizes[chrom])
    return (chrom, s, max(s, e))


def build_blocks(c1: str, p1: int, c2: str, p2: int, flank: int,
                 chrom_sizes: dict[str, int]) -> tuple[list[Block], list[dict]]:
    if c1 == c2:
        lo, hi = sorted((p1, p2))
        blocks = [clamp_block(c1, lo - flank, hi + flank, chrom_sizes)]
        bps = [{"block_idx": 0, "pos": lo}, {"block_idx": 0, "pos": hi}]
    else:
        blocks = [clamp_block(c1, p1 - flank, p1 + flank, chrom_sizes),
                  clamp_block(c2, p2 - flank, p2 + flank, chrom_sizes)]
        bps = [{"block_idx": 0, "pos": p1}, {"block_idx": 1, "pos": p2}]
    return blocks, bps


# ------------------------------- Windowing ------------------------------------
def make_windows_df(blocks: list[Block], bins_per_block: int) -> pd.DataFrame:
    rows = []
    for bi, (chrom, start, end) in enumerate(blocks):
        if end <= start:
            continue
        edges = np.maximum.accumulate(np.linspace(start, end, bins_per_block + 1, dtype=int))
        for wi in range(bins_per_block):
            beg, fin = int(edges[wi]), int(edges[wi + 1])
            if fin > beg:
                rows.append((chrom, beg, fin, bi, wi, f"b{bi}w{wi}"))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "block_idx", "win_idx", "name"])
    df["xidx"] = np.arange(len(df))
    return df


def write_windows_bed(windows_df: pd.DataFrame, bed_path: Path) -> None:
    windows_df[["chrom", "start", "end", "name"]].to_csv(
        bed_path, sep="\t", header=False, index=False)


# ------------------------------- mosdepth + caching ---------------------------
def run_mosdepth(prefix: str, bam: str, threads: int, bed: str | None = None) -> None:
    cmd = ["mosdepth", "--threads", str(threads), "--no-per-base", "--fast-mode"]
    if bed is not None:
        cmd += ["--by", bed]
    cmd += [prefix, bam]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"mosdepth failed (exit {p.returncode})\nCMD: {' '.join(cmd)}\nSTDERR:\n{p.stderr}")


def _stamp_path(prefix: Path) -> Path:
    return Path(f"{prefix}.cache.json")


def _bam_fingerprint(bam: str) -> dict:
    st = os.stat(bam)
    return {"bam": bam, "size": st.st_size, "mtime_ns": st.st_mtime_ns}


def cache_is_valid(prefix: Path, bam: str, outputs: list[Path]) -> bool:
    """Valid only if all outputs exist and were produced from this exact BAM file."""
    stamp = _stamp_path(prefix)
    if not stamp.exists() or not all(o.exists() for o in outputs):
        return False
    try:
        return json.loads(stamp.read_text()) == _bam_fingerprint(bam)
    except (json.JSONDecodeError, OSError):
        return False


def write_stamp(prefix: Path, bam: str) -> None:
    _stamp_path(prefix).write_text(json.dumps(_bam_fingerprint(bam)))


def get_wg_summary(sample_id: str, bam: str, mosdepth_dir: Path, threads: int) -> Path:
    prefix = mosdepth_dir / sample_id
    summary = Path(f"{prefix}.mosdepth.summary.txt")
    if not cache_is_valid(prefix, bam, [summary]):
        run_mosdepth(str(prefix), bam, threads)
        write_stamp(prefix, bam)
    return summary


def _read_summary_rows(summary_path: Path) -> list[tuple[str, int, float]]:
    rows = []
    with open(summary_path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)  # header
        for row in rdr:
            if row and not row[0].startswith("#"):
                rows.append((row[0], int(row[1]), float(row[3])))
    return rows


def chrom_sizes_from_summary(summary_path: Path) -> dict[str, int]:
    sizes = {c: L for c, L, _ in _read_summary_rows(summary_path)
             if c.lower() != "total" and not c.endswith("_region")}
    if not sizes:
        raise RuntimeError(f"No contig sizes in mosdepth summary: {summary_path}")
    return sizes


def wg_mean_from_summary(summary_path: Path) -> float:
    rows = _read_summary_rows(summary_path)
    for c, _, m in rows:
        if c.lower() == "total":
            return m
    contigs = [(L, m) for c, L, m in rows if not c.endswith("_region")]
    if not contigs:
        raise RuntimeError(f"No rows in mosdepth summary: {summary_path}")
    total_len = sum(L for L, _ in contigs)
    return sum(L * m for L, m in contigs) / max(1, total_len)


def _match_regions(regions_gz: Path, windows_df: pd.DataFrame) -> np.ndarray | None:
    """
    Read mosdepth regions keyed by window name and verify coordinates.
    Returns means in windows_df order, or None if any window is missing or differs. 
    """
    found: dict[str, tuple[str, int, int, float]] = {}
    with gzip.open(regions_gz, "rt") as fh:
        for ln in fh:
            if ln.startswith("#") or not ln.strip():
                continue
            f = ln.rstrip("\n").split("\t")
            found[f[3]] = (f[0], int(f[1]), int(f[2]), float(f[4]))
    vals = np.empty(len(windows_df), dtype=float)
    for i, r in enumerate(windows_df.itertuples(index=False)):
        hit = found.get(r.name)
        if hit is None or hit[:3] != (r.chrom, int(r.start), int(r.end)):
            return None
        vals[i] = hit[3]
    return vals


def get_window_coverage(sv_tag: str, sid: str, bam: str, bed_path: Path,
                        windows_df: pd.DataFrame, mosdepth_dir: Path,
                        threads: int) -> np.ndarray:
    prefix = mosdepth_dir / f"{sv_tag}__{sid}"
    regions = Path(f"{prefix}.regions.bed.gz")
    for attempt in range(2):
        if attempt == 1 or not cache_is_valid(prefix, bam, [regions]):
            run_mosdepth(str(prefix), bam, threads, bed=str(bed_path))
            write_stamp(prefix, bam)
        vals = _match_regions(regions, windows_df)
        if vals is not None:
            return vals
    raise RuntimeError(f"mosdepth regions do not match windows for {sid} / {sv_tag}")


# ------------------------------- Statistics -----------------------------------
def compute_log2fc_target_vs_bg(mat: np.ndarray) -> np.ndarray:
    """Row 0 is the target; background = per-window median of the other rows.
    With no background samples, the target is compared with 1 (genome-wide mean)"""
    target = mat[0, :]
    bg = np.nanmedian(mat[1:, :], axis=0) if mat.shape[0] > 1 else np.ones_like(target)
    return np.log2((target + EPS) / (bg + EPS))


def masks_for_intra(block: Block, bps: list[dict], n_bins: int,
                    flank_gap_bins: int) -> tuple[np.ndarray, np.ndarray] | None:
    if len(bps) != 2 or n_bins == 0:
        return None
    _, start, end = block
    span = max(1, end - start)

    def xpos(pos: int) -> float:
        return float(np.clip((pos - start) / span, 0.0, 1.0)) * (n_bins - 1)

    x0, x1 = sorted(xpos(bp["pos"]) for bp in bps)
    left_i, right_i = int(np.floor(x0)), int(np.ceil(x1))
    idx = np.arange(n_bins)
    m_event = (idx >= left_i) & (idx <= right_i)

    left_cand = idx[idx <= left_i - flank_gap_bins - 1]
    right_cand = idx[idx >= right_i + flank_gap_bins + 1]

    # Take up to as many flank bins as event bins; if one side is short,
    # compensate from the other side.
    ev_len = int(m_event.sum())
    left_take = min(ev_len, left_cand.size)
    right_take = min(ev_len, right_cand.size)
    if left_take < ev_len:
        right_take = min(2 * ev_len - left_take, right_cand.size)
    if right_take < ev_len:
        left_take = min(2 * ev_len - right_take, left_cand.size)

    m_flank = np.zeros(n_bins, dtype=bool)
    if left_take > 0:
        m_flank[left_cand[-left_take:]] = True
    if right_take > 0:
        m_flank[right_cand[:right_take]] = True
    return m_event, m_flank


def exact_circular_shift_test(effect_fn, values: np.ndarray):
    """
    Exact two-sided permutation test over all n-1 non-identity circular shifts.
    p = (#{shifts with |null| >= |obs|} + 1) / n, so the minimum achievable p is 1/n.
    Returns (obs, p, p_min, shifts, null).
    """
    n = values.size
    obs = effect_fn(values)
    if not np.isfinite(obs) or n < 2:
        return obs, np.nan, np.nan, np.array([], int), np.array([], float)
    shifts = np.arange(1, n)
    null = np.array([effect_fn(np.roll(values, s)) for s in shifts], dtype=float)
    b = int(np.sum(np.abs(null) >= np.abs(obs) - TIE_TOL))
    return float(obs), (b + 1) / n, 1.0 / n, shifts, null


def compute_sv_stats_intra(blocks, bps, windows_df, mat, flank_gap_bins, detrend_window):
    log2fc = compute_log2fc_target_vs_bg(mat)
    results = []
    for bi in sorted(windows_df["block_idx"].unique()):
        block_bps = [bp for bp in bps if bp["block_idx"] == bi]
        if len(block_bps) != 2:
            continue
        sub = windows_df[windows_df["block_idx"] == bi].sort_values("win_idx")
        series = log2fc[sub["xidx"].to_numpy()]
        masks = masks_for_intra(blocks[bi], block_bps, series.size, flank_gap_bins)
        if masks is None:
            continue
        m_event, m_flank = masks
        if detrend_window and detrend_window > 1:
            trend = pd.Series(series).rolling(int(detrend_window), center=True,
                                              min_periods=1).median().to_numpy()
            series = series - trend
        good = np.isfinite(series)
        if not (m_event & good).any() or not (m_flank & good).any():
            continue

        def eff(v: np.ndarray) -> float:
            return float(np.nanmedian(v[m_event]) - np.nanmedian(v[m_flank]))

        obs, p, p_min, shifts, null = exact_circular_shift_test(eff, series)
        results.append({
            "block_idx": int(bi),
            "mode": "intra",
            "effect_median_diff_log2": obs,
            "p_perm_two_sided": p,
            "p_min_achievable": p_min,
            "n_windows": int(series.size),
            "n_event_bins": int(m_event.sum()),
            "n_flank_bins": int(m_flank.sum()),
            "_shifts": shifts,
            "_null": null,
        })
    return results


def bh_adjust(pvals: pd.Series) -> pd.Series:
    p = pvals.astype(float).to_numpy()
    q = np.full_like(p, np.nan, dtype=float)
    idx = np.where(np.isfinite(p))[0]
    if idx.size:
        order = idx[np.argsort(p[idx], kind="mergesort")]
        ranks = np.arange(1, order.size + 1)
        tmp = p[order] * idx.size / ranks
        q[order] = np.clip(np.minimum.accumulate(tmp[::-1])[::-1], 0.0, 1.0)
    return pd.Series(q, index=pvals.index, name="p_adj_bh")


# ------------------------------- Plotting -------------------------------------
def smooth_per_block(mat: np.ndarray, block_idx: np.ndarray, k: int) -> np.ndarray:
    """Centered moving average within each block; no zero-padding, no cross-block mixing"""
    out = mat.copy()
    for b in np.unique(block_idx):
        cols = block_idx == b
        for i in range(mat.shape[0]):
            out[i, cols] = (pd.Series(mat[i, cols])
                            .rolling(k, center=True, min_periods=1).mean().to_numpy())
    return out


def plot_coverage_lines(out_svg, windows_df, blocks, bps, mat, sample_ids, title,
                        ymin=None, ymax=None, smooth=0, show_mid_ticks=False):
    n_samples, n_windows = mat.shape
    x = np.arange(n_windows)
    block_idx_arr = windows_df["block_idx"].to_numpy()
    if smooth and smooth > 1:
        mat = smooth_per_block(mat, block_idx_arr, int(smooth))

    fig = plt.figure(figsize=(max(8, n_windows / 12), 4 + 0.28 * n_samples))
    ax = plt.gca()
    for i, sid in enumerate(sample_ids):
        ax.plot(x, mat[i, :], label=sid, linewidth=1.2)

    for cut in np.where(block_idx_arr[1:] != block_idx_arr[:-1])[0] + 1:
        ax.axvline(cut - 0.5, linewidth=0.6, alpha=0.5)

    x_start = windows_df.groupby("block_idx")["xidx"].min().to_dict()
    n_bins = windows_df.groupby("block_idx").size().to_dict()

    for bi, (chrom, _, _) in enumerate(blocks):
        if bi in x_start:
            ax.text(x_start[bi] + (n_bins[bi] - 1) / 2.0, 1.02, chrom,
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=8, clip_on=False)

    for bp in bps:
        bi = bp["block_idx"]
        _, s, e = blocks[bi]
        if e <= s or bi not in x_start:
            continue
        frac = min(1.0, max(0.0, (bp["pos"] - s) / float(e - s)))
        ax.axvline(x_start[bi] + frac * (n_bins[bi] - 1), linestyle="--", linewidth=1.0)

    is_inter = len(blocks) == 2 and blocks[0][0] != blocks[1][0]
    xticks, labels, rots = [], [], []
    for bi, (_, s, e) in enumerate(blocks):
        if bi not in x_start:
            continue
        x0, nb = x_start[bi], n_bins[bi]
        rot = (45 if bi == 0 else -45) if is_inter else 45
        ticks = [(x0, s)]
        if show_mid_ticks:
            ticks.append((x0 + (nb - 1) / 2.0, (s + e) // 2))
        ticks.append((x0 + nb - 1, e))
        for xt, lab in ticks:
            xticks.append(xt); labels.append(f"{lab:,}"); rots.append(rot)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, fontsize=8)
    for lab, rot in zip(ax.get_xticklabels(), rots):
        lab.set_rotation(rot)
        lab.set_ha("right" if rot > 0 else "left")

    ax.set_xlim(0, n_windows - 1)
    if ymin is not None or ymax is not None:
        lo, hi = ax.get_ylim()
        ax.set_ylim(ymin if ymin is not None else lo, ymax if ymax is not None else hi)
    ax.set_xlabel("Position (bp)", labelpad=2)
    ax.set_ylabel("Relative coverage (genome-wide mean = 1)")
    ax.set_title(title, fontsize=10, pad=25, loc="center")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), fontsize=8,
              frameon=False, ncol=1, borderaxespad=0)
    plt.tight_layout(rect=(0, 0, 0.82, 0.95))
    fig.savefig(out_svg)
    plt.close(fig)


# ------------------------------- PDF ----------------------------
def build_pdf(svg_paths: list[Path], out_pdf: Path, margin_mm: float = 6.0) -> None:
    """Portrait A4, two figures per page (top/bottom), in the given order"""
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import mm
        from reportlab.pdfgen import canvas
        from reportlab.graphics import renderPDF
        from svglib.svglib import svg2rlg
    except ImportError:
        warn("reportlab/svglib not installed; skipping PDF generation")
        return

    W, H = A4
    m = margin_mm * mm
    cell_w, cell_h = W - 2 * m, (H - 2 * m) / 2
    c = canvas.Canvas(str(out_pdf), pagesize=(W, H))
    for i, p in enumerate(svg_paths):
        if i and i % 2 == 0:
            c.showPage()
        d = svg2rlg(str(p))
        if d is None or not getattr(d, "width", 0) or not getattr(d, "height", 0):
            warn(f"could not render {p} into PDF; skipped")
            continue
        s = min(cell_w / d.width, cell_h / d.height)
        d.scale(s, s)
        y = m + (1 - i % 2) * cell_h
        renderPDF.draw(d, c, m + (cell_w - d.width * s) / 2, y + (cell_h - d.height * s) / 2)
    c.save()


# ------------------------------- CLI ------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Coverage panels and an exact circular-shift coverage test for SVs."
                    "See README.md for details.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    io = p.add_argument_group("inputs/outputs")
    io.add_argument("--sv_table", required=True,
                    help="SVs_changing.tab (tab-separated; see README for columns)")
    io.add_argument("--paths_table", required=True,
                    help="Tab-separated table with columns sampleID, sorted_bam")
    io.add_argument("--outdir", required=True, help="Output directory (created if missing)")
    io.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="mosdepth threads")
    io.add_argument("--all_targets", action="store_true",
                    help="Process every (INFO_variantID, target_sampleID) pair. By default "
                         "only the first target per INFO_variantID is kept")
    io.add_argument("--no_pdf", action="store_true", help="Do not build the PDF")

    fl = p.add_argument_group("windows")
    fl.add_argument("--intra_flank_pct", type=float, default=DEFAULT_INTRA_FLANK_PCT,
                    help="Intra-chr flank as a fraction of SV length")
    fl.add_argument("--intra_min_flank_bp", type=int, default=DEFAULT_INTRA_MIN_FLANK_BP,
                    help="Minimum intra-chr flank (bp)")
    fl.add_argument("--intra_max_flank_bp", type=int, default=-1,
                    help="Maximum intra-chr flank (bp); -1 disables the cap")
    fl.add_argument("--inter_flank_bp", type=int, default=DEFAULT_INTER_FLANK_BP,
                    help="Flank around each inter-chr breakend (bp)")
    fl.add_argument("--bins_per_block", type=int, default=DEFAULT_BINS_PER_BLOCK,
                    help="Bins per block. Also sets test resolution: minimum p = 1/bins")

    st = p.add_argument_group("statistics")
    st.add_argument("--flank_gap_bins", type=int, default=1,
                    help="Bins just outside the event excluded from the flanks")
    st.add_argument("--detrend_window", type=int, default=0,
                    help="If >1, subtract a rolling-median trend of this many bins before testing")

    pl = p.add_argument_group("plotting")
    pl.add_argument("--smooth", type=int, default=0,
                    help="Moving-average window in bins (plots only; 0 = off)")
    pl.add_argument("--ymin", type=float, default=None, help="Y-axis minimum")
    pl.add_argument("--ymax", type=float, default=None, help="Y-axis maximum")
    pl.add_argument("--show_mid", action="store_true", help="Add a midpoint tick per block")

    a = p.parse_args()
    if a.bins_per_block < 5:
        p.error("--bins_per_block must be >= 5")
    if a.intra_max_flank_bp != -1 and a.intra_max_flank_bp < 1:
        p.error("--intra_max_flank_bp must be >= 1, or -1 to disable")
    if a.intra_min_flank_bp < 1 or a.inter_flank_bp < 1:
        p.error("flank sizes must be >= 1 bp")
    return a


def load_tables(args):
    sv_df = pd.read_csv(expand_abs_path(args.sv_table), sep="\t")
    missing = {"INFO_variantID", "target_sampleID", "type_background"} - set(sv_df.columns)
    if missing:
        sys.exit(f"SV table is missing required columns: {sorted(missing)}")
    sv_df = sv_df[sv_df["type_background"].astype(str).str.lower() == "low_confidence"]
    key = ["INFO_variantID", "target_sampleID"] if args.all_targets else ["INFO_variantID"]
    sv_df = sv_df.drop_duplicates(subset=key, keep="first").reset_index(drop=True)

    paths_df = pd.read_csv(expand_abs_path(args.paths_table), sep="\t")
    missing = {"sampleID", "sorted_bam"} - set(paths_df.columns)
    if missing:
        sys.exit(f"paths table is missing required columns: {sorted(missing)}")
    if paths_df.empty:
        sys.exit("paths table has no rows")
    paths_df["sampleID"] = paths_df["sampleID"].astype(str)
    dups = paths_df["sampleID"][paths_df["sampleID"].duplicated()].unique()
    if len(dups):
        sys.exit(f"duplicated sampleIDs in paths table: {list(dups)}")
    paths_df["sorted_bam"] = paths_df["sorted_bam"].map(expand_abs_path)
    for bam in paths_df["sorted_bam"]:
        if not os.path.exists(bam):
            sys.exit(f"BAM not found: {bam}")
    return sv_df, paths_df.set_index("sampleID")


# ------------------------------- Main -----------------------------------------
def main() -> None:
    args = parse_args()
    sv_df, paths_df = load_tables(args)
    if sv_df.empty:
        info("no SVs with type_background == low_confidence; nothing to do")
        return

    outdir = Path(expand_abs_path(args.outdir))
    bed_dir = outdir / "beds"
    mosdepth_dir = outdir / "mosdepth"
    figs_real, figs_like = outdir / "figs" / "real", outdir / "figs" / "like"
    stats_dir = outdir / "stats"
    cov_dir = stats_dir / "coverage_stats"
    by_sample_dir = cov_dir / "by_sample"
    perms_dir = stats_dir / "permutations"
    for d in (bed_dir, mosdepth_dir, figs_real, figs_like, by_sample_dir, perms_dir):
        d.mkdir(parents=True, exist_ok=True)

    intra_max = None if args.intra_max_flank_bp == -1 else args.intra_max_flank_bp

    # Contig sizes from the genome-wide summary of the first sample.
    first = paths_df.index[0]
    chrom_sizes = chrom_sizes_from_summary(
        get_wg_summary(first, paths_df.loc[first, "sorted_bam"], mosdepth_dir, args.threads))

    wg_means: dict[str, float] = {}
    per_sample_rows: defaultdict[str, list[dict]] = defaultdict(list)
    all_rows: list[dict] = []
    real_figs: list[str] = []
    n_skipped = 0

    for _, row in sv_df.iterrows():
        sv_id = str(row["INFO_variantID"])
        target_id = str(row["target_sampleID"])
        diff_tag = row.get("type_differetial_variant", "NA")
        diff_tag = "NA" if pd.isna(diff_tag) else str(diff_tag)
        bg_raw = row.get("background_sampleIDs", "")
        bg_list = ([s.strip() for s in str(bg_raw).split(",") if s.strip()]
                   if pd.notna(bg_raw) else [])

        if target_id not in paths_df.index:
            warn(f"{sv_id}: target sample '{target_id}' not in paths table; SV skipped")
            n_skipped += 1
            continue

        try:
            c1, p1, c2, p2 = parse_sv_coords(sv_id)
            flank, is_inter = compute_flank(c1, p1, c2, p2, args.intra_flank_pct,
                                            args.intra_min_flank_bp, intra_max,
                                            args.inter_flank_bp)
            blocks, bps = build_blocks(c1, p1, c2, p2, flank, chrom_sizes)
        except SVSkip as e:
            warn(f"{e}; SV skipped")
            n_skipped += 1
            continue

        windows_df = make_windows_df(blocks, args.bins_per_block)
        if windows_df.empty:
            warn(f"{sv_id}: no valid windows after clamping; SV skipped")
            n_skipped += 1
            continue

        sv_tag = sanitize_sv_id(sv_id)
        bed_path = bed_dir / f"{sv_tag}.windows.bed"
        write_windows_bed(windows_df, bed_path)

        sample_order = [target_id] + [s for s in bg_list if s != target_id]
        sample_ids, rows_norm = [], []
        for sid in sample_order:
            if sid not in paths_df.index:
                warn(f"{sv_id}: background sample '{sid}' not in paths table; ignored")
                continue
            bam = paths_df.loc[sid, "sorted_bam"]
            if sid not in wg_means:
                wg_means[sid] = wg_mean_from_summary(
                    get_wg_summary(sid, bam, mosdepth_dir, args.threads))
            wg = wg_means[sid]
            if not np.isfinite(wg) or wg <= 0:
                warn(f"invalid genome-wide mean for '{sid}'; sample ignored")
                continue
            vals = get_window_coverage(sv_tag, sid, bam, bed_path, windows_df,
                                       mosdepth_dir, args.threads)
            rows_norm.append(vals / wg)
            sample_ids.append(sid)

        # The target must be row 0; otherwise a background sample would be tested as the target.
        if not sample_ids or sample_ids[0] != target_id:
            warn(f"{sv_id}: target sample '{target_id}' has no usable coverage; SV skipped")
            n_skipped += 1
            continue
        mat = np.vstack(rows_norm)

        like = is_like_sv(sv_id)
        base = f"{sv_tag}_{target_id}_{diff_tag}"
        fig_path = (figs_like if like else figs_real) / f"{base}.coverage.svg"
        plot_coverage_lines(str(fig_path), windows_df, blocks, bps, mat, sample_ids, base,
                            ymin=args.ymin, ymax=args.ymax, smooth=args.smooth,
                            show_mid_ticks=args.show_mid)
        fig_rel = str(fig_path.relative_to(outdir))
        if not like:
            real_figs.append(fig_rel)

        if like or is_inter:
            continue
        stats = compute_sv_stats_intra(blocks, bps, windows_df, mat,
                                       args.flank_gap_bins, args.detrend_window)
        if not stats:
            continue

        public = [{k: v for k, v in s.items() if not k.startswith("_")} for s in stats]
        cov_path = cov_dir / f"{base}.coverage.stats.tsv"
        pd.DataFrame(public).to_csv(cov_path, sep="\t", index=False)

        long_rows = []
        for s, pub in zip(stats, public):
            rec = {"sv_id": sv_id, "target_sampleID": target_id,
                   "type_differetial_variant": diff_tag, **pub, "figure": fig_rel}
            per_sample_rows[target_id].append(rec)
            all_rows.append(rec | {"cov_stats_file": str(cov_path.relative_to(outdir))})
            long_rows.append({"sv_id": sv_id, "target_sampleID": target_id,
                              "block_idx": s["block_idx"], "kind": "obs", "shift": 0,
                              "effect": s["effect_median_diff_log2"]})
            long_rows += [{"sv_id": sv_id, "target_sampleID": target_id,
                           "block_idx": s["block_idx"], "kind": "perm", "shift": int(sh),
                           "effect": float(v)} for sh, v in zip(s["_shifts"], s["_null"])]
        pd.DataFrame(long_rows).to_csv(perms_dir / f"{base}.permutations.tsv",
                                       sep="\t", index=False)

    # Per-sample tables (sorted by raw p).
    for sample_id, rows in per_sample_rows.items():
        (pd.DataFrame(rows)
         .sort_values("p_perm_two_sided", kind="stable", na_position="last")
         .to_csv(by_sample_dir / f"{sample_id}.coverage.stats.tsv", sep="\t", index=False))

    # Aggregated table with BH adjustment; always written so no stale file survives.
    cols = ["sv_id", "target_sampleID", "type_differetial_variant", "block_idx", "mode",
            "effect_median_diff_log2", "p_perm_two_sided", "p_min_achievable",
            "n_windows", "n_event_bins", "n_flank_bins", "figure", "cov_stats_file", "p_adj_bh"]
    all_df = pd.DataFrame(all_rows)
    if all_df.empty:
        all_df = pd.DataFrame(columns=cols)
    else:
        all_df["p_adj_bh"] = bh_adjust(all_df["p_perm_two_sided"])
        all_df = all_df.sort_values(["p_adj_bh", "p_perm_two_sided"],
                                    na_position="last", kind="stable")[cols]
    all_df.to_csv(stats_dir / "coverage.stats.tsv", sep="\t", index=False)

    # PDF: tested real SVs in ranked order, then untested real SVs in input order.
    if not args.no_pdf:
        ranked = list(dict.fromkeys(all_df["figure"].astype(str)))
        ordered = ranked + [f for f in dict.fromkeys(real_figs) if f not in set(ranked)]
        svgs = [outdir / f for f in ordered if (outdir / f).exists()]
        if svgs:
            build_pdf(svgs, outdir / "figs" / "coverage_plots_REAL_SVs.pdf")
        else:
            info("no real-SV figures; PDF not built")

    n_sig = int((all_df["p_adj_bh"] < 0.05).sum()) if len(all_df) else 0
    info(f"done: {len(sv_df) - n_skipped} SVs processed, {n_skipped} skipped, "
         f"{len(all_df)} blocks tested, {n_sig} with BH-adjusted p < 0.05")


if __name__ == "__main__":
    main()