"""
Final Project: Motif Discovery & Repeated Substring Mining
"""

import os
import sys
import time
import tracemalloc
import collections
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend for servers
import matplotlib.pyplot as plt

# Configurable parameters
L_VALUES = [30, 70, 120]   # motif lengths to inspect
MIN_FREQ = 3               # minimum frequency for reporting motifs
MINIMAL_MAXIMAL_LEN = 100
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "output")
OUTPUT_DIRS = {
    'summary': os.path.join(OUTPUT_DIR, 'summaries'),
    'detailed': os.path.join(OUTPUT_DIR, 'detailed'),
    'figures': os.path.join(OUTPUT_DIR, 'figures'),
    'comparison': os.path.join(OUTPUT_DIR, 'comparison'),
}

for d in OUTPUT_DIRS.values():
    os.makedirs(d, exist_ok=True)


def read_sequence(path):
    """Read a sequence file (plain text/FASTA-ish) and return a cleaned uppercase sequence"""
    with open(path, 'r', encoding='utf-8') as f:
        txt = f.read()
    # drop FASTA headers and whitespace
    lines = [ln.strip() for ln in txt.splitlines() if ln and not ln.startswith(">")]
    seq = "".join(lines).replace(" ", "").replace("\r", "").upper()
    return seq


# ---------------------------
# Suffix array and LCP utils
# ---------------------------
def build_suffix_array(s):
    """
    Build suffix array using doubling algorithm: O(n log n)
    Returns list of start indices sorted by suffix.
    """
    n = len(s)
    if n == 0:
        return []
    k = 1
    sa = list(range(n))
    rank = [ord(c) for c in s]
    tmp = [0] * n
    while True:
        sa.sort(key=lambda x: (rank[x], rank[x + k] if x + k < n else -1))
        tmp[sa[0]] = 0
        for i in range(1, n):
            prev, cur = sa[i - 1], sa[i]
            prev_key = (rank[prev], rank[prev + k] if prev + k < n else -1)
            cur_key = (rank[cur], rank[cur + k] if cur + k < n else -1)
            tmp[cur] = tmp[prev] + (cur_key != prev_key)
        rank, tmp = tmp, rank
        if rank[sa[-1]] == n - 1:
            break
        k <<= 1
    return sa


def build_lcp(s, sa):
    """
    Kasai's algorithm to compute LCP array in O(n).
    lcp[i] = LCP(sa[i], sa[i+1]) for i in [0, n-2]
    """
    n = len(s)
    if n < 2:
        return []
    rank = [0] * n
    for i, p in enumerate(sa):
        rank[p] = i
    lcp = [0] * (n - 1)
    h = 0
    for i in range(n):
        r = rank[i]
        if r == 0:
            continue
        j = sa[r - 1]
        while i + h < n and j + h < n and s[i + h] == s[j + h]:
            h += 1
        lcp[r - 1] = h
        if h:
            h -= 1
    return lcp


def find_actual_lrs(seq):
    """Find the exact longest repeated substring using SA+LCP. Returns (substring, length)."""
    if not seq:
        return "", 0
    sa = build_suffix_array(seq)
    lcp = build_lcp(seq, sa)
    if not lcp:
        return "", 0
    max_len = max(lcp)
    if max_len == 0:
        return "", 0
    idx = lcp.index(max_len)
    start = sa[idx]
    return seq[start:start + max_len], max_len


def extract_maximal_repeats(seq, min_len, sa=None, lcp=None):
    """
    Extract repeated substrings of length >= min_len using SA+LCP.
    """
    n = len(seq)
    repeats = {}
    if n < 2 or min_len <= 0:
        return repeats
    if sa is None:
        sa = build_suffix_array(seq)
    if lcp is None:
        lcp = build_lcp(seq, sa)
    i = 0
    while i < len(lcp):
        if lcp[i] >= min_len:
            j = i
            cur_min = lcp[i]
            while j + 1 < len(lcp) and lcp[j + 1] >= min_len:
                j += 1
                cur_min = min(cur_min, lcp[j])
            occ = sorted(sa[k] for k in range(i, j + 2))
            substr = seq[sa[i]:sa[i] + cur_min]
            if substr not in repeats or len(repeats[substr]) < len(occ):
                repeats[substr] = occ
            i = j + 1
        else:
            i += 1
    return repeats


def extract_L_repeats_from_lcp(seq, sa, lcp, L, min_freq=2):
    """
    Extract all distinct substrings of length L that occur at least min_freq times
    """
    repeats = {}
    n = len(seq)
    if n < L or not sa or not lcp:
        return repeats
    i = 0
    while i < len(lcp):
        if lcp[i] >= L:
            j = i
            while j + 1 < len(lcp) and lcp[j + 1] >= L:
                j += 1
            idxs = [sa[k] for k in range(i, j + 2)]
            patterns = {}
            for p in idxs:
                if p + L <= len(seq):
                    pat = seq[p:p + L]
                    patterns.setdefault(pat, []).append(p)
            for pat, pos in patterns.items():
                if len(pos) >= min_freq:
                    if pat in repeats:
                        repeats[pat].extend(pos)
                    else:
                        repeats[pat] = list(pos)
            i = j + 1
        else:
            i += 1
    for pat in list(repeats.keys()):
        pos = sorted(set(repeats[pat]))
        if len(pos) < min_freq:
            del repeats[pat]
        else:
            repeats[pat] = pos
    return repeats


# ---------------------------
# Simple motif-transform helpers
# ---------------------------
def to_parameterized(s):
    mapping = {}
    next_id = 0
    out = []
    for ch in s:
        if ch not in mapping:
            mapping[ch] = str(next_id)
            next_id += 1
        out.append(mapping[ch])
    return "".join(out)


def to_abelian(s):
    return "".join(sorted(s))


def to_order_preserving(s):
    mapping = {}
    next_id = 0
    out = []
    for ch in s:
        if ch not in mapping:
            mapping[ch] = str(next_id)
            next_id += 1
        out.append(mapping[ch])
    return "".join(out)


def to_cartesian(s):
    if len(s) <= 1:
        return s
    out = []
    for a, b in zip(s, s[1:]):
        if a < b: out.append("<")
        elif a > b: out.append(">")
        else: out.append("=")
    return "".join(out)


def to_swap_robust(s):
    lst = list(s)
    for i in range(len(lst) - 1):
        if lst[i] > lst[i + 1]:
            lst[i], lst[i + 1] = lst[i + 1], lst[i]
    return "".join(lst)


# ---------------------------
# Core analysis
# ---------------------------
def analyze_species(file_path, name, l_values=L_VALUES, min_freq=MIN_FREQ):
    print(f"Analyzing {name} ({file_path}) ...")
    t0 = time.time()
    tracemalloc.start()

    seq = read_sequence(file_path)
    total_len = len(seq)

    sa = build_suffix_array(seq)
    lcp = build_lcp(seq, sa)

    exact_data = {}
    for L in l_values:
        if total_len < L:
            exact_data[L] = {}
            continue
        threshold = min_freq if isinstance(min_freq, int) else min_freq.get(L, MIN_FREQ)
        bk = extract_L_repeats_from_lcp(seq, sa, lcp, L, threshold)
        exact_data[L] = bk

    if sa and lcp:
        max_len = max(lcp) if lcp else 0
        if max_len > 0:
            idx = lcp.index(max_len)
            start = sa[idx]
            lrs_seq, lrs_len = seq[start:start + max_len], max_len
        else:
            lrs_seq, lrs_len = "", 0
    else:
        lrs_seq, lrs_len = "", 0

    maximal_at_min = extract_maximal_repeats(seq, MINIMAL_MAXIMAL_LEN, sa=sa, lcp=lcp)

    processed = {}
    for L in l_values:
        processed[L] = {
            'Exact': exact_data[L],
            'Parameterized': collections.defaultdict(list),
            'Cartesian': collections.defaultdict(list),
            'OrderPreserving': collections.defaultdict(list),
            'Abelian': collections.defaultdict(list),
            'Swaps': collections.defaultdict(list),
            'Maximal': {}
        }
        for pat, pos_list in exact_data[L].items():
            processed[L]['Parameterized'][to_parameterized(pat)].extend(pos_list)
            processed[L]['Cartesian'][to_cartesian(pat)].extend(pos_list)
            processed[L]['OrderPreserving'][to_order_preserving(pat)].extend(pos_list)
            processed[L]['Abelian'][to_abelian(pat)].extend(pos_list)
            processed[L]['Swaps'][to_swap_robust(pat)].extend(pos_list)

            if L >= MINIMAL_MAXIMAL_LEN and pat in maximal_at_min:
                processed[L]['Maximal'][pat] = maximal_at_min[pat]
            else:
                processed[L]['Maximal'][pat] = pos_list

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    runtime = time.time() - t0
    stats = {'runtime_sec': runtime, 'peak_mem_kb': peak / 1024.0, 'total_len': total_len}
    return processed, total_len, (lrs_seq, lrs_len), stats


def save_reports(name, data, total_len, lrs_info, stats=None):
    lrs_seq, lrs_len = lrs_info
    stats = stats or {}
    sum_path = os.path.join(OUTPUT_DIRS['summary'], f"{name}_Executive_Summary.txt")
    with open(sum_path, "w", encoding='utf-8') as f:
        f.write(f"EXECUTIVE REPORT: {name}\n" + "=" * 60 + "\n")
        f.write(f"DATA: Total Length: {total_len}\n")
        f.write(f"ACTUAL LRS: {lrs_seq}\nLength: {lrs_len}\n\n")
        if stats:
            f.write(f"ANALYSIS RUNTIME: {stats.get('runtime_sec', 0):.2f} s\n")
            f.write(f"PEAK MEMORY: {stats.get('peak_mem_kb', 0):.0f} KB\n\n")
        for L in sorted(data.keys()):
            exact = data[L]['Exact']
            f.write(f"--- L = {L} ---\n")
            if exact:
                counts = [len(v) for v in exact.values()]
                most = max(exact.items(), key=lambda x: len(x[1]))
                f.write(f"Total Unique Motifs: {len(exact)}\nAvg Frequency: {np.mean(counts):.2f}\n")
                f.write(f"Most Frequent Motif: {most[0]} ({len(most[1])} times)\n")
                topk = sorted(exact.items(), key=lambda x: len(x[1]), reverse=True)[:5]
                f.write("Top-5 motifs:\n")
                for pat, pos in topk: f.write(f"  {pat} | {len(pos)}\n")
            f.write("-" * 30 + "\n\n")

    for L in sorted(data.keys()):
        det_path = os.path.join(OUTPUT_DIRS['detailed'], f"{name}_L{L}_Detailed_Inventory.txt")
        with open(det_path, "w", encoding='utf-8') as f:
            f.write(f"DETAILED INVENTORY: {name} (L={L})\n" + "=" * 60 + "\n")
            for algo in ['Exact', 'Parameterized', 'Cartesian', 'OrderPreserving', 'Abelian', 'Swaps', 'Maximal']:
                f.write(f"\n[ALGORITHM: {algo}]\n")
                items = sorted(data[L][algo].items(), key=lambda x: len(x[1]), reverse=True)
                buf = []
                for pattern, pos in items:
                    buf.append(f"Pattern: {pattern} | Count: {len(pos)} | Pos: {pos}\n")
                    if len(buf) >= 1000:
                        f.write(''.join(buf))
                        buf = []
                f.write(''.join(buf))


def discover_data_files(data_dir=DATA_DIR):
    files = {}
    if not os.path.exists(data_dir): return files
    for root, _, filenames in os.walk(data_dir):
        for fn in filenames:
            if fn.lower().endswith(('.txt', '.fa', '.fasta')):
                files[os.path.splitext(fn)[0]] = os.path.join(root, fn)
    return files


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", "-d", default=DATA_DIR)
    parser.add_argument("--min-freq", type=int, default=MIN_FREQ)
    parser.add_argument("--L-values", nargs="+", type=int, default=L_VALUES)
    args = parser.parse_args(argv)

    files = discover_data_files(args.data_dir)
    master_data = {}
    for name, path in files.items():
        data, length, lrs_info, stats = analyze_species(path, name, l_values=args.L_values, min_freq=args.min_freq)
        master_data[name] = data
        save_reports(name, data, length, lrs_info, stats)

    if master_data:
        species = list(master_data.keys())
        all_shared_path = os.path.join(OUTPUT_DIRS['comparison'], 'Shared_By_ALL_Species.txt')
        with open(all_shared_path, "w", encoding='utf-8') as f_all:
            f_all.write("GENOMIC CORE ANALYSIS\n" + "="*55 + "\n")
            for L_cmp in args.L_values:
                all_sets = [set(master_data[s][L_cmp]['Exact'].keys()) for s in species]
                shared_by_all = set.intersection(*all_sets) if all_sets else set()
                f_all.write(f"\n[ AT LENGTH L = {L_cmp} ]\nShared by ALL: {len(shared_by_all)}\n")
                for i in range(len(species)):
                    for j in range(i + 1, len(species)):
                        s1, s2 = species[i], species[j]
                        num_shared = len(all_sets[i] & all_sets[j])
                        f_all.write(f"  - {s1} vs {s2}: {num_shared} motifs\n")
    print("All done.")

if __name__ == "__main__":
    main()