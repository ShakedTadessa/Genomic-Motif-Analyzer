#!/usr/bin/env python3
"""
Compare Detailed Inventory files and write:
- motifs shared by ALL species per L
- motifs shared by each PAIR of species per L

Each result is written to its own file.
NO summary file is created.
"""

import os
import re
import argparse
from collections import defaultdict
from pathlib import Path

# Try to reuse OUTPUT_DIRS from finalproject.py if available
try:
    from finalproject import OUTPUT_DIRS
    COMP_DIR = OUTPUT_DIRS.get('comparison', os.path.join('output', 'comparison'))
    DETAILED_DIR = OUTPUT_DIRS.get('detailed', os.path.join('output', 'detailed'))
except Exception:
    COMP_DIR = os.path.join('output', 'comparison')
    DETAILED_DIR = os.path.join('output', 'detailed')

PAT_RE = re.compile(r'^Pattern:\s*(.*?)\s*\|\s*Count:', re.IGNORECASE)
FNAME_RE = re.compile(r'(.+)_L(\d+)_Detailed_Inventory', re.IGNORECASE)


def parse_exact_patterns(path: Path):
    """Extract Exact motifs from a Detailed Inventory file."""
    patterns = set()
    in_exact = False

    with path.open('r', encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')

            if line.startswith("[ALGORITHM:"):
                in_exact = line.upper().startswith("[ALGORITHM: EXACT")
                continue

            if not in_exact:
                continue

            m = PAT_RE.match(line)
            if m:
                patterns.add(m.group(1).strip())

    return patterns


def collect_by_L(detailed_dir: str):
    """Build mapping: data[L][species] = set of motifs"""
    data = defaultdict(lambda: defaultdict(set))

    for fn in os.listdir(detailed_dir):
        if not fn.endswith('_Detailed_Inventory.txt'):
            continue

        m = FNAME_RE.match(fn)
        if not m:
            continue

        species, L = m.group(1), int(m.group(2))
        path = Path(detailed_dir) / fn

        motifs = parse_exact_patterns(path)
        data[L][species].update(motifs)

    return data


def write_shared_outputs(data, out_dir):
    """Write shared motif files (ALL species + pairwise)."""
    os.makedirs(out_dir, exist_ok=True)

    for L in sorted(data.keys()):
        species_list = sorted(data[L].keys())

        # ===== Shared by ALL species =====
        sets = [data[L][s] for s in species_list]
        shared_all = set.intersection(*sets) if sets else set()

        all_out = Path(out_dir) / f"shared_all_L{L}.txt"
        all_out.write_text("\n".join(sorted(shared_all)), encoding="utf-8")

        # ===== Pairwise shared =====
        for i in range(len(species_list)):
            for j in range(i + 1, len(species_list)):
                a, b = species_list[i], species_list[j]
                shared_pair = data[L][a] & data[L][b]

                pair_out = Path(out_dir) / f"shared_L{L}_{a}_vs_{b}.txt"
                pair_out.write_text("\n".join(sorted(shared_pair)), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Compare motif inventories across species")
    parser.add_argument("--detailed-dir", default=DETAILED_DIR)
    parser.add_argument("--out-dir", default=COMP_DIR)
    args = parser.parse_args()

    data = collect_by_L(args.detailed_dir)

    if not data:
        print(f"No detailed inventory files found in: {args.detailed_dir}")
        return

    write_shared_outputs(data, args.out_dir)
    print(f"Done. Shared motif files written to: {args.out_dir}")


if __name__ == "__main__":
    main()

