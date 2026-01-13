import matplotlib.pyplot as plt
import os
import re

INPUT_FILE = "output/comparison/Shared_By_ALL_Species.txt"
OUTPUT_DIR = "output/figures"

def parse_genomic_data(file_path):
    lengths = []
    shared_all = []
    chimp_human = []
    dog_human = []
    chimp_dog = []

    if not os.path.exists(file_path):
        print(f"Error: File {file_path} not found.")
        return None

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    blocks = content.split('[ AT LENGTH L = ')
    for block in blocks[1:]:
        # 1. Find the length (L) value
        l_val = int(block.split(' ]')[0])
        lengths.append(l_val)

        # 2. Find 'Shared by ALL' count
        all_match = re.search(r'shared by ALL 3 species: (\d+)', block)
        if all_match:
            shared_all.append(int(all_match.group(1)))

        # 3. Find 'chimpanzee vs human' count
        ch_match = re.search(r'chimpanzee vs human: (\d+)', block)
        if ch_match:
            chimp_human.append(int(ch_match.group(1)))
            
        # 4. Find 'dog vs human' count
        dh_match = re.search(r'dog vs human: (\d+)', block)
        if dh_match:
            dog_human.append(int(dh_match.group(1)))
        
        # 5. Find 'chimpanzee vs dog' count
        cd_match = re.search(r'chimpanzee vs dog: (\d+)', block)
        if cd_match:
            chimp_dog.append(int(cd_match.group(1)))

    return lengths, chimp_human, dog_human, shared_all, chimp_dog

def parse_summaries(summaries_dir):
    """Parse Executive Summary files and return mapping species -> (Ls, counts)"""
    if not os.path.exists(summaries_dir):
        return {}
    summary_files = [fn for fn in os.listdir(summaries_dir) if fn.endswith('_Executive_Summary.txt')]
    data = {}
    import re
    l_re = re.compile(r'^--- L = (\d+) ---')
    total_re = re.compile(r'^Total Unique Motifs:\s*(\d+)', re.IGNORECASE)
    for fn in summary_files:
        species = fn.replace('_Executive_Summary.txt', '')
        path = os.path.join(summaries_dir, fn)
        Ls = []
        totals = []
        with open(path, 'r', encoding='utf-8') as fh:
            for ln in fh:
                m = l_re.match(ln.strip())
                if m:
                    Ls.append(int(m.group(1)))
                    # next lines will contain Total Unique Motifs - scan ahead a few lines
                    for _ in range(10):
                        nxt = fh.readline()
                        if not nxt:
                            break
                        tm = total_re.match(nxt.strip())
                        if tm:
                            totals.append(int(tm.group(1)))
                            break
        if Ls and totals and len(Ls) == len(totals):
            # sort by L
            pairs = sorted(zip(Ls, totals), key=lambda x: x[0])
            Ls_s, totals_s = zip(*pairs)
            data[species] = (list(Ls_s), list(totals_s))
    return data


def create_visualizations():
    data = parse_genomic_data(INPUT_FILE)
    if not data: return
    
    lengths, ch_data, dh_data, all_data, cd_data = data

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    # Create the plot
    plt.figure(figsize=(12, 7))
    
    # Add the lines (use log scale because values range widely)
    plt.plot(lengths, ch_data, marker='o', label='Chimpanzee vs Human', linewidth=2, color='#1f77b4')
    plt.plot(lengths, dh_data, marker='^', label='Dog vs Human', linewidth=2, color='#ff7f0e')
    plt.plot(lengths, all_data, marker='s', label='Shared by ALL 3', linewidth=2, color='#d62728', linestyle='--')
    plt.plot(lengths, cd_data, marker='d', label='Chimpanzee vs Dog', linewidth=2, color='#9467bd')

    plt.yscale('log')
    plt.title('Genomic Similarity Across Species\n(Shared Motifs by Sequence Length)', fontsize=15)
    plt.xlabel('Motif Length (L)', fontsize=12)
    plt.ylabel('Number of Shared Motifs (Log Scale)', fontsize=12)
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(fontsize=11)

    # Add labels to the points for readability
    for i, txt in enumerate(ch_data):
        plt.annotate(f"{txt:,}", (lengths[i], ch_data[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9)

    output_file = os.path.join(OUTPUT_DIR, 'evolutionary_decay_complete.png')
    plt.savefig(output_file)
    plt.show()
    print(f"\n--- DONE! ---\nGraph saved to: {output_file}")

    # Additionally, create per-species motif-count graphs from Executive Summaries
    summaries_dir = os.path.join('output', 'summaries')
    species_data = parse_summaries(summaries_dir)
    if species_data:
        for species, (Ls, counts) in species_data.items():
            plt.figure(figsize=(10,6))
            plt.plot(Ls, counts, marker='o', linewidth=2)
            plt.yscale('log')
            plt.title(f'Motif Counts vs Length for {species.capitalize()}')
            plt.xlabel('Motif Length (L)')
            plt.ylabel('Total Unique Motifs (Log scale)')
            plt.grid(True, which='both', ls='-', alpha=0.2)
            for i, v in enumerate(counts):
                plt.annotate(f"{v:,}", (Ls[i], counts[i]), textcoords='offset points', xytext=(0,8), ha='center', fontsize=9)
            outf = os.path.join(OUTPUT_DIR, f'{species}_motif_counts_vs_L.png')
            plt.tight_layout()
            plt.savefig(outf)
            plt.close()
            print(f'Per-species graph saved: {outf}')
    else:
        print('No executive summaries found for per-species plots.')

if __name__ == "__main__":
    create_visualizations()