#!/usr/bin/env python3

import sys
import glob
import re
from collections import defaultdict


GFF_FILE = "/fs/cbcb-lab/mfritz13/FritzLab_RawData_Archive/North_etal_2023_reanalysis/data_files/other/Liftoff_Analysis/test1.gff"


def parse_outlier_file(filename):
    """
    Parse a single outlier file to extract outlier regions as intervals per chromosome.
    Returns a dict: chrom -> list of (start, end)
    """
    regions = defaultdict(list)
    chrom_pattern = re.compile(r'^(chrom[^\:]+):(\d+)-(\d+)$')

    with open(filename) as f:
        in_outlier_section = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Outlier positions"):
                in_outlier_section = True
                continue
            if not in_outlier_section:
                continue
            if line.startswith('chrom'):
                m = chrom_pattern.match(line)
                if m:
                    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
                    regions[chrom].append((start, end))
            else:
                continue

    # Merge overlapping intervals per chromosome
    for chrom in regions:
        regions[chrom] = merge_intervals(sorted(regions[chrom]))

    return regions


def merge_intervals(intervals):
    """
    Merge overlapping or adjacent intervals in a list.
    Intervals: list of (start, end)
    """
    if not intervals:
        return []

    merged = []
    current_start, current_end = intervals[0]

    for s, e in intervals[1:]:
        if s <= current_end + 1:  # overlapping or adjacent
            current_end = max(current_end, e)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = s, e
    merged.append((current_start, current_end))
    return merged


def coverage_array(intervals, length):
    """
    Create a coverage boolean array of given length marking positions covered by intervals.
    0-based index, length is max coordinate + 1 for safety.
    """
    cov = [False] * length
    for start, end in intervals:
        for pos in range(start, end + 1):
            if pos < length:
                cov[pos] = True
    return cov


def find_windows(chrom, cov_arrays, min_window_size=950, step_size=500, min_bp_fraction=0.8, min_sample_fraction=0.50):
    """
    Find conserved high depth windows on one chromosome.

    cov_arrays: list of boolean coverage arrays, one per sample, same length.
    For each window:
      - Check if >= min_bp_fraction of bases in window are covered in each sample (sample window coverage)
      - Check if >= min_sample_fraction of samples have this window coverage
    Merge adjacent/overlapping windows meeting the criteria.

    Returns list of tuples:
    (chrom, start, end, window_bp_fraction, sample_fraction, window_length)
    """
    n_samples = len(cov_arrays)
    length = len(cov_arrays[0])

    windows = []

    # Sliding windows, overlapping by half (step_size)
    window_starts = range(0, length, step_size)

    def calc_cov_frac(cov, s, e):
        length_w = e - s + 1
        if length_w <= 0:
            return 0
        covered = sum(cov[s:e + 1])
        return covered / length_w

    for start in window_starts:
        end = min(start + min_window_size - 1, length - 1)
        window_length = end - start + 1
        if window_length < min_window_size:
            continue

        sample_cov_flags = []
        bp_cov_sum = 0

        for cov in cov_arrays:
            frac = calc_cov_frac(cov, start, end)
            sample_cov_flags.append(frac >= min_bp_fraction)
            if frac >= min_bp_fraction:
                bp_cov_sum += frac

        n_samples_cov = sum(sample_cov_flags)
        sample_fraction = n_samples_cov / n_samples

        if n_samples_cov > 0:
            window_bp_fraction = bp_cov_sum / n_samples_cov
        else:
            window_bp_fraction = 0

        if sample_fraction >= min_sample_fraction and window_bp_fraction >= min_bp_fraction:
            windows.append([chrom, start, end, window_bp_fraction, sample_fraction, window_length])

    # Merge overlapping/adjacent windows that meet criteria
    merged = []
    if not windows:
        return merged

    current = windows[0]
    for w in windows[1:]:
        if w[1] <= current[2] + step_size:
            current[2] = max(current[2], w[2])
            total_len = current[5] + w[5]
            current[3] = (current[3] * current[5] + w[3] * w[5]) / total_len
            current[4] = (current[4] * current[5] + w[4] * w[5]) / total_len
            current[5] = total_len
        else:
            merged.append(current)
            current = w
    merged.append(current)
    return merged


def filter_chromosomes(chrom_list):
    """
    Keep chromosome names that start with 'chrom' only.
    """
    return [chrom for chrom in chrom_list if chrom.startswith('chrom')]


chrom_lengths = {
    f'chrom{i}': (15512169)
    for i in range(1, 32)
}

chrom_lengths['chrom1'] = 15512169
chrom_lengths['chrom2'] = 15061584
chrom_lengths['chrom3'] = 14695159
chrom_lengths['chrom4'] = 14165716
chrom_lengths['chrom5'] = 14051263
chrom_lengths['chrom6'] = 14023251
chrom_lengths['chrom7'] = 13649923
chrom_lengths['chrom8'] = 13649271
chrom_lengths['chrom9'] = 13615264
chrom_lengths['chrom10'] = 13258385
chrom_lengths['chrom11'] = 13088840
chrom_lengths['chrom12'] = 12972780
chrom_lengths['chrom13'] = 12841304
chrom_lengths['chrom14'] = 12757947
chrom_lengths['chrom15'] = 12605921
chrom_lengths['chrom16'] = 12332656
chrom_lengths['chrom17'] = 12254882
chrom_lengths['chrom18'] = 12245473
chrom_lengths['chrom19'] = 11629846
chrom_lengths['chrom20'] = 11599754
chrom_lengths['chrom21'] = 11461517
chrom_lengths['chrom22'] = 10871041
chrom_lengths['chrom23'] = 10422137
chrom_lengths['chrom24'] = 9742508
chrom_lengths['chrom25'] = 9280487
chrom_lengths['chrom26'] = 9278413
chrom_lengths['chrom27'] = 7355033
chrom_lengths['chrom28'] = 7211800
chrom_lengths['chrom29'] = 7203478
chrom_lengths['chrom30'] = 6316813
chrom_lengths['chrom31'] = 18805280


def overlaps_telomere(chrom, start, end):
    """
    Return True if [start, end] overlaps first or last 1% telomere regions,
    using known chromosome length from chrom_lengths dict.
    """
    chrom_length = chrom_lengths.get(chrom)
    if chrom_length is None:
        return False

    one_percent = int(chrom_length * 0.01)

    tel_start = (0, one_percent - 1)
    tel_end = (chrom_length - one_percent, chrom_length - 1)

    if end >= tel_start[0] and start <= tel_start[1]:
        return True
    if end >= tel_end[0] and start <= tel_end[1]:
        return True
    return False


def chrom_to_genbank(chrom):
    """
    Convert chrom# to GenBank accession format.
    E.g., chrom7 -> CM040209.1
    """
    chrom_num = int(chrom.replace('chrom', ''))
    return f"CM040{202 + chrom_num}.1"


def genbank_to_chrom(genbank):
    """
    Convert GenBank accession to chrom# format.
    E.g., CM040209.1 -> chrom7
    """
    match = re.match(r'CM040(\d+)\.1', genbank)
    if match:
        num = int(match.group(1)) - 202
        return f"chrom{num}"
    return None


def parse_gff_genes(gff_file):
    """
    Parse GFF file and extract gene features with their coordinates and names.
    Returns dict: chrom# -> list of (start, end, gene_name)
    """
    genes = defaultdict(list)
    
    with open(gff_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type != 'gene':
                continue
            
            genbank_chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            
            # Extract Note field
            note_match = re.search(r'Note=([^;]+)', attributes)
            if not note_match:
                continue
            
            gene_name = note_match.group(1)
            
            # Convert GenBank accession to chrom# format
            chrom = genbank_to_chrom(genbank_chrom)
            if chrom:
                genes[chrom].append((start, end, gene_name))
    
    return genes


def find_overlapping_genes(chrom, window_start, window_end, genes_dict):
    """
    Find all genes that overlap with the given window.
    Returns list of gene names.
    """
    overlapping_genes = []
    
    if chrom not in genes_dict:
        return overlapping_genes
    
    for gene_start, gene_end, gene_name in genes_dict[chrom]:
        # Check for overlap: window and gene overlap if they don't NOT overlap
        # NOT overlap means: window_end < gene_start OR window_start > gene_end
        # So overlap means: NOT (window_end < gene_start OR window_start > gene_end)
        # Which simplifies to: window_end >= gene_start AND window_start <= gene_end
        if window_end >= gene_start and window_start <= gene_end:
            overlapping_genes.append(gene_name)
    
    return overlapping_genes


def main(file_pattern, output_file):
    print(f"Reading input files matching: {file_pattern}")
    files = sorted(glob.glob(file_pattern))
    if not files:
        print("No input files found. Exiting.")
        sys.exit(1)

    print(f"Found {len(files)} files.")

    # Parse GFF file
    print(f"Parsing GFF file: {GFF_FILE}")
    genes_dict = parse_gff_genes(GFF_FILE)
    print(f"Found genes on {len(genes_dict)} chromosomes")

    sample_regions = []
    all_chromosomes = set()

    for f in files:
        regions = parse_outlier_file(f)
        regions = {chrom: intervals for chrom, intervals in regions.items() if chrom.startswith('chrom')}
        sample_regions.append(regions)
        all_chromosomes.update(regions.keys())

    def chromosome_sort_key(chrom):
        return int(chrom.replace("chrom", ""))

    all_chromosomes = filter_chromosomes(sorted(all_chromosomes, key=chromosome_sort_key))
    print(f"Chromosomes to analyze: {all_chromosomes}")

    results = []

    for chrom in all_chromosomes:
        max_pos = 0
        for regions in sample_regions:
            if chrom in regions:
                max_pos = max(max_pos, max(end for _, end in regions[chrom]))
        max_pos += 1

        cov_arrays = []
        for regions in sample_regions:
            intervals = regions.get(chrom, [])
            cov = coverage_array(intervals, max_pos)
            cov_arrays.append(cov)

        merged_windows = find_windows(chrom, cov_arrays)

        # Filter windows: skip telomeres and windows without genes
        for w in merged_windows:
            chrom_, start, end = w[0], w[1], w[2]
            
            if overlaps_telomere(chrom_, start, end):
                continue
            
            # Find overlapping genes
            overlapping_genes = find_overlapping_genes(chrom_, start, end, genes_dict)
            
            if not overlapping_genes:
                continue
            
            # Add gene names to result
            gene_names = ','.join(overlapping_genes)
            results.append(w + [gene_names])

    with open(output_file, 'w') as out:
        header = ['chromosome', 'start', 'end', 'percent_bp_covered', 'percent_samples', 'window_length', 'gene_name']
        out.write('\t'.join(header) + '\n')
        for r in results:
            out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]*100:.1f}\t{r[4]*100:.1f}\t{r[5]}\t{r[6]}\n")

    print(f"Done. Results written to {output_file}")
    print(f"Total windows with genes: {len(results)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} '<input_files_pattern>' <output_file>")
        print("Example:")
        print(f"  {sys.argv[0]} 'trimmed_*_outliers12x.txt' conserved_windows.txt")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])