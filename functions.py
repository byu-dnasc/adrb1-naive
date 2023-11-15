import os

from adrb1 import HAPLOTYPES, MAPQ_MIN_THRESHOLD

def count_mapped_reads(bamfile):
    return sum(1 for _ in bamfile.fetch())

def get_aligned_pos(pairs, pos):
    '''
    Get the read position corresponding to a reference position
    Returns none if the reference position is not found or if
    the read position does not exist
    '''
    for read_pos, ref_pos in pairs:
        if ref_pos == pos:
            return read_pos # read_pos is None if the read position does not exist
    return None # ref_pos not found

def get_qc_metrics(bamfile) -> tuple:
    passing_reads = []
    reads_w_mapq_below_threshold = 0
    target_locus_not_found = 0
    for alignment in bamfile.fetch():
        passed = True
        if alignment.mapping_quality < MAPQ_MIN_THRESHOLD:
            passed = False
            reads_w_mapq_below_threshold += 1
        pairs = alignment.get_aligned_pairs()
        aligned_pos_1 = get_aligned_pos(pairs, 144)
        aligned_pos_2 = get_aligned_pos(pairs, 1164)
        if aligned_pos_1 is None or aligned_pos_2 is None:
            passed = False
            target_locus_not_found += 1
        if passed:
            passing_reads.append(alignment)
    return passing_reads, reads_w_mapq_below_threshold, target_locus_not_found

def get_haplotyped_read_names(passing_reads) -> tuple:
    reads_by_haplotype = {ht: [] for ht in HAPLOTYPES + ('invalid',)}

    for read in passing_reads:
        pairs = read.get_aligned_pairs()
        aligned_pos_1 = get_aligned_pos(pairs, 144)
        aligned_pos_2 = get_aligned_pos(pairs, 1164)
        snp_1 = read.query_sequence[aligned_pos_1]
        snp_2 = read.query_sequence[aligned_pos_2]
        haplotype = (snp_1+snp_2).lower()
        if haplotype in HAPLOTYPES:
            reads_by_haplotype[haplotype].append(read.query_name)
        else:
            reads_by_haplotype['invalid'].append(read.query_name)

    return reads_by_haplotype

def write_reads_by_haplotype(sample_name, reads_by_haplotype, failed_haplotyping):
    '''Write the read names for each haplotype to a file'''
    os.makedirs(f'reads_by_haplotype/{sample_name}', exist_ok=True)
    for ht in HAPLOTYPES:
        with open(f'reads_by_haplotype/{sample_name}/haplotype_{ht}_reads.txt', 'w') as f:
            for read_name in reads_by_haplotype[ht]:
                f.write(read_name + '\n')
    with open(f'reads_by_haplotype/{sample_name}/haplotype_unknown_reads.txt', 'w') as f:
        for read_name in failed_haplotyping:
            f.write(read_name + '\n')

