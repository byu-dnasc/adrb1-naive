import os
import glob

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

def write_reads_by_haplotype(sample_name, reads_by_haplotype):
    '''Write the read names for each haplotype to a file'''
    path = f'reads_by_haplotype/{sample_name}'
    os.makedirs(path, exist_ok=True)
    with open(f'{path}/reads_by_haplotype.csv', 'w') as f:
        for ht in reads_by_haplotype:
            for read_name in reads_by_haplotype[ht]:
                f.write(f'{ht},{read_name}\n')

def get_bam_path(sample_name):
    '''Get bam file or handle missing or too many bam files'''
    sample_bam_dir = f'bam/{sample_name}/'
    bam_files = glob.glob(f'{sample_bam_dir}/*.bam')
    if not bam_files:
        print(f'No bam files found in {sample_bam_dir}')
        return None
    else:
        assert len(bam_files) == 1, f'Found more than one bam file in {sample_bam_dir}'
    return bam_files[0]

