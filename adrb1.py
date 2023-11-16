import pysam
import os

from functions import *

HAPLOTYPES = ('ag', 'ac', 'gc', 'gg')
MAPQ_MIN_THRESHOLD = 60
INVALID_HAPLOTYPE = 'NA'

if __name__ == '__main__':

    sample_names = [fn for fn in os.listdir('bam') if os.path.isdir(f'bam/{fn}')]
    assert sample_names, 'No samples found in bam directory'
    rows = []
    for sample_name in sample_names:
        bam_path = get_bam_path(sample_name)
        if bam_path is None: continue
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
        # note: here, "read" refers to a sequence which may represent a single- or double-stranded molecule
        # "strand" refers to a single-stranded molecule
        passing_reads, num_reads_w_mapq_below_threshold, times_target_locus_not_found = get_qc_metrics(bamfile)
        reads_by_haplotype = get_haplotyped_read_names(passing_reads)
        write_reads_by_haplotype(sample_name, reads_by_haplotype)
        total_mapped_strands = count_mapped_strands(bamfile)
        strands_passing_qc = sum(map(lambda alignment: num_strands(alignment.query_name), passing_reads))
        ht_counts = [sum(map(num_strands, reads_by_haplotype[ht])) for ht in (INVALID_HAPLOTYPE,) + HAPLOTYPES]
        rows.append((sample_name, total_mapped_strands, strands_passing_qc, num_reads_w_mapq_below_threshold, times_target_locus_not_found) + tuple(ht_counts))

    # write the haplotype counts to a file
    header = ('sample_name', 'total_mapped_strands', 'passing_qc', 'mapq_below_threshold', 'missing_target_locus', 'haplotype_invalid') + tuple(f'haplotype_{ht}' for ht in HAPLOTYPES)
    with open('haplotype_counts.tsv', 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')