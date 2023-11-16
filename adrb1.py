import pysam
import os

import functions

HAPLOTYPES = ('ag', 'ac', 'gc', 'gg')
MAPQ_MIN_THRESHOLD = 60

if __name__ == '__main__':

    sample_names = [fn for fn in os.listdir('bam') if os.path.isdir('bam/'+fn)]
    assert sample_names, 'No samples found in bam directory'
    rows = []
    for sample_name in sample_names:
        bam_path = functions.get_bam_path(sample_name)
        if bam_path is None: continue
        bamfile = pysam.AlignmentFile(bam_path, 'rb')
        passing_reads, num_reads_w_mapq_below_threshold, times_target_locus_not_found = functions.get_qc_metrics(bamfile)
        reads_by_haplotype = functions.get_haplotyped_read_names(passing_reads)
        functions.write_reads_by_haplotype(sample_name, reads_by_haplotype)
        total_mapped_reads = functions.count_mapped_reads(bamfile)
        reads_passing_qc = total_mapped_reads - num_reads_w_mapq_below_threshold - times_target_locus_not_found - len(reads_by_haplotype['invalid'])
        ht_counts = tuple([len(reads_by_haplotype[ht]) for ht in HAPLOTYPES])
        rows.append((sample_name, total_mapped_reads, reads_passing_qc, len(reads_by_haplotype['invalid'])) + ht_counts)

    # write the haplotype counts to a file
    header = ('sample_name', 'total_mapped_reads', 'reads_passing_qc', 'invalid_haplotype_reads') + tuple('haplotype_' + ht + '_reads' for ht in HAPLOTYPES)
    with open('haplotype_counts.tsv', 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rows:
            f.write('\t'.join(str(x) for x in row) + '\n')