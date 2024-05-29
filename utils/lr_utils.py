import sys

import pandas as pd
import pysam
from genomicranges import GenomicRanges
import cigar

def warn_diff_read_ids(gr1_ids, gr2_ids):
    if len(gr1_ids) != len(gr2_ids):
        sys.stderr.write('Warning: Reads are not the same between the two haplotypes\n')
        #     list the read IDs that are not in both haplotypes
        for read_id in gr1_ids.difference(gr2_ids):
            sys.stderr.write('Read ID: ' + read_id + ' not in haplotype 2\n')
        for read_id in gr2_ids.difference(gr1_ids):
            sys.stderr.write('Read ID: ' + read_id + ' not in haplotype 1\n')

def extract_reads(alignment_file):
    a = pysam.AlignmentFile(alignment_file, "rb")
    reads = []
    for read in a:
        reads.append(read)
    return reads


def create_pandas_df(reads):
    d = []
    for read in reads:
        seq_name = read.reference_name
        # if the sequence name is a non-standard chromosome, add chr0 to the sequence name There appears to be a bug
        # in the reduce function that causes it to fail when the sequence name is not a standard chromosome
        if not seq_name.startswith('chr'):
            # replace "_" with "." in the sequence name
            seq_name = seq_name.replace("_", ".")
            seq_name = seq_name + '_0'
        d.append({'seqnames': seq_name, 'strand': read.is_reverse, 'cigar': read.cigarstring,
                  'read_name': read.query_name, 'starts': read.reference_start, 'ends': read.reference_end,
                  'sequence': read.query_sequence})

    return pd.DataFrame(d)


def get_chr_start_stop(gr_row):
    return gr_row.get_seqnames()[0], gr_row.get_start()[0], gr_row.get_end()[0]


# checks that the reads represent a continous region of the reference genome
def convert_to_range(reads):
    grr = GenomicRanges.from_pandas(create_pandas_df(reads))
    if len(set(grr.get_seqnames())) > 1:
        sys.stderr.write('Error: Reads are not from the same chromosome\n')
        sys.exit(1)
    sys.stderr.write('Reads are from the same chromosome\n')
    reduce = grr.reduce(ignore_strand=True)
    if len(reduce) > 1:
        sys.stderr.write('Error: Reads span  discontinuous range\n')
        sys.exit(1)
    return grr, reduce

def add_cigar_span(gr):
    cigar_span = []
    for i in range(len(gr)):
        cigar_ops = cigar.Cigar(gr.mcols.get_column('cigar')[i])
        span = get_cigar_span_of_read(cigar_ops)
        # if gr.get_strand()[i] == 1:
        #    span = reverse_span(span, gr.get_width()[i])
        cigar_span.append(span)
    gr.mcols.set_column('cigar_span', cigar_span, in_place=True)
    return gr

def consolidate_cigar_spans_to_read_id(gr):
    #     get the cigar spans for each read
    gr = add_cigar_span(gr)
    #     get the cigar spans for each read
    read_name_to_cigar_span = {}
    read_name_ref_span = {}
    for i in range(len(gr)):
        read_name = gr.mcols.get_column('read_name')[i]
        if read_name not in read_name_to_cigar_span:
            read_name_to_cigar_span[read_name] = []
        if read_name not in read_name_ref_span:
            read_name_ref_span[read_name] = []
        read_name_to_cigar_span[read_name].append(gr.mcols.get_column('cigar_span')[i])
        read_name_ref_span[read_name].append(get_chr_start_stop(gr[i]))
    return read_name_to_cigar_span, read_name_ref_span



# find the indices of the first and last non-clipped bases in the read
def get_cigar_span_of_read(cigar_s):
    start = 0
    end = 0
    finding_start = True
    for c in cigar_s.items():
        if c[1] in ['M', 'I', 'X', '=']:
            finding_start = False
            end += c[0]
        elif c[1] in ['S', 'H'] and finding_start:
            start += c[0]
    return start, end + start