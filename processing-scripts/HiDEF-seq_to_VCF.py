#!/usr/bin/env python

################################################
#
#   Convert HIDEF-seq TSV output to VCF
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################


################################################
#   Libraries
################################################
import sys, argparse, os


################################################
#   Global
################################################
# VCF VERSION
VERSION = '##fileformat=VCFv4.2'

# INFO DEFINITIONS
STRANDTYPE = '##INFO=<ID=STRANDTYPE,Number=1,Type=String,Description="ssDNA (ssDNA mismatch) or dsDNA (dsDNA mutation)">'
SYNTHESIZEDSTRAND = '##INFO=<ID=SYNTHESIZEDSTRAND,Number=1,Type=String,Description="The strand of the reference genome to which the strand synthesized during sequencing aligned (+ or -)">'
TEMPLATESTRAND = '##INFO=<ID=TEMPLATESTRAND,Number=1,Type=String,Description="The strand of the reference genome to which the \'template strand being copied during sequencing\' aligned (+ or -). This is the reverse of SYNTHESIZEDSTRAND">'
REF_REFPLUSSTRAND = '##INFO=<ID=REF_REFPLUSSTRAND,Number=1,Type=String,Description="Sequence of the reference genome on the reference plus strand">'
ALT_REFPLUSSTRAND = '##INFO=<ID=ALT_REFPLUSSTRAND,Number=1,Type=String,Description="Sequence of the call on the reference plus strand">'
TNC_REFPLUSSTRAND = '##INFO=<ID=TNC_REFPLUSSTRAND,Number=1,Type=String,Description="Trinucleotide context of the reference genome on the reference plus strand">'
REF32_REFPLUSSTRAND = '##INFO=<ID=REF32_REFPLUSSTRAND,Number=1,Type=String,Description="Sequence of the reference genome on the reference plus strand, then transformed to the strand where the reference base is a pyrimidine">'
ALT32_REFPLUSSTRAND = '##INFO=<ID=ALT32_REFPLUSSTRAND,Number=1,Type=String,Description="Sequence of the call on the reference plus strand, then transformed to the strand where the reference base is a pyrimidine">'
TNC32_REFPLUSSTRAND = '##INFO=<ID=TNC32_REFPLUSSTRAND,Number=1,Type=String,Description="Trinucleotide context of the reference genome on the reference plus strand, then transformed to the strand where the reference base is a pyrimidine">'
REF_SSDNASYNTHESIZEDSTRAND = '##INFO=<ID=REF_SSDNASYNTHESIZEDSTRAND,Number=1,Type=String,Description="Sequence of the reference genome on the strand synthesized during sequencing">'
ALT_SSDNASYNTHESIZEDSTRAND = '##INFO=<ID=ALT_SSDNASYNTHESIZEDSTRAND,Number=1,Type=String,Description="Sequence of the call on the strand synthesized during sequencing">'
TNC_SSDNASYNTHESIZEDSTRAND = '##INFO=<ID=TNC_SSDNASYNTHESIZEDSTRAND,Number=1,Type=String,Description="Trinucleotide context of the reference genome on the strand synthesized during sequencing">'
REF_SSDNATEMPLATESTRAND = '##INFO=<ID=REF_SSDNATEMPLATESTRAND,Number=1,Type=String,Description="Sequence of the reference genome on the template strand being copied during sequencing">'
ALT_SSDNATEMPLATESTRAND = '##INFO=<ID=ALT_SSDNATEMPLATESTRAND,Number=1,Type=String,Description="Sequence of the call on the template strand being copied during sequencing">'
TNC_SSDNATEMPLATESTRAND = '##INFO=<ID=TNC_SSDNATEMPLATESTRAND,Number=1,Type=String,Description="Trinucleotide context of the reference genome on the template strand being copied during sequencing">'

# FORMAT DEFINITIONS
RQ = '##FORMAT=<ID=RQ,Number=2,Type=Float,Description="Read quality (\'rq\' tag) of the consensus sequence [FORWARD_STRAND, REVERSE_STRAND]">'
ECOV = '##FORMAT=<ID=ECOV,Number=2,Type=Float,Description="Effective coverage (\'ec\' tag; i.e., average number of subreads across consensus sequence) of the consensus sequence [FORWARD_STRAND, REVERSE_STRAND]">'
MAPQ = '##FORMAT=<ID=MAPQ,Number=2,Type=Integer,Description="Mapping quality of the consensus sequence [FORWARD_STRAND, REVERSE_STRAND]">'
QQ = '##FORMAT=<ID=QQ,Number=2,Type=Integer,Description="Base quality of the call position on the consensus sequence [FORWARD_STRAND, REVERSE_STRAND]">'
BPD = '##FORMAT=<ID=BPD,Number=1,Type=Integer,Description="Distance of the call in bases from the start or end of the consensus sequence alignment, whichever is closer. Calculated for dsDNA calls on each end using the strand that is closer to the call">'
MPS = '##FORMAT=<ID=MPS,Number=1,Type=Integer,Description="Number of post-filtering ssDNA or dsDNA calls called on the strand (for ssDNA calls) or duplex DNA molecule (for dsDNA calls), respectively, on which this call was made">'
VR = '##FORMAT=<ID=VR,Number=3,Type=Integer,Description="Number of reads with the same sequence as the call, at the position of the call [GERMLINE_READS, FORWARD_STRAND_SUBREADS, REVERSE_STRAND_SUBREADS]">'
TOT = '##FORMAT=<ID=TOT,Number=1,Type=Integer,Description="Total number of germline reads at the position of the call">'
VAF = '##FORMAT=<ID=VAF,Number=3,Type=Float,Description="Variant allele frequency of reads with the same sequence as the call, at the position of the call [GERMLINE_READS, FORWARD_STRAND_SUBREADS, REVERSE_STRAND_SUBREADS]">'
ALN = '##FORMAT=<ID=ALN,Number=2,Type=Float,Description="Fraction of subreads aligning to the position of the call [FORWARD_STRAND_SUBREADS, REVERSE_STRAND_SUBREADS]">'


################################################
#   HIDEFSeqVariant
################################################
class HIDEFSeqVariant(object):

    def __init__(self, tsv_line):
        '''
        '''
        self.sampleid, \
        self.genome, \
        self.filterid, \
        self.zmw_name, \
        self.zmw_hole, \
        self.isizefwd, \
        self.isizerev, \
        self.zmw_isize, \
        self.chrom, \
        self.pos, \
        self.strandtype, \
        self.synthesizedstrand, \
        self.templatestrand, \
        self.ref_refplusstrand, \
        self.alt_refplusstrand, \
        self.tnc_refplusstrand, \
        self.ref32_refplusstrand, \
        self.alt32_refplusstrand, \
        self.tnc32_refplusstrand, \
        self.ref_ssDNAsynthesizedstrand, \
        self.alt_ssDNAsynthesizedstrand, \
        self.tnc_ssDNAsynthesizedstrand, \
        self.ref_ssDNAtemplatestrand, \
        self.alt_ssDNAtemplatestrand, \
        self.tnc_ssDNAtemplatestrand, \
        self.zmw_rqfwd, \
        self.zmw_rqrev, \
        self.zmw_ecfwd, \
        self.zmw_ecrev, \
        self.zmw_mapqfwd, \
        self.zmw_mapqrev, \
        self.qqfwd, \
        self.qqrev, \
        self.bpenddist, \
        self.mutationsperstrandorzmw, \
        self.GermlineVariantReads, \
        self.GermlineTotalReads, \
        self.GermlineVAF, \
        self.fwdsubreadsVariantReads, \
        self.fwdsubreadsVAF, \
        self.fwdsubreadscvgfraction, \
        self.revsubreadsVariantReads, \
        self.revsubreadsVAF, \
        self.revsubreadscvgfraction = tsv_line.rstrip().split()

    def format_INFO(self):
        ''' Generates INFO column values for variant
        '''
        INFO = []
        tags = [
            'strandtype', 'synthesizedstrand', 'templatestrand',
            'ref_refplusstrand', 'alt_refplusstrand', 'tnc_refplusstrand',
            'ref32_refplusstrand', 'alt32_refplusstrand', 'tnc32_refplusstrand',
            'ref_ssDNAsynthesizedstrand', 'alt_ssDNAsynthesizedstrand', 'tnc_ssDNAsynthesizedstrand',
            'ref_ssDNAtemplatestrand', 'alt_ssDNAtemplatestrand', 'tnc_ssDNAtemplatestrand'
            ]

        for tag in tags:
            val = getattr(self, tag)
            if val != 'NA':
                INFO.append(f'{tag.upper()}={val}')

        return ';'.join(INFO)

    def format_FORMAT_GENOTYPE(self):
        ''' Generates FORMAT and GENOTYPE columns values for variant
        '''
        FORMAT = []
        GENOTYPE = []
        ids_tags = {
            'RQ': ['zmw_rqfwd', 'zmw_rqrev'],
            'ECOV': ['zmw_ecfwd', 'zmw_ecrev'],
            'MAPQ': ['zmw_mapqfwd', 'zmw_mapqrev'],
            'QQ': ['qqfwd', 'qqrev'],
            'BPD': ['bpenddist'],
            'MPS': ['mutationsperstrandorzmw'],
            'VR': ['GermlineVariantReads', 'fwdsubreadsVariantReads', 'revsubreadsVariantReads'],
            'TOT': ['GermlineTotalReads'],
            'VAF': ['GermlineVAF', 'fwdsubreadsVAF', 'revsubreadsVAF'],
            'ALN': ['fwdsubreadscvgfraction', 'revsubreadscvgfraction']
            }

        for id, tags in ids_tags.items():
            FORMAT.append(id)
            vals = []
            for tag in tags:
                val = getattr(self, tag)
                if val != 'NA':
                    vals.append(val)
                else:
                    vals.append('.')
            GENOTYPE.append(vals)

        return ':'.join(FORMAT), ':'.join(','.join(vals) for vals in GENOTYPE)

    def to_VCF(self):
        ''' Converts variant to VCF
        '''
        m_ = '.' # placeholder for missing value

        # REF & ALT
        REF, ALT = m_, m_
        if self.strandtype == 'dsDNA':
            REF = self.ref_refplusstrand
            ALT = self.alt_refplusstrand
        else:
            REF = self.ref_ssDNAtemplatestrand
            ALT = self.alt_ssDNAtemplatestrand
        # FORMAT & GENOTYPE
        FORMAT, GENOTYPE = self.format_FORMAT_GENOTYPE()

        return f'{self.chrom}\t{self.pos}\t{m_}\t{REF}\t{ALT}\t{m_}\t{m_}\t{self.format_INFO()}\t{FORMAT}\t{GENOTYPE}\n'


################################################
#   Functions
################################################
def format_header(VCF_version, INFO_tags, FORMAT_tags, sampleid, genome):
    '''
    '''
    reference = f'##reference={genome}'
    columns = f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sampleid}\n'
    header = [VCF_version, reference] + INFO_tags + FORMAT_tags
    return '\n'.join(header) + '\n' + columns

def main(args):
    '''
    '''
    # Args
    inputfile = args['inputfile']
    basename = args['outprefix'] if args['outprefix'] else os.path.splitext(os.path.basename(inputfile))[0]

    output_buffers = {}

    with open(inputfile) as fi:
        for i, line in enumerate(fi):
            if i != 0: # Removing the header
                var = HIDEFSeqVariant(line)
                samplefilter = f'{var.sampleid}_{var.filterid}'
                if samplefilter in output_buffers:
                    output_buffers[samplefilter].write(var.to_VCF())
                else:
                    filename_ = f'{basename}_{samplefilter}.vcf'
                    buffer_ = open(filename_, 'w')
                    output_buffers.setdefault(samplefilter, buffer_)
                    header = format_header(
                            VERSION,
                            [
                                STRANDTYPE, SYNTHESIZEDSTRAND, TEMPLATESTRAND,
                                REF_REFPLUSSTRAND, ALT_REFPLUSSTRAND, TNC_REFPLUSSTRAND,
                                REF32_REFPLUSSTRAND, ALT32_REFPLUSSTRAND, TNC32_REFPLUSSTRAND,
                                REF_SSDNASYNTHESIZEDSTRAND, ALT_SSDNASYNTHESIZEDSTRAND, TNC_SSDNASYNTHESIZEDSTRAND,
                                REF_SSDNATEMPLATESTRAND, ALT_SSDNATEMPLATESTRAND, TNC_SSDNATEMPLATESTRAND
                            ],
                            [RQ, ECOV, MAPQ, QQ, BPD, MPS, VR, TOT, VAF, ALN],
                            var.sampleid, var.genome
                            )
                    buffer_.write(header)
                    buffer_.write(var.to_VCF())

    for _, buffer_ in output_buffers.items():
        buffer_.close()


################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert HIDEF-seq TSV output to VCF')

    parser.add_argument('-i','--inputfile', help='HIDEF-seq output TSV file', required=True)
    parser.add_argument('--outprefix', help='prefix to use for the outputfile', required=False)

    args = vars(parser.parse_args())

    main(args)
