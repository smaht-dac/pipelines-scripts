#!/usr/bin/env python

########################################################################
#
#   Author: Michele Berselli
#       Harvard Medical School
#       berselli.michele@gmail.com
#
#   Script to parse the output from several QC software to generate
#       portal-compatible Quality Metric (QM) objects
#
########################################################################

########################################################################
# Libraries
########################################################################
import sys, os
import json
import zipfile

########################################################################
# Global Variables
########################################################################
BAMQM_NAME='BAM Quality Metrics'
VALUES_FILENAME='qc_values.json'
ARCHIVE_NAME='metrics.zip'
# Samtools
SAMTOOLS_stats_OUTPUT='output.samtools.stats'
SAMTOOLS_flagstat_OUTPUT='output.samtools.flagstat'
SAMTOOLS_idxstats_OUTPUT='output.samtools.idxstats'
# Picard
PICARD_CollectAlignmentSummaryMetrics_OUTPUT='output.picard.CollectAlignmentSummaryMetrics'
PICARD_CollectBaseDistributionByCycle_OUTPUT='output.picard.CollectBaseDistributionByCycle'
PICARD_CollectBaseDistributionByCycle_PDF='collect_base_dist_by_cycle.pdf'
PICARD_CollectGcBiasMetrics_OUTPUT='output.picard.CollectGcBiasMetrics'
PICARD_CollectGcBiasMetrics_SUMMARY='summary.picard.CollectGcBiasMetrics'
PICARD_CollectGcBiasMetrics_PDF='gc_bias_metrics.pdf'
PICARD_CollectInsertSizeMetrics_OUTPUT='output.picard.CollectInsertSizeMetrics'
PICARD_CollectInsertSizeMetrics_PDF='insert_size_histogram.pdf'
PICARD_CollectWgsMetrics_OUTPUT='output.picard.CollectWgsMetrics'
PICARD_MeanQualityByCycle_OUTPUT='output.picard.MeanQualityByCycle'
PICARD_MeanQualityByCycle_PDF='mean_qual_by_cycle.pdf'

########################################################################
# Tooltips
########################################################################
samtools_stats = {
    'raw total sequences': ('Total Sequences [Samtools]', 'Total number of reads in a file, excluding supplementary and secondary reads'),
    'sequences': ('Processed Sequences [Samtools]', 'Number of processed reads'),
    '1st fragments': ('1st Fragments [Samtools]', 'Number of first fragment reads'),
    'last fragments': ('Last Fragments [Samtools]', 'Number of last fragment reads'),
    'reads paired': ('Reads Paired [Samtools]', 'Number of paired reads, mapped or unmapped, that are neither secondary nor supplementary (paired-end technology bit set)'),
    'reads mapped': ('Reads Mapped [Samtools]', 'Number of reads, paired or single, that are mapped'),
    'reads mapped and paired': ('Reads Mapped and Paired [Samtools]', 'Number of mapped paired reads (paired-end technology bit set with both mates mapped)'),
    'reads properly paired': ('Reads Properly Paired [Samtools]', 'Number of mapped paired reads (proper-pair bit set, both mates mapped within the expected distance)'),
    'reads unmapped': ('Reads Unmapped [Samtools]', 'Number of unmapped reads'),
    'reads duplicated': ('Reads Duplicated [Samtools]', 'Number of duplicate reads'),
    'reads MQ0': ('Reads MQ0 [Samtools]', 'Number of mapped reads with mapping quality 0'),
    'reads QC failed': ('Reads QC-Failed [Samtools]', 'Number of reads that failed the quality checks'),
    'non-primary alignments': ('Non-Primary Alignments [Samtools]', 'Number of secondary reads'),
    'supplementary alignments': ('Supplementary Alignments [Samtools]', 'Number of supplementary reads'),
    'pairs on different chromosomes': ('Pairs on Different Chromosomes [Samtools]', 'Number of pairs where one read is on one chromosome and the mate read is on a different chromosome'),
    'percentage of properly paired reads (%)': ('Percentage of Properly Paired Reads [Samtools]', None)
}

picard_CollectAlignmentSummaryMetrics = {
    'PF_ALIGNED_BASES': ('Aligned Bases [Picard]', 'The total number of aligned bases'),
    'PF_HQ_ALIGNED_BASES': ('Aligned Bases (High Quality) [Picard]', 'The number of aligned bases in reads that were mapped at high quality'),
    'PF_MISMATCH_RATE': ('Aligned Bases Mismatch Rate [Picard]', 'The fraction of bases mismatching the reference for all aligned bases'),
    'PF_HQ_ERROR_RATE': ('Aligned Bases Mismatch Rate (High Quality) [Picard]', 'The fraction of bases mismatching the reference in reads that were mapped at high quality'),
    'PF_INDEL_RATE': ('Indel Rate [Picard]', 'The number of insertion and deletion events per 100 aligned bases'),
    'MEAN_READ_LENGTH': ('Mean Read Length [Picard]', 'The mean length of the set of reads examined'),
    'SD_READ_LENGTH': ('Read Length Standard Deviation [Picard]', 'The standard deviation for the length of the set of reads examined')
}

picard_CollectInsertSizeMetrics = {
    'MEAN_INSERT_SIZE': ('Mean Insert Size', 'The mean insert size for the pair orientation'), # [5]
    'STANDARD_DEVIATION': ('Insert Size Standard Deviation', 'Standard deviation of insert sizes for the pair orientation'), # [6]
    'READ_PAIRS': ('Total Number of Read Pairs', 'The total number of read pairs that were examined for the pair orientation') # [7]
    # 'PAIR_ORIENTATION' [8]
}

picard_CollectWgsMetrics = {
    'GENOME_TERRITORY': ('Effective Genome Size [Picard]', 'The number of non-N bases in the genome'),
    'MEAN_COVERAGE': ('Mean Coverage [Picard]', 'The mean coverage of the genome'),
    'SD_COVERAGE': ('Coverage Standard Deviation [Picard]', 'The standard deviation for the coverage')
}

########################################################################
# QMGeneric
#   -> QMValue
########################################################################
class QMGeneric:
    """
    """

    def __init__(self, name):
        """
        """
        self.name = name
        # List of QMValue objects
        self.qm_values = []

    class QMValue:
        """
        """

        def __init__(
                    self,
                    key,
                    value,
                    visible=None,
                    derived_from=None,
                    flag=None,
                    tooltip=None
                    ):
            """
            """
            # Required
            self.key = key
            self.value = value
            # Optional
            self.visible = visible
            self.derived_from = derived_from
            self.flag = flag
            self.tooltip = tooltip

        def to_dict(self):
            """
            """
            return { k:v for (k,v) in self.__dict__.items() if v != None }

    def add_value(self, QMValue_):
        """
        """
        self.qm_values.append(QMValue_)

    def to_dict(self):
        """
        """
        return {
            'name': self.name,
            'qc_values': [ QMValue_.to_dict() for QMValue_ in self.qm_values ]
        }

    def create_archive(self, filenames, archive_name=ARCHIVE_NAME):
        """
        """
        with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for filename in filenames:
                archive.write(filename)

    def write_json(self, filename, dict=None):
        """
        """
        if not dict:
            dict = self.to_dict()
        with open(filename, 'w') as fo:
            json.dump(dict, fo, sort_keys=True, indent=2)

########################################################################
# BamQM
#
#   BAM Quality Metrics [BamQM]
#       Samtools stats:
#           - SAMTOOLS_stats_OUTPUT
#       Samtools flagstat:
#           - SAMTOOLS_flagstat_OUTPUT
#       Samtools idxstats:
#           - SAMTOOLS_idxstats_OUTPUT
#       Picard CollectAlignmentSummaryMetrics:
#           - PICARD_CollectAlignmentSummaryMetrics_OUTPUT
#       Picard CollectBaseDistributionByCycle:
#           - PICARD_CollectBaseDistributionByCycle_OUTPUT
#           - PICARD_CollectBaseDistributionByCycle_PDF
#       Picard CollectGcBiasMetrics:
#           - PICARD_CollectGcBiasMetrics_OUTPUT
#           - PICARD_CollectGcBiasMetrics_SUMMARY
#           - PICARD_CollectGcBiasMetrics_PDF
#       Picard CollectInsertSizeMetrics:
#           - PICARD_CollectInsertSizeMetrics_OUTPUT
#           - PICARD_CollectInsertSizeMetrics_PDF
#       Picard CollectWgsMetrics:
#           - PICARD_CollectWgsMetrics_OUTPUT
#       Picard MeanQualityByCycle:
#           - PICARD_MeanQualityByCycle_OUTPUT
#           - PICARD_MeanQualityByCycle_PDF
#
########################################################################
class BamQM(QMGeneric):
    """
    """
    # Variables
    FILENAMES = [SAMTOOLS_stats_OUTPUT,
                 SAMTOOLS_flagstat_OUTPUT,
                 SAMTOOLS_idxstats_OUTPUT,
                 PICARD_CollectAlignmentSummaryMetrics_OUTPUT,
                 PICARD_CollectBaseDistributionByCycle_OUTPUT,
                 PICARD_CollectBaseDistributionByCycle_PDF,
                 PICARD_CollectGcBiasMetrics_OUTPUT,
                 PICARD_CollectGcBiasMetrics_SUMMARY,
                 PICARD_CollectGcBiasMetrics_PDF,
                 PICARD_CollectInsertSizeMetrics_OUTPUT,
                 PICARD_CollectInsertSizeMetrics_PDF,
                 PICARD_CollectWgsMetrics_OUTPUT,
                 PICARD_MeanQualityByCycle_OUTPUT,
                 PICARD_MeanQualityByCycle_PDF]

    def __init__(self):
        """
        """
        super().__init__(BAMQM_NAME)
        # Run Parsers
        self.parse_samtools_stats()
        self.parse_picard_CollectAlignmentSummaryMetrics()
        self.parse_picard_CollectInsertSizeMetrics()
        self.parse_picard_CollectWgsMetrics()
        # Write JSON object
        self.write_json(VALUES_FILENAME)
        #Create Archive
        self.create_archive(self.FILENAMES)

    def parse_samtools_stats(self):
        """
        """
        # Parse file and save values
        with open(SAMTOOLS_stats_OUTPUT) as fi:
            for line in fi:
                if line.startswith('SN'):
                    line = line.rstrip().split('\t')
                    field, value = line[1].replace(':', ''), line[2]
                    if field in samtools_stats:
                        field, tooltip_ = samtools_stats[field]
                        self.add_value(self.QMValue(field, value, tooltip=tooltip_))

    def parse_picard_CollectAlignmentSummaryMetrics(self):
        """
        """
        header, pair = [], []
        with open(PICARD_CollectAlignmentSummaryMetrics_OUTPUT) as fi:
            for line in fi:
                if line.startswith('CATEGORY'):
                    header = line.rstrip().split('\t')
                elif line.startswith('PAIR'):
                    pair = line.rstrip().split('\t')

        for i, field in enumerate(header):
            if field in picard_CollectAlignmentSummaryMetrics:
                field, tooltip_ = picard_CollectAlignmentSummaryMetrics[field]
                self.add_value(self.QMValue(field, pair[i], tooltip=tooltip_))

    def parse_picard_CollectInsertSizeMetrics(self):
        """
        """
        header, pairs, is_block = [], [], False
        with open(PICARD_CollectInsertSizeMetrics_OUTPUT) as fi:
            for line in fi:
                line = line.rstrip()
                if line:
                    if line.startswith('## HISTOGRAM'):
                        break
                    elif line.startswith('MEDIAN_INSERT_SIZE'):
                        is_block = True
                        header = line.split('\t')
                    elif is_block:
                        line = line.split('\t')
                        pairs.append(line)

        for pair in pairs:
            orientation = pair[8]
            for i, field in enumerate(header):
                if field in picard_CollectInsertSizeMetrics:
                    field, tooltip_ = picard_CollectInsertSizeMetrics[field]
                    field += f' ({orientation}) [Picard]'
                    self.add_value(self.QMValue(field, pair[i], tooltip=tooltip_))

    def parse_picard_CollectWgsMetrics(self):
        """
        """
        header, stats, is_block = [], [], False
        with open(PICARD_CollectWgsMetrics_OUTPUT) as fi:
            for line in fi:
                line = line.rstrip()
                if line:
                    if line.startswith('## HISTOGRAM'):
                        break
                    elif line.startswith('GENOME_TERRITORY'):
                        is_block = True
                        header = line.split('\t')
                    elif is_block:
                        stats = line.split('\t')

        for i, field in enumerate(header):
            if field in picard_CollectWgsMetrics:
                field, tooltip_ = picard_CollectWgsMetrics[field]
                self.add_value(self.QMValue(field, stats[i], tooltip=tooltip_))
