########################################################################
# Metrics to extract for each tool
########################################################################
samtools = {
    'raw total sequences': {
        'key': 'Total Sequences [Samtools]',
        'tooltip': 'Total number of reads in a file, excluding supplementary and secondary reads',
        'type': int
    },
    'sequences': {
        'key': 'Processed Sequences [Samtools]',
        'tooltip': 'Number of processed reads',
        'type': int
    },
    '1st fragments': {
        'key': '1st Fragments [Samtools]',
        'tooltip': 'Number of first fragment reads',
        'type': int
    },
    'last fragments': {
        'key': 'Last Fragments [Samtools]',
        'tooltip': 'Number of last fragment reads',
        'type': int
    },
    'reads paired': {
        'key': 'Reads Paired [Samtools]',
        'tooltip': 'Number of paired reads, mapped or unmapped, that are neither secondary nor supplementary (paired-end technology bit set)',
        'type': int
    },
    'reads mapped': {
        'key': 'Reads Mapped [Samtools]',
        'tooltip': 'Number of reads, paired or single, that are mapped',
        'type': int
    },
    'reads mapped and paired': {
        'key': 'Reads Mapped and Paired [Samtools]',
        'tooltip': 'Number of mapped paired reads (paired-end technology bit set with both mates mapped)',
        'type': int
    },
    'reads properly paired': {
        'key': 'Reads Properly Paired [Samtools]',
        'tooltip': 'Number of mapped paired reads (proper-pair bit set, both mates mapped within the expected distance)',
        'type': int
    },
    'reads unmapped': {
        'key': 'Reads Unmapped [Samtools]',
        'tooltip': 'Number of unmapped reads',
        'type': int
    },
    'reads duplicated': {
        'key': 'Reads Duplicated [Samtools]',
        'tooltip': 'Number of duplicate reads',
        'type': int
    },
    'reads MQ0': {
        'key': 'Reads MQ0 [Samtools]',
        'tooltip': 'Number of mapped reads with mapping quality 0',
        'type': int
    },
    'reads QC failed': {
        'key': 'Reads QC-Failed [Samtools]',
        'tooltip': 'Number of reads that failed the quality checks',
        'type': int
    },
    'non-primary alignments': {
        'key': 'Non-Primary Alignments [Samtools]',
        'tooltip': 'Number of secondary reads',
        'type': int
    },
    'supplementary alignments': {
        'key': 'Supplementary Alignments [Samtools]',
        'tooltip': 'Number of supplementary reads',
        'type': int
    },
    'pairs on different chromosomes': {
        'key': 'Pairs on Different Chromosomes [Samtools]',
        'tooltip': 'Number of pairs where one read is on one chromosome and the mate read is on a different chromosome',
        'type': int
    },
    'percentage of properly paired reads (%)': {
        'key': 'Percentage of Properly Paired Reads [Samtools]',
        'tooltip': None,
        'type': float
    },
}

picard_CollectAlignmentSummaryMetrics = {
    'PF_ALIGNED_BASES': {
        'key': 'Aligned Bases [Picard]',
        'tooltip': 'The total number of aligned bases',
        'type': int
    },
    'PF_HQ_ALIGNED_BASES': {
        'key': 'Aligned Bases (High Quality) [Picard]',
        'tooltip': 'The number of aligned bases in reads that were mapped at high quality',
        'type': int
    },
    'PF_MISMATCH_RATE': {
        'key': 'Aligned Bases Mismatch Rate [Picard]',
        'tooltip': 'The fraction of bases mismatching the reference for all aligned bases',
        'type': float
    },
    'PF_HQ_ERROR_RATE': {
        'key': 'Aligned Bases Mismatch Rate (High Quality) [Picard]',
        'tooltip': 'The fraction of bases mismatching the reference in reads that were mapped at high quality',
        'type': float
    },
    'PF_INDEL_RATE': {
        'key': 'Indel Rate [Picard]',
        'tooltip':  'The number of insertion and deletion events per 100 aligned bases',
        'type': float
    },
    'MEAN_READ_LENGTH': {
        'key': 'Mean Read Length [Picard]',
        'tooltip': 'The mean length of the set of reads examined',
        'type': int
    },
    'SD_READ_LENGTH': {
        'key': 'Read Length Standard Deviation [Picard]',
        'tooltip': 'The standard deviation for the length of the set of reads examined',
        'type': int
    },
}

picard_CollectInsertSizeMetrics = {
    # [5]
    'MEAN_INSERT_SIZE': {
        'key': 'Mean Insert Size',
        'tooltip': 'The mean insert size for the pair orientation',
        'type': float
    },
    # [6]
    'STANDARD_DEVIATION': {
        'key': 'Insert Size Standard Deviation',
        'tooltip': 'Standard deviation of insert sizes for the pair orientation',
        'type': float
    },
    # [7]
    'READ_PAIRS': {
        'key': 'Total Number of Read Pairs',
        'tooltip': 'The total number of read pairs that were examined for the pair orientation',
        'type': int
    },
    # 'PAIR_ORIENTATION' [8]
}

picard_CollectWgsMetrics = {
    'GENOME_TERRITORY': {
        'key': 'Effective Genome Size [Picard]',
        'tooltip': 'The number of non-N bases in the genome',
        'type': int
    },
    'MEAN_COVERAGE': {
        'key': 'Mean Coverage [Picard]',
        'tooltip': 'The mean coverage of the genome',
        'type': float
    },
    'SD_COVERAGE': {
        'key': 'Coverage Standard Deviation [Picard]',
        'tooltip': 'The standard deviation for the coverage',
        'type': float
    },
}

bamstats = {
    'Estimate_Average_Coverage': {
        'key': 'Average coverage (estimated)',
        'tooltip': 'Estimated average coverage',
        'type': float
    }
}


metrics = {
    'samtools': samtools,
    'bamstats': bamstats,
    'picard_CollectAlignmentSummaryMetrics': picard_CollectAlignmentSummaryMetrics,
    'picard_CollectInsertSizeMetrics': picard_CollectInsertSizeMetrics,
    'picard_CollectWgsMetrics': picard_CollectWgsMetrics
}
