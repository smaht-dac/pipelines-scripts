########################################################################
# Metrics to extract for each tool
########################################################################
samtools = {
    'raw total sequences': {
        'name': 'Total Sequences [Samtools]',
        'tooltip': 'Total number of reads in a file, excluding supplementary and secondary reads'
    },
    'sequences': {
        'name': 'Processed Sequences [Samtools]',
        'tooltip': 'Number of processed reads'
    },
    '1st fragments': {
        'name': '1st Fragments [Samtools]', 
        'tooltip': 'Number of first fragment reads',
    },
    'last fragments': {
        'name': 'Last Fragments [Samtools]', 
        'tooltip': 'Number of last fragment reads',
    },
    'reads paired': {
        'name': 'Reads Paired [Samtools]', 
        'tooltip': 'Number of paired reads, mapped or unmapped, that are neither secondary nor supplementary (paired-end technology bit set)',
    },
    'reads mapped': {
        'name': 'Reads Mapped [Samtools]', 
        'tooltip': 'Number of reads, paired or single, that are mapped',
    },
    'reads mapped and paired': {
        'name': 'Reads Mapped and Paired [Samtools]', 
        'tooltip': 'Number of mapped paired reads (paired-end technology bit set with both mates mapped)',
    },
    'reads properly paired': {
        'name': 'Reads Properly Paired [Samtools]', 
        'tooltip': 'Number of mapped paired reads (proper-pair bit set, both mates mapped within the expected distance)',
    },
    'reads unmapped': {
        'name': 'Reads Unmapped [Samtools]', 
        'tooltip': 'Number of unmapped reads',
    },
    'reads duplicated': {
        'name': 'Reads Duplicated [Samtools]', 
        'tooltip': 'Number of duplicate reads',
    },
    'reads MQ0': {
        'name': 'Reads MQ0 [Samtools]', 
        'tooltip': 'Number of mapped reads with mapping quality 0',
    },
    'reads QC failed': {
        'name': 'Reads QC-Failed [Samtools]', 
        'tooltip': 'Number of reads that failed the quality checks',
    },
    'non-primary alignments': {
        'name': 'Non-Primary Alignments [Samtools]', 
        'tooltip': 'Number of secondary reads',
    },
    'supplementary alignments': {
        'name': 'Supplementary Alignments [Samtools]', 
        'tooltip': 'Number of supplementary reads',
    },
    'pairs on different chromosomes': {
        'name': 'Pairs on Different Chromosomes [Samtools]', 
        'tooltip': 'Number of pairs where one read is on one chromosome and the mate read is on a different chromosome',
    },
    'percentage of properly paired reads (%)': {
        'name': 'Percentage of Properly Paired Reads [Samtools]',
        'tooltip': None
    },
}

picard_CollectAlignmentSummaryMetrics = {
    'PF_ALIGNED_BASES': {
        'name': 'Aligned Bases [Picard]', 
        'tooltip': 'The total number of aligned bases',
    },
    'PF_HQ_ALIGNED_BASES': {
        'name': 'Aligned Bases (High Quality) [Picard]', 
        'tooltip': 'The number of aligned bases in reads that were mapped at high quality',
    },
    'PF_MISMATCH_RATE': {
        'name': 'Aligned Bases Mismatch Rate [Picard]', 
        'tooltip': 'The fraction of bases mismatching the reference for all aligned bases',
    },
    'PF_HQ_ERROR_RATE': {
        'name': 'Aligned Bases Mismatch Rate (High Quality) [Picard]', 
        'tooltip': 'The fraction of bases mismatching the reference in reads that were mapped at high quality',
    },
    'PF_INDEL_RATE': {
        'name': 'Indel Rate [Picard]', 
        'tooltip':  'The number of insertion and deletion events per 100 aligned bases',
    },
    'MEAN_READ_LENGTH': {
        'name': 'Mean Read Length [Picard]', 
        'tooltip': 'The mean length of the set of reads examined',
    },
    'SD_READ_LENGTH': {
        'name': 'Read Length Standard Deviation [Picard]', 
        'tooltip': 'The standard deviation for the length of the set of reads examined'
    },
}

picard_CollectInsertSizeMetrics = {
    # [5]
    'MEAN_INSERT_SIZE': {
        'name': 'Mean Insert Size', 
        'tooltip': 'The mean insert size for the pair orientation',
    },
    # [6]
    'STANDARD_DEVIATION': {
        'name': 'Insert Size Standard Deviation', 
        'tooltip': 'Standard deviation of insert sizes for the pair orientation',
    },
    # [7]
    'READ_PAIRS': {
        'name': 'Total Number of Read Pairs', 
        'tooltip': 'The total number of read pairs that were examined for the pair orientation'
    },
    # 'PAIR_ORIENTATION' [8]
}

picard_CollectWgsMetrics = {
    'GENOME_TERRITORY': {
        'name': 'Effective Genome Size [Picard]', 
        'tooltip': 'The number of non-N bases in the genome',
    },
    'MEAN_COVERAGE': {
        'name': 'Mean Coverage [Picard]', 
        'tooltip': 'The mean coverage of the genome',
    },
    'SD_COVERAGE': {
        'name': 'Coverage Standard Deviation [Picard]', 
        'tooltip': 'The standard deviation for the coverage'
    },
}


metrics = {
    'samtools': samtools,
    'picard_CollectAlignmentSummaryMetrics': picard_CollectAlignmentSummaryMetrics,
    'picard_CollectInsertSizeMetrics': picard_CollectInsertSizeMetrics,
    'picard_CollectWgsMetrics': picard_CollectWgsMetrics
}
