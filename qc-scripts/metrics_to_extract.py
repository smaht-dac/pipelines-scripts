########################################################################
# Metrics to extract for each tool
########################################################################
samtools = {
    'raw total sequences': {
        'key': 'Total Sequences [Samtools]',
        'tooltip': 'Total number of reads in a file, excluding supplementary and secondary reads'
    },
    'sequences': {
        'key': 'Processed Sequences [Samtools]',
        'tooltip': 'Number of processed reads'
    },
    '1st fragments': {
        'key': '1st Fragments [Samtools]', 
        'tooltip': 'Number of first fragment reads',
    },
    'last fragments': {
        'key': 'Last Fragments [Samtools]', 
        'tooltip': 'Number of last fragment reads',
    },
    'reads paired': {
        'key': 'Reads Paired [Samtools]', 
        'tooltip': 'Number of paired reads, mapped or unmapped, that are neither secondary nor supplementary (paired-end technology bit set)',
    },
    'reads mapped': {
        'key': 'Reads Mapped [Samtools]', 
        'tooltip': 'Number of reads, paired or single, that are mapped',
    },
    'reads mapped and paired': {
        'key': 'Reads Mapped and Paired [Samtools]', 
        'tooltip': 'Number of mapped paired reads (paired-end technology bit set with both mates mapped)',
    },
    'reads properly paired': {
        'key': 'Reads Properly Paired [Samtools]', 
        'tooltip': 'Number of mapped paired reads (proper-pair bit set, both mates mapped within the expected distance)',
    },
    'reads unmapped': {
        'key': 'Reads Unmapped [Samtools]', 
        'tooltip': 'Number of unmapped reads',
    },
    'reads duplicated': {
        'key': 'Reads Duplicated [Samtools]', 
        'tooltip': 'Number of duplicate reads',
    },
    'reads MQ0': {
        'key': 'Reads MQ0 [Samtools]', 
        'tooltip': 'Number of mapped reads with mapping quality 0',
    },
    'reads QC failed': {
        'key': 'Reads QC-Failed [Samtools]', 
        'tooltip': 'Number of reads that failed the quality checks',
    },
    'non-primary alignments': {
        'key': 'Non-Primary Alignments [Samtools]', 
        'tooltip': 'Number of secondary reads',
    },
    'supplementary alignments': {
        'key': 'Supplementary Alignments [Samtools]', 
        'tooltip': 'Number of supplementary reads',
    },
    'pairs on different chromosomes': {
        'key': 'Pairs on Different Chromosomes [Samtools]', 
        'tooltip': 'Number of pairs where one read is on one chromosome and the mate read is on a different chromosome',
    },
    'percentage of properly paired reads (%)': {
        'key': 'Percentage of Properly Paired Reads [Samtools]',
        'tooltip': None
    },
}

picard_CollectAlignmentSummaryMetrics = {
    'PF_ALIGNED_BASES': {
        'key': 'Aligned Bases [Picard]', 
        'tooltip': 'The total number of aligned bases',
    },
    'PF_HQ_ALIGNED_BASES': {
        'key': 'Aligned Bases (High Quality) [Picard]', 
        'tooltip': 'The number of aligned bases in reads that were mapped at high quality',
    },
    'PF_MISMATCH_RATE': {
        'key': 'Aligned Bases Mismatch Rate [Picard]', 
        'tooltip': 'The fraction of bases mismatching the reference for all aligned bases',
    },
    'PF_HQ_ERROR_RATE': {
        'key': 'Aligned Bases Mismatch Rate (High Quality) [Picard]', 
        'tooltip': 'The fraction of bases mismatching the reference in reads that were mapped at high quality',
    },
    'PF_INDEL_RATE': {
        'key': 'Indel Rate [Picard]', 
        'tooltip':  'The number of insertion and deletion events per 100 aligned bases',
    },
    'MEAN_READ_LENGTH': {
        'key': 'Mean Read Length [Picard]', 
        'tooltip': 'The mean length of the set of reads examined',
    },
    'SD_READ_LENGTH': {
        'key': 'Read Length Standard Deviation [Picard]', 
        'tooltip': 'The standard deviation for the length of the set of reads examined'
    },
}

picard_CollectInsertSizeMetrics = {
    # [5]
    'MEAN_INSERT_SIZE': {
        'key': 'Mean Insert Size', 
        'tooltip': 'The mean insert size for the pair orientation',
    },
    # [6]
    'STANDARD_DEVIATION': {
        'key': 'Insert Size Standard Deviation', 
        'tooltip': 'Standard deviation of insert sizes for the pair orientation',
    },
    # [7]
    'READ_PAIRS': {
        'key': 'Total Number of Read Pairs', 
        'tooltip': 'The total number of read pairs that were examined for the pair orientation'
    },
    # 'PAIR_ORIENTATION' [8]
}

picard_CollectWgsMetrics = {
    'GENOME_TERRITORY': {
        'key': 'Effective Genome Size [Picard]', 
        'tooltip': 'The number of non-N bases in the genome',
    },
    'MEAN_COVERAGE': {
        'key': 'Mean Coverage [Picard]', 
        'tooltip': 'The mean coverage of the genome',
    },
    'SD_COVERAGE': {
        'key': 'Coverage Standard Deviation [Picard]', 
        'tooltip': 'The standard deviation for the coverage'
    },
}

bamstats = {
    'Estimate_Average_Coverage' : {
        'key': 'Average coverage (estimated)', 
        'tooltip': 'Estimated average coverage',
    }
}


metrics = {
    'samtools': samtools,
    'bamstats': bamstats,
    'picard_CollectAlignmentSummaryMetrics': picard_CollectAlignmentSummaryMetrics,
    'picard_CollectInsertSizeMetrics': picard_CollectInsertSizeMetrics,
    'picard_CollectWgsMetrics': picard_CollectWgsMetrics
}
