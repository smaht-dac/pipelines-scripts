import sys
from metrics_to_extract import metrics
from QMGeneric import QMValue
from typing import List
import json


class Parser:

    def __init__(self, tool : str, metrics_file : str):

        self.tool = tool
        self.path = metrics_file

    def parse(self):

        if self.tool == "samtools":
            return self.parse_samtools_stats()
        elif self.tool == "picard_CollectAlignmentSummaryMetrics":
            return self.parse_picard_CollectAlignmentSummaryMetrics()
        elif self.tool == "picard_CollectInsertSizeMetrics":
            return self.parse_picard_CollectInsertSizeMetrics()
        elif self.tool == "picard_CollectWgsMetrics":
            return self.parse_picard_CollectWgsMetrics()
        elif self.tool == "bamstats":
            return self.parse_bamstats()
        else:
            sys.exit(f"{self.tool} is not supported. Please add a parser to Parser.py")

    def parse_samtools_stats(self) -> List[QMValue]:
        qm_values = []
        # Parse file and save values
        with open(self.path) as fi:
            for line in fi:
                if line.startswith('SN'):
                    line = line.rstrip().split('\t')
                    field, value = line[1].replace(':', ''), line[2]
                    if field in metrics['samtools']:
                        m = metrics['samtools'][field]
                        m_type = m["type"]
                        qmv = QMValue(
                            m["key"], m_type(value), tooltip=m["tooltip"])
                        qm_values.append(qmv)
        return qm_values
    
    def parse_bamstats(self) -> List[QMValue]:
        qm_values = []
        # Parse file and save values
        fi = open(self.path)
        bamstats_res = json.load(fi)
        for key in bamstats_res.keys():
            if key in metrics['bamstats']:
                m = metrics['bamstats'][key]
                value = bamstats_res[key]
                m_type = m["type"]
                qmv = QMValue(
                    m["key"], m_type(value), tooltip=m["tooltip"])
                qm_values.append(qmv)
        fi.close()
        return qm_values

    def parse_picard_CollectAlignmentSummaryMetrics(self) -> List[QMValue]:
        qm_values = []
        header, pair = [], []
        with open(self.path) as fi:
            for line in fi:
                if line.startswith('CATEGORY'):
                    header = line.rstrip().split('\t')
                elif line.startswith('PAIR'):
                    pair = line.rstrip().split('\t')

        for i, field in enumerate(header):
            if field in metrics['picard_CollectAlignmentSummaryMetrics']:
                m = metrics['picard_CollectAlignmentSummaryMetrics'][field]
                m_type = m["type"]
                qmv = QMValue(
                    m["key"], m_type(pair[i]), tooltip=m["tooltip"])
                qm_values.append(qmv)
        return qm_values

    def parse_picard_CollectInsertSizeMetrics(self) -> List[QMValue]:
        qm_values = []
        header, pairs, is_block = [], [], False
        with open(self.path) as fi:
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
                if field in metrics['picard_CollectInsertSizeMetrics']:
                    m = metrics['picard_CollectInsertSizeMetrics'][field]
                    m_type = m["type"]
                    qmv = QMValue(
                        m["key"]+f' ({orientation}) [Picard]', m_type(pair[i]), tooltip=m["tooltip"])
                    qm_values.append(qmv)
        return qm_values

    def parse_picard_CollectWgsMetrics(self) -> List[QMValue]:
        qm_values = []
        header, stats, is_block = [], [], False
        with open(self.path) as fi:
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
            if field in metrics['picard_CollectWgsMetrics']:
                m = metrics['picard_CollectWgsMetrics'][field]
                m_type = m["type"]
                qmv = QMValue(
                    m["key"], m_type(stats[i]), tooltip=m["tooltip"])
                qm_values.append(qmv)
        return qm_values
