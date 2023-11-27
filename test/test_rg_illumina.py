import sys, os
import pytest
import pysam

# setting path for pytest
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

# import code to test
from processing_scripts.AddReadGroups import get_read_group_illumina
from processing_scripts.AddReadGroups import main

# TEST func get_read_group_illumina(...)
def test_get_read_group_illumina_NEW():
    expected_QNAME = set([
    'MG01HX02_343_HTCVYCCXX.2',
    'MG01HX02_341_HTCVYCCXX.2',
    'MG01HX02_343_HTCVYCCXX.3'
    ])
    samfile = pysam.AlignmentFile('test/files/rg_illumina_NEW_test.sam', 'r')
    for read in samfile:
        QNAME = get_read_group_illumina(read.query_name)
        assert QNAME in expected_QNAME

def test_get_read_group_illumina_OLD():
    expected_QNAME = set([
    'HWUSI-EAS100R.5',
    'HWUSI-EAS100R.6',
    'HWUSI-EAS100R.7'
    ])
    samfile = pysam.AlignmentFile('test/files/rg_illumina_OLD_test.sam', 'r')
    for read in samfile:
        QNAME = get_read_group_illumina(read.query_name)
        assert QNAME in expected_QNAME

def test_get_read_group_illumina_MOD():
    expected_QNAME = set([
    'READGROUP'
    ])
    samfile = pysam.AlignmentFile('test/files/rg_illumina_MOD_test.sam', 'r')
    for read in samfile:
        QNAME = get_read_group_illumina(read.query_name)
        assert QNAME in expected_QNAME

def test_get_read_group_illumina_ERR():
        samfile = pysam.AlignmentFile('test/files/rg_illumina_ERROR_test.sam', 'r')
        with pytest.raises(SystemExit) as e:
            for read in samfile:
                QNAME = get_read_group_illumina(read.query_name)
                assert str(e.value) == '\nFORMAT ERROR: read name format MG01HX02:343:HTCVYCCXX:2:1208:20354:69168 1:N:18:1 not recognized\n'
