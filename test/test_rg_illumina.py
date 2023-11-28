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

# TEST end to end illumina
#   this requires samtools installed and will fail otherwise
def test_illumina():
    args = {
        'samplename': 'SAMPLE',
        'library': 'LIB',
        'inputfile': 'test/files/rg_illumina_NEW_test.sam',
        'directory': 'test/files/',
        'platform': None,
        'threads': None,
        'index': None
    }
    expected_header_RG = [
        {'ID': 'SAMPLE.MG01HX02_343_HTCVYCCXX.2', 'SM': 'SAMPLE', 'PL': 'ILLUMINA', 'PU': 'MG01HX02_343_HTCVYCCXX.2', 'LB': 'SAMPLE.LIB'},
        {'ID': 'SAMPLE.MG01HX02_343_HTCVYCCXX.3', 'SM': 'SAMPLE', 'PL': 'ILLUMINA', 'PU': 'MG01HX02_343_HTCVYCCXX.3', 'LB': 'SAMPLE.LIB'},
        {'ID': 'SAMPLE.MG01HX02_341_HTCVYCCXX.2', 'SM': 'SAMPLE', 'PL': 'ILLUMINA', 'PU': 'MG01HX02_341_HTCVYCCXX.2', 'LB': 'SAMPLE.LIB'}
    ]
    expected_RG = {
        'MG01HX02:343:HTCVYCCXX:2:1208:20354:69168': 'SAMPLE.MG01HX02_343_HTCVYCCXX.2',
        'MG01HX02:343:HTCVYCCXX:2:2101:26392:21966': 'SAMPLE.MG01HX02_343_HTCVYCCXX.2',
        'MG01HX02:341:HTCVYCCXX:2:1215:5365:64615': 'SAMPLE.MG01HX02_341_HTCVYCCXX.2',
        'MG01HX02:343:HTCVYCCXX:3:2124:28442:1080': 'SAMPLE.MG01HX02_343_HTCVYCCXX.3'
    }
    main(args) # AddReadGroups
    samfile = pysam.AlignmentFile('test/files/rg_illumina_NEW_test_rg.bam', 'rb')
    header_dict = dict(samfile.header)
    assert len(expected_header_RG) == len(header_dict['RG'])
    for rg in expected_header_RG:
        assert rg in header_dict['RG']
    for read in samfile:
        RG = read.get_tag('RG')
        QNAME = read.query_name
        assert expected_RG[QNAME] == RG
    os.remove('test/files/rg_illumina_NEW_test_rg.bam')
