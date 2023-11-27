#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################
#
#   Split FASTQ input file by lanes.
#   Reads per each lane will be written in a different file.
#   The script will generate compressed FASTQ files.
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################


################################################
#   Libraries
################################################
import sys, argparse, os
import gzip


################################################
#   Objects
################################################
class FASTQParser(object):

    def __init__(self, inputfile):
        '''
        '''
        self.inputfile = inputfile
        self.parse_reads()

    class FASTQRead(object):

        def __init__(self, id, sequence, description, quality):
            '''
            '''
            self.id = id
            self.sequence = sequence
            self.description = description
            self.quality = quality

        def to_string(self):
            ''' return the read in FASTQ format as string
            '''
            return f'{self.id}\n{self.sequence}\n{self.description}\n{self.quality}\n'

        def get_lane_illumina(self):
            ''' get lane information for illumina data in the format:
                    - <instrument>_<run number>_<flowcell ID>_<lane> [new]
                    - <flowcell>_<lane> [legacy]
            '''
            query_name = self.id.split()[0].replace('@', '')
            qn_as_list = query_name.split(':')
            l = len(qn_as_list) # Number of fields in query_name
            if l == 7:
                # New style header
                return '_'.join(qn_as_list[:4])
            elif l == 5:
                # Old style header
                return '_'.join(qn_as_list[:2])
            else:
                sys.exit('\nFORMAT ERROR: read format {0} not recognized\n'
                            .format(query_name))

    def read_fastq(self):
        ''' read FASTQ file, gzipped or uncompressed,
        return a generator '''
        if self.inputfile.endswith('.gz') or \
           self.inputfile.endswith('.bgz'):
            with gzip.open(self.inputfile, 'rb') as fz:
                for byteline in fz:
                    yield byteline.decode('utf-8')
        else:
            with open(self.inputfile, encoding='utf-8') as fi:
                for line in fi:
                    yield line

    def parse_reads(self):
        ''' parse FASTQ multi-line format
        '''
        i = 0
        for line in self.read_fastq():
            i += 1
            # Parse FASTQ read structure
            if i == 1:
                id = line.rstrip()
            elif i == 2:
                sequence = line.rstrip()
            elif i == 3:
                description = line.rstrip()
            elif i == 4:
                quality = line.rstrip()
                yield self.FASTQRead(id, sequence, description, quality)
                i = 0


################################################
#   Functions
################################################
def main(args):

    inputfile = args['inputfile']
    output_lanes = []
    output_buffers = []

    fasqt_parser = FASTQParser(inputfile)

    for i, read_obj in enumerate(fasqt_parser.parse_reads()):
        lane = read_obj.get_lane_illumina()
        if lane in output_lanes:
            idx = output_lanes.index(lane)
            output_buffers[idx].write(read_obj.to_string().encode(encoding='utf-8'))
        else:
            filename_ = inputfile.split('/')[-1].split('.')[0] + '_' + lane.replace('@', '') + '.fastq.gz'
            buffer_ = gzip.open(filename_, 'wb')
            output_lanes.append(lane)
            output_buffers.append(buffer_)
            buffer_.write(read_obj.to_string().encode(encoding='utf-8'))

    for buffer_ in output_buffers:
        buffer_.close()


################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Split FASTQ file by lanes')

    parser.add_argument('-i','--inputfile', help='input FASTQ file', required=True)

    args = vars(parser.parse_args())

    main(args)
