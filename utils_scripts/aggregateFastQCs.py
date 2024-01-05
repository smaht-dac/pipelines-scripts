#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################
#
#   Aggregate FastQC metrics.
#   Requires multiqc.
#
#   Michele Berselli
#   berselli.michele@gmail.com
#
################################################


################################################
#   Libraries
################################################
import sys, os
import argparse
import json
import boto3
import zipfile
import subprocess
from dcicutils import ff_utils, s3_utils

################################################
#   Functions
################################################
def get_credentials(keydicts_json, ff_env):
    """Get auth credentials.
    """
    # Get portal credentials
    if os.path.exists(keydicts_json):
        with open(os.path.expanduser(keydicts_json)) as keyfile:
            keys = json.load(keyfile)
        ff_key = keys.get(ff_env)
    elif os.environ.get('GLOBAL_ENV_BUCKET') and os.environ.get('S3_ENCRYPT_KEY'):
        s3 = s3_utils.s3Utils(env=ff_env)
        ff_key = s3.get_access_keys('access_key_admin')
    else:
        raise Exception('Required deployment vars GLOBAL_ENV_BUCKET and/or S3_ENCRYPT_KEY not set, and no entry for specified enivornment exists in keydicts file.')

    # Get encryption key
    kms_key_id = os.environ.get('S3_ENCRYPT_KEY_ID', None)

    return ff_key, kms_key_id

def get_download_link(qualitymetric_uuid, ff_key):
    """
    """
    return ff_utils.get_metadata(qualitymetric_uuid, add_on='frame=raw&datastore=database', key=ff_key)['url']

def aws_dowload_link(url, target_dir):
    """
    """
    bucket = url.split('.s3')[0].split('//')[-1]
    filename = '/'.join(url.split('/')[-2:])
    target_name = f'{target_dir}/{os.path.basename(url)}'
    s3 = boto3.client('s3')
    s3.download_file(bucket, filename, target_name)

    sys.stderr.write(f'Downloaded {target_name}\n')

    return target_name

def unzip(file_zip, target_dir):
    """
    """
    with zipfile.ZipFile(file_zip, 'r') as zip_ref:
        zip_ref.extractall(target_dir)
        sys.stderr.write(f'Unzipped {file_zip}\n')

    os.remove(file_zip)

def main(args):

    # <rdm_name>
    #   -> fastqc.report.zip
    #   -> fastqc.summary.txt
    #
    # fastqc.report.zip
    #   -> <ACCID>_fastqc
    #       -> fastqc_data.txt
    #       -> ...

    working_dir = args['working_dir']
    ff_key, _ = get_credentials(args['keydicts_json'], args['ff_env'])

    for uuid in args['input_uuids']:
        url = get_download_link(uuid, ff_key)
        download_zip = aws_dowload_link(url, working_dir)
        unzip(download_zip, working_dir)

    for rndm in os.listdir(working_dir):
        rndm_ = f'{working_dir}/{rndm}'
        unzip(f'{rndm_}/fastqc.report.zip', rndm_)

    # Run multiqc
    subprocess.run([ 'multiqc', working_dir, '--title', args['output_name'], '--no-data-dir', '--outdir', working_dir ])


################################################
#   MAIN
################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Aggregate FastQC metrics in a multiqc report')

    parser.add_argument('-i','--input_uuids', nargs='+', help='List of FastQC QualityMetric uuid', required=True)
    parser.add_argument('--ff_env', help='Environment to use', required=True)
    parser.add_argument('-o','--output_name', default='metrics', help='Name of the output [metrics]', required=False)
    parser.add_argument('-d','--working_dir', default='.', help='Working directory [.]', required=False)
    parser.add_argument('--keydicts_json', default='~/.smaht-keys.json',
                        help='Path to file with keys for portal auth in JSON format [~/.smaht-keys.json]', required=False)

    args = vars(parser.parse_args())

    main(args)
