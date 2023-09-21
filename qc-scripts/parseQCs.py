#!/usr/bin/env python

########################################################################
#
#   Authors:
#       Michele Berselli
#       Harvard Medical School
#       berselli.michele@gmail.com
#
#       Alexander Veit
#       Harvard Medical School
#       alexander_veit@hms.harvard.edu
#
#   Script to parse the output from several QC tools to generate
#       portal-compatible Quality Metric (QM) objects
#
########################################################################

import click
from QMGeneric import QMGeneric
from MetricsParser import Parser

########################################################################
# Global Variables
########################################################################
VALUES_FILENAME = 'qc_values.json'
ARCHIVE_NAME = 'metrics.zip'


@click.command()
@click.help_option("--help", "-h")
@click.option(
    "-n",
    "--qm-name",
    required=True,
    type=str,
    help="Name of the Quality Metric",
)
@click.option(
    "-m",
    "--metrics",
    required=True,
    type=str,
    nargs=2,
    multiple=True,
    help="QC tool and File that will be parsed",
)
@click.option(
    "-a",
    "--additional-files",
    required=False,
    type=str,
    multiple=True,
    help="Files that will be added to the zip archive",
)
def main(qm_name, metrics, additional_files):
    """
    This script gathers metrics from different tools and creates a QualityMetricGeneric Item for the portal.

    Example usage:
    python parseQCs.py \
        -n 'BAM Quality Metrics' \
        --metrics samtools /PATH/samtools.stats.txt \
        --metrics picard_CollectInsertSizeMetrics /PATH/picard_cis_metrics.txt \
        --additional-files /PATH/additional_output_1.pdf \
        --additional-files /PATH/additional_output_2.tsv \

    This command will parse samtools and picard_CollectInsertSizeMetrics metrics from the files provided for each
    tool and create a qc.json file. It will also create a zip archive with all 4 provided files.

    """
    all_files = additional_files if additional_files else []
    qmg = QMGeneric(qm_name)

    for m in metrics:
        tool_name, metrics_file = m
        all_files.append(metrics_file)

        parser = Parser(tool_name, metrics_file)
        qm_values = parser.parse()
        qmg.add_values(qm_values)

    # Write JSON object
    qmg.write_json(VALUES_FILENAME)
    # Create Archive
    qmg.create_archive(all_files, ARCHIVE_NAME)


if __name__ == "__main__":
    main()
