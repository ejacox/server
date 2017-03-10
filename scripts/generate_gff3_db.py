"""
Generate a SQLite DB file from a sequence annotation file in GFF3 format
(http://sequenceontology.org/resources/gff3.html)
The resulting database file will be mapped to a FeatureSet in the current
GA4GH server implementation, and each row in the 'feature' table will be
mapped to a Feature.
"""

from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys

import generate_feature_db
import ga4gh.common.utils as utils
import glue

glue.ga4ghImportGlue()


@utils.Timed()
def main():
    parser = argparse.ArgumentParser(
        description="Script to generate SQLite database corresponding to "
        "an input GFF3 sequence annotations file.")
    parser.add_argument(
        "--outputFile", "-o", default="annotations.db",
        help="The file to output the server-ready database to.")
    parser.add_argument(
        "--inputFile", "-i",
        help="Path to input GFF3 file.",
        default='.')
    parser.add_argument('--verbose', '-v', action='count', default=0)
    args = parser.parse_args()
    g2d = generate_feature_db.file2Db(args.inputFile, args.outputFile, 'gff3')
    g2d.run()


if __name__ == "__main__":
    main()
