"""
Module responsible for translating bed sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import random
import re
import math

# no step/span; requires numpy
import pyBigWig

# for running bigwig tool externally
import subprocess

import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions
import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol

TODO
1. write bed parser and create gff3 test file
2. do bigBed

class BedSet(object):
    """
    A set of GFF3 sequence annotations
    """
    def __init__(self, fileName=None):
        self.fileName = fileName
        self.roots = set()     # root nodes (those with out parents)
        # index of features by id. GFF3 allows disjoint features with
        # the same id.  None is used to store features without ids
        self.byFeatureName = collections.defaultdict(list)

class BedParser(object):
    """
    Parses a BED or BigBed file into a BedSet. Performs basic validation.
    """

    def __init__(self, fileName):
        """
        :param str fileName: Name of GFF3 file to parse,
            compressed files (.gz or .bz2) will be automatically decompressed.
        """
        self.fileName = fileName
        self.lineNumber = 0

    def _open(self):
        """
        open input file, optionally with decompression
        """
        if self.fileName.endswith(".gz"):
            return gzip.open(self.fileName)
        elif self.fileName.endswith(".bz2"):
            return bz2.BZ2File(self.fileName)
        else:
            return open(self.fileName)

    def _parseRecord(self, gff3Set, line):
        """
        Parse one record.
        """
        row = line.split("\t")
        if len(row) != self.GFF3_NUM_COLS:
            raise GFF3Exception(
                "Wrong number of columns, expected {}, got {}".format(
                    self.GFF3_NUM_COLS, len(row)),
                self.fileName, self.lineNumber)
        feature = Feature(
            urllib.unquote(row[0]),
            urllib.unquote(row[1]),
            urllib.unquote(row[2]),
            int(row[3]), int(row[4]),
            row[5], row[6], row[7],
            self._parseAttrs(row[8]))
        gff3Set.add(feature)

    # spaces or comment line
    IGNORED_LINE_RE = re.compile("(^[ ]*$)|(^[ ]*#.*$)")

    @staticmethod
    def _isIgnoredLine(line):
        return Gff3Parser.IGNORED_LINE_RE.search(line) is not None

    def _checkHeader(self, line):
        # split to allow multiple spaces and tabs
        if line.split() != GFF3_HEADER.split():
            raise GFF3Exception(
                "First line is not GFF3 header ({}), got: {}".format(
                    GFF3_HEADER, line), self.fileName, self.lineNumber)

    def _parseLine(self, gff3Set, line):
        if self.lineNumber == 1:
            self._checkHeader(line)
        elif not self._isIgnoredLine(line):
            self._parseRecord(gff3Set, line)

    def parse(self):
        """
        Run the parse and return the resulting Gff3Set object.
        """
        fh = self._open()
        try:
            bedSet = bedSet(self.fileName)
            for line in fh:
                self.lineNumber += 1
                self._parseLine(gff3Set, line[0:-1])
        finally:
            fh.close()
        return bedSet
