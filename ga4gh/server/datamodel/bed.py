"""
Module responsible for translating bed sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import urllib
import gzip
import bz2

# no step/span; requires numpy
import pyBigWig


import ga4gh.server.datamodel as datamodel
import ga4gh.server.exceptions as exceptions
import ga4gh.server.gff3 as gff3
import ga4gh.schemas.pb as pb
import ga4gh.schemas.protocol as protocol

#TODO
#1. write bed parser and create bed parser test file
#1a. Find several example files (merge oddities for a test file)
#2. do read from bigBed 
#3. gtf
# - should be able to handle a bigBed store, but would need
#   to change access patterns and store type of sequence annotation
#   (out of scope?)


class ValidationException(Exception):
    """
    Exception caused by unexpected formatting in data files.
    """
    def __init__(self, message, fileName=None, lineNumber=None):
        """
        :param message: The error message to emit
        :param fileName: File being processed
        :param lineNumber: line in file, can be int or string
        """
        if fileName is not None:
            if lineNumber is not None:
                message = "{}:{}: {}".format(fileName, lineNumber, message)
            else:
                message = "{}: {}".format(fileName, message)
        super(ValidationException, self).__init__(message)


class Features(object):
    """
    A set of sequence annotations, accessed through byFeatureName, which
    is a dictionary of Feature, keyed by Feature.name.
    """
    def __init__(self):
        # index of features by id. GFF3 allows disjoint features with
        # the same id.  None is used to store features without ids
        self.byFeatureName = collections.defaultdict(list)

    def add(self, feature):
        """
        Add a feature record by featureName (which may be None)

        :param feature: Feature object being added.
        """
        self.byFeatureName[feature.featureName].append(feature)

class BedParser(object):
    """
    Parses a BED file into a BedSet. Performs basic validation.
    """

    def __init__(self, fileName, entryType='mRNA', blockType='CDS',
                 fivePrimeType='five_prime_UTR', 
                 threePrimeType='three_primeUTR'):
        """
        :param str fileName: Name of GFF3 file to parse,
            compressed files (.gz or .bz2) will be automatically decompressed.
        """
        self._fileName = fileName
        self._lineNumber = 0
        self._features = Features()
        self._type = entryType
        self._blockType = blockType
        self._fivePrimeType = fivePrimeType
        self._threePrimeType = threePrimeType

    def _open(self):
        """
        open input file, optionally with decompression
        """
        if self._fileName.endswith(".gz"):
            return gzip.open(self._fileName)
        elif self._fileName.endswith(".bz2"):
            return bz2.BZ2File(self._fileName)
        else:
            return open(self._fileName)

    def _addBlocks(self, blockCount, blockSizesStr, blockStartsStr, 
                   thickStart, thickEnd, feature):
        """
        Parse bed blocks and add children features for each block
        """
        blockStarts = blockStartsStr.split(',')
        if blockStarts.length != blockCount:
            raise ValidationException(
                "Number of block starts doesn't match block count"
                self._fileName, self._lineNumber)
        blockEnds = blockEndsStr.split(',')
        if blockEnds.length != blockCount:
            raise ValidationException(
                "Number of block ends doesn't match block count"
                self._fileName, self._lineNumber)
        for i in range(0, blockCount):
            try:
                start = feature.start + int(blockStarts[i])
                end = start + int(blockSizes[i])
            except ValueError:
                raise ValidationException(
                    "Block start or size is not an integer",
                    self._fileName, self._lineNumber)
#thickness

            attributes = dict()
            childFeature = gff3.Feature(feature.refName, feature.source, 
                                        self._blockType, start, end,
                                        feature.score, feature.strand, 
                                        feature.strand, attributes)
            childFeature.parents.add(feature) 
            feature.children.add(childFeature)
            self._features.add(childFeature)

    def _parseLine(self, line):
        """ 
        Parse one line of the bed file.
        """
        if line.startswith('track') or line.startswith('browser'):
            return

        fields = line.split()
        if fields.length < 3:
            raise ValidationException(
                "Line has less than three fields",
                self._fileName, self._lineNumber)
        (refName, start, end, name, score, strand, thickStart, thickEnd,
                itemRGB, blockCount, blockSizesStr, blockStartsStr) = fields
        refName = urlib.unquote(refName)
        try:
            start = int(start)
            end = int(end)
        except ValueError:
            raise ValidationException(
                "Start or end is not an integer",
                self._fileName, self._lineNumber)
        source = None 
        frame = 0
        attributes = dict()
        if name:
            attributes['ID'] = name
        feature = gff3.Feature(refName, source, self._type, start, end,
                               score, strand, frame, attributes)
        if blockCount:
            self._addBlocks(blockCount, blockSizesStr, blockStartsStr,
                            thickStart, thickEnd, feature)
        self._features.add(feature)

    def parse(self):
        """
        Run the parse and return the resulting Features object.
        """
        fh = self._open()
        try:
            for line in fh:
                self._lineNumber += 1
                line = line.rstrip()
                if line:
                    self._parseLine(line)
        finally:
            fh.close()
        return _features
