"""
Module responsible for translating sequence annotation data
into GA4GH native objects.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import json
import random
import re
import math

# bx-python: bigwig reader missing from pypi version
#bx-python==0.7.3
#from bx.arrays.array_tree import *
#import bx.wiggle
#import bx.bbi.bigwig_file

# not being maintained; difficult to install (.so path problem, currently)
#ngslib==1.1.0
#from ngslib import BigWigFile

# no step/span; requires numpy
import pyBigWig

# run bigwig tool externally
import subprocess

import ga4gh.server.protocol as protocol
import ga4gh.server.datamodel as datamodel
import ga4gh.server.datamodel.sequence_annotations as sequence_annotations
import ga4gh.server.exceptions as exceptions

import ga4gh.schemas.pb as pb



"""
COPY FROM SCHEMAS

Assuming 0-based for everything (convert wiggle). CHECK CODE AGAIN

Doesn't assume wiggle covers chromosome so that queries outside of 
range of file return 0 rather than an error (like pyBigWig).

THIS?
Could also return a "no data" message.
"""



class WiggleDataSource:
    """
    Class for reading from wiggle files. 
    """

    # is this in a good python style?
    _VARIABLE_STEP = 1
    _FIXED_STEP = 2

    def __init__(self, reference, start, end):
        self._mode = None
        self._span = 1
        self._step = 1
        self._data = protocol.Continuous()
        #self._data.start = start
        #self._pos = start 
        self._pos = None
        self._queryReference = reference
        self._queryStart = start
        self._queryEnd = end
   
    def getData(self):
        return self._data

    def parseStep(self, line):
        fields = dict( [ field.split( '=' ) for field in line.split()[1:] ] )

        if 'chrom' in fields:
            self._reference = fields['chrom']
        else:
            raise ValueError("Missing chrom field in %s" % line.strip()) 

        if line.startswith("fixedStep"):
            if 'start' in fields:
                self._start = int(fields['start']) - 1  # to 0-based
            else:
                raise ValueError("Missing start field in %s" % line.strip()) 

        if 'span' in fields:
            self._span = int(fields['span'])
        if 'step' in fields:
            self._step = int(fields['step'])

    def readWiggleLine(self, line):
        if(line.isspace() or line.startswith("#") 
                or line.startswith("browser") or line.startswith("track")):
            return 
        elif line.startswith("variableStep"):
            self._mode = self._VARIABLE_STEP
            self.parseStep(line)
            return 
        elif line.startswith("fixedStep"):
            self._mode = self._FIXED_STEP
            self.parseStep(line)
            return 
        elif self._mode is None:
            raise ValueError("Unexpected input line: %s" % line.strip())

        if self._queryReference != self._reference:
            return 

        # read data lines
        fields = line.split()
        if self._mode == self._VARIABLE_STEP:
            start = int(fields[0])-1  # to 0-based
            val = float(fields[1])
        else:
            start = self._start
            self._start += self._step
            val = float(fields[0])


        if start < self._queryEnd and start > self._queryStart:
            if self._pos is None:
                self._pos = start
                self._data.start = start

            # fill gap
            while self._pos<start:
                self._data.values.append(float('NaN'))
                self._pos += 1
            for _ in xrange(self._span):
                self._data.values.append(val)
            self._pos += self._span

        # for gapped data
        #data = protocol.Continuous()
        #data.start = start
        #data.values.append(val)
        #yield data

    def fillEnd(self):
        """
        Pad end values with NaN to fill query range
        """
        while self._pos<self._queryEnd:
                self._data.values.append(float('NaN'))
                self._pos += 1

    def wiggleFileHandleToProtocol(self, fileHandle):
        """
        Return a continuous protocol object satsifiying the given query
        parameters from the given wiggle file handle.
        """
        for line in fileHandle:
            self.readWiggleLine(line)
        #self.fillEnd()
        return self._data

    def wiggleFileToProtocol(self, fileName):
        """
        Return a continuous protocol object satsifiying the given query
        parameters from the given wiggle file
        """
        with open(fileName,'r') as f:
            return self.wiggleFileHandleToProtocol(f)

class BigWigDataSource:
    """
    Class for reading from bigwig files. 
    """

    def __init__(self, sourceFile):
        self._sourceFile = sourceFile

    def checkReference(self, reference):
        """
        Check the reference for security
        """
        pattern = re.compile(r'[\s,;"\'&\\]')
        if pattern.findall(reference.strip()):
            return False
        return True

    def readValuesPyBigWig(self, reference, start, end):
        """
        Use pyBigWig package to read a BigWig file for the
        given range and return a protocol object.
        """
        if not self.checkReference(reference):
            raise exceptions.ReferenceNameNotFoundException(reference)
        if start < 0 :
            raise exceptions.ReferenceRangeErrorException(
                reference, start, end)
        bw = pyBigWig.open(self._sourceFile)
        refLen = bw.chroms(reference)
        if refLen is None:
            raise exceptions.ReferenceNameNotFoundException(reference)
        if start > refLen or end > refLen:
            raise exceptions.ReferenceRangeErrorException(
                reference, start, end)

        data = protocol.Continuous()
        data.start = start
        startNan = True
        skipNanCount = 0
        for val in bw.values(reference, start, end):
            if math.isnan(val):
               skipNanCount += 1
            else:
                if startNan:
                    data.start += skipNanCount 
                    startNan = False
                    skipNanCount = 0
                for i in range(skipNanCount):
                    data.values.append(float('NaN'))
                skipNanCount = 0
                data.values.append(val)
        bw.close()
        return data

    def readValuesBigWigToWig(self, reference, start, end):
        """
        The a bigwig file and return a protocol object with values
        within the query range.

        This call uses the bigWigToWig command line from golden path.
        There could be memory issues if the returned results are large.

        The input reference can be a security problem. Ideally, it should be
        checked against a list of known chromosomes. Start and end should
        not be problems if they are integers.
        """
        if not self.checkReference(reference):
            raise exceptions.ReferenceNameNotFoundException(reference)
        if start < 0 :
            raise exceptions.ReferenceRangeErrorException(
                reference, start, end)
        # NEEDS TO CHECK IF QUERY IS BEYOND END

        cmd = ["bigWigToWig", self._sourceFile, "stdout", "-chrom="+reference,
                "-start="+str(start), "-end="+str(end)]
        wiggleReader = WiggleDataSource(reference, start, end)
        try:
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            while True:
                line = process.stdout.readline()
                if line == '' and process.poll() is not None:
                    break
                wiggleReader.readWiggleLine(line.strip())
        except ValueError:
            raise
        except:
            raise Exception("bigWigToWig failed to run")

        #wiggleReader.fillEnd()
        return wiggleReader.getData()
        
    def bigWigToProtocol(self, reference, start, end):
        #return self.readValuesBigWigToWig(reference, start, end)
        return self.readValuesPyBigWig(reference, start, end)


class SimulatedContinuousSet(sequence_annotations.AbstractFeatureSet):
    """
    Simulated data backend for ContinuousSet, used for internal testing.
    """
    def __init__(self, parentContainer, localId, randomSeed=1):
        self._randomSeed = randomSeed
        super(SimulatedContinuousSet, self).__init__(parentContainer, localId)

    def _generateSimulatedContinuous(self, randomNumberGenerator):
        continuous = protocol.Continuous()
        continuous.start = randomNumberGenerator.randint(1000, 2000)
        continuous.values = [100, 200.3, 400]

    def getContinuousData(self, referenceName=None, start=None, end=None):
        """
        Returns a set number of simulated continuous data.

        :param referenceName: name of reference to "search" on
        :param start: start coordinate of query
        :param end: end coordinate of query
        :return: Yields continuous list
        """
        randomNumberGenerator = random.Random()
        randomNumberGenerator.seed(self._randomSeed)
        for i in range(numResults):
            gaContinuous = self._generateSimulatedContinuous(
                                    randomNumberGenerator)
            match = (
                gaContinuous.start < end and
                gaContinuous.end > start and
                gaContinuous.reference_name == referenceName )
            if match:
                yield gaContinuous
