"""
Unit tests for continuous objects. This is used for all tests
that can be performed in isolation from input data.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import unittest

from nose.tools import raises

import ga4gh.server.datarepo as datarepo
import ga4gh.server.datamodel.continuous as continuous
import ga4gh.server.datamodel.datasets as datasets
import ga4gh.server.exceptions as exceptions

import tests.paths as paths


class TestContinuous(unittest.TestCase):
    """
    Unit tests for continuous data
    """
    def _createContinuousSet(self):
        """
        Creates a ContinuousSet from the specified directory.
        """
        self._continuousSetName = "testContinuous"
        self._repo = datarepo.SqlDataRepository(paths.testDataRepo)
        self._repo.open(datarepo.MODE_READ)
        self._dataset = datasets.Dataset("testDs")
        self._continuousSet = continuous.readSet(
            self._dataset, self._continuousSetName)

    def setUp(self):
        dataDir = "tests/data/datasets/dataset1/continuous"
        self._wiggleFile = dataDir + "/wiggle_2.txt"
        self._bigWigFile = dataDir + "/bigwig_1.bw"

    def testReadWiggle(self):
        continuousObj = continuous.WiggleReader(
                            'chr19', 49307698, 49308020)
        obj = continuousObj.wiggleFileToProtocol(self._wiggleFile)
        self.assertEqual(obj.start, 49307700)
        self.assertEqual(obj.values[0], 900)
        self.assertEqual(obj.values[300], 800)
        self.assertEqual(len(obj.values), 302)

    def testReadBigWig(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        obj = continuousObj.bigWigToProtocol("chr19", 49305897, 49306090)
        self.assertEqual(obj.start, 49305900)
        self.assertEqual(obj.values[0], 20.0)
        self.assertEqual(obj.values[183], 17.5)
        self.assertEqual(len(obj.values), 185)

    def testReadBigWigAllNan(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        obj = continuousObj.bigWigToProtocol("chr19", 49305927, 49305997)
        self.assertEqual(len(obj.values), 0)

    @raises(exceptions.ReferenceRangeErrorException)
    def testReadBigWigOutsideRange(self):
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        continuousObj.bigWigToProtocol("chr19", 493059030, 493059034)

    @raises(exceptions.ReferenceNameNotFoundException)
    def testReadBigWigChromsomeException(self):
        """
        Test for catching bad chromosome names.
        """
        continuousObj = continuous.BigWigDataSource(self._bigWigFile)
        continuousObj.bigWigToProtocol("chr&19", 49305602, 49308000)
