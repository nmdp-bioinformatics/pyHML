#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    pyhml pyHML.
#    Copyright (c) 2017 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#



"""
test_pyhml
----------------------------------

Tests for `pyhml` module.
"""

import os
import sys
import unittest

from pyhml.pyhml import HmlParser
from pyhml.models.hml import HML
from pyhml.models.sample import Sample
from pyhml.models.typing import Typing
from pyhml.models.haploid import Haploid
from pyhml.models.consensus import Consensus
from pyhml.models.ref_database import RefDatabase
from pyhml.models.ref_sequence import RefSequence
from pyhml.models.reporting_center import ReportingCenter
from pyhml.models.allele_assignment import AlleleAssignment
from pyhml.models.consensus_seq_block import ConsensusSeqBlock

from Bio import SeqIO
from pandas import DataFrame


class TestPyhml(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        pass

    def test_001_gzip(self):
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        hml_file = self.data_dir + "/1111.hml101.xml.gz"
        hml = self.hmlparser.parse(hml_file)
        hml_df = hml.toPandas()
        self.assertIsInstance(hml, HML)
        self.assertIsInstance(hml_df, DataFrame)
        hml_file1 = self.data_dir + "/1111.hml101-1.xml"
        hml_unzipped = self.data_dir + "/1111.hml101.xml"
        cmd_zip = "gzip " + hml_unzipped
        cmd_cp = "cp " + hml_file1 + " " + hml_unzipped
        os.system(cmd_cp)
        os.system(cmd_zip)
        pass

    def test_002_gzip(self):
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        hml_file = self.data_dir + "/2222.hml101.xml.gz"
        hml = self.hmlparser.parse(hml_file)
        hml_df = hml.toPandas()
        hml.tobiotype(self.data_dir)
        fasta_file = self.data_dir + "/ABDR.fasta"
        seqs = list(SeqIO.parse(fasta_file, "fasta"))
        self.assertEqual(len(seqs), 20)
        self.assertIsInstance(hml, HML)
        self.assertIsInstance(hml_df, DataFrame)
        hml_file1 = self.data_dir + "/2222.hml101-1.xml"
        hml_unzipped = self.data_dir + "/2222.hml101.xml"
        cmd_rm = "rm " + fasta_file
        cmd_zip = "gzip " + hml_unzipped
        cmd_cp = "cp " + hml_file1 + " " + hml_unzipped
        os.system(cmd_rm)
        os.system(cmd_cp)
        os.system(cmd_zip)
        pass

    def test_003_sample(self):
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        sample = Sample(center_code='000', id='100000',
                        collection_method='DNA')
        self.assertIsInstance(sample, Sample)
        self.assertEqual(sample.id, '100000')
        self.assertEqual(sample.center_code, '000')
        self.assertEqual(sample.collection_method, 'DNA')
        pass

    def test_004_sampleTyping(self):
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        sample = Sample(center_code='000', id='100000',
                        collection_method='DNA')
        typing = Typing(date='04-11-2011', gene_family='HLA')
        sample.typing = [typing]
        self.assertIsInstance(typing, Typing)
        self.assertIsInstance(sample, Sample)
        self.assertEqual(typing.gene_family, 'HLA')
        self.assertIsInstance(sample.typing[0], Typing)
        self.assertEqual(sample.id, '100000')
        self.assertEqual(sample.center_code, '000')
        self.assertEqual(sample.collection_method, 'DNA')
        pass

    def test_005_haploid(self):
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        haploid = Haploid(locus="HLA-A",
                          method="DNA",
                          type="HLA-A*01:01")
        allele_assignment = AlleleAssignment(allele_db="3.31.0",
                                             allele_version="3.31.0",
                                             date="04-11-2011",
                                             haploid=[haploid])
        self.assertIsInstance(haploid, Haploid)
        self.assertEqual(haploid.method, 'DNA')
        self.assertEqual(haploid.type, "HLA-A*01:01")
        self.assertEqual(allele_assignment.allele_db, "3.31.0")
        self.assertIsInstance(allele_assignment, AlleleAssignment)
        self.assertEqual(allele_assignment.allele_version, "3.31.0")
        self.assertIsInstance(allele_assignment.haploid[0], Haploid)
        pass

    def test_006_hmlversion(self):
        hmlparser = HmlParser(verbose=True, hmlversion='1.0.1')
        self.assertIsInstance(hmlparser, HmlParser)
        self.assertEqual(hmlparser.hmlversion, '1.0.1')
        self.assertEqual(len(hmlparser.schemas), 1)

    def test_007_hml(self):
        hml = HML(project_name="HML TEST", version='1.0.1')
        self.assertIsInstance(hml, HML)
        self.assertEqual(hml.project_name, "HML TEST")
        self.assertEqual(hml.version, '1.0.1')
