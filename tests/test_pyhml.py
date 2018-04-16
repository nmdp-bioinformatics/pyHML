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

from Bio import SeqIO
from pandas import DataFrame


class TestPyhml(unittest.TestCase):

    def setUp(self):
        self.data_dir = os.path.dirname(__file__) + "/resources"
        self.hmlparser = HmlParser(verbose=True)
        self.assertIsInstance(self.hmlparser, HmlParser)
        pass

    def test_001_gzip(self):
        hml_file = self.data_dir + "/3054.hml101.xml.gz"
        hml = self.hmlparser.parse(hml_file)
        hml_df = hml.toPandas()
        self.assertIsInstance(hml, HML)
        self.assertIsInstance(hml_df, DataFrame)
        hml_file1 = self.data_dir + "/3054.hml101-1.xml"
        hml_unzipped = self.data_dir + "/3054.hml101.xml"
        cmd_zip = "gzip " + hml_unzipped
        cmd_cp = "cp " + hml_file1 + " " + hml_unzipped
        os.system(cmd_cp)
        os.system(cmd_zip)
        pass

    def test_002_gzip(self):
        hml_file = self.data_dir + "/2609.hml101.xml.gz"
        hml = self.hmlparser.parse(hml_file)
        hml_df = hml.toPandas()
        self.assertIsInstance(hml, HML)
        self.assertIsInstance(hml_df, DataFrame)
        hml_file1 = self.data_dir + "/2609.hml101-1.xml"
        hml_unzipped = self.data_dir + "/2609.hml101.xml"
        cmd_zip = "gzip " + hml_unzipped
        cmd_cp = "cp " + hml_file1 + " " + hml_unzipped
        os.system(cmd_cp)
        os.system(cmd_zip)
        pass


