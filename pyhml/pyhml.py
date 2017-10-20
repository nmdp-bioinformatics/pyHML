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
import os
import re
import xmlschema
import xmltodict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
from Bio.Alphabet import IUPAC

from pyhml.models.hml import HML
from pyhml.models.reporting_center import ReportingCenter
from pyhml.models.sample import Sample
from pyhml.models.typing import Typing
from pyhml.models.allele_assignment import AlleleAssignment
from pyhml.models.consensus import Consensus
from pyhml.models.typing_method import TypingMethod
from pyhml.models.consensus_seq_block import ConsensusSeqBlock
from pyhml.models.ref_database import RefDatabase
from pyhml.models.ref_sequence import RefSequence
import pandas as pd


class HmlParser(object):
    """
    hmlparser = HmlParser()
    hml = hmlparser.parse(hml)
    print(summary(hml, by='subject'))

    printHml(hml, by='subject', type='fasta')

    hml_filters = hml.filterby('date', cmp='lt', value='')
    hml_filters = hml.filterby('locus', value='HLA-A')
    df = toDF(hml)

    sequences = hml_filters.sequence_records
    SeqIO.write(sequence, 'seqfile.fasta', 'fasta')
    SeqIO.write(sequence, 'seqfile.imgt', 'imgt')

    """
    def __init__(self, hmlversion=None):
        """
        HmlParser - a model
        """
        data_dir = os.path.dirname(__file__)
        self.schemas = {}
        self.versions = ['1.0.1', '1.0', '0.9.4', '0.9.5', '0.9.6', '0.9.7', '1.0.2']
        if not hmlversion:
            for ver in self.versions:
                xsd_file = data_dir + '/data/hml-' + ver + '.xsd'
                self.schemas.update({ver: xmlschema.XMLSchema(xsd_file)})

    def parse(self, hml_file):
        hml_version = self._get_version(hml_file)
        schema = self.schemas[hml_version]
        hml_data = schema.to_dict(hml_file)
        rpc = ReportingCenter(reporting_center_context=hml_data['hmlns:reporting-center']['@reporting-center-context'],
                              reporting_center_id=hml_data['hmlns:reporting-center']['@reporting-center-id'])
        hml = HML(project_name=hml_data['@project-name'],
                  version=hml_data['@version'],
                  schema_location=hml_data['@xsi:schemaLocation'],
                  reporting_center=rpc)
        samples = []
        for i in range(0, len(hml_data['hmlns:sample'])):
            sample_id = hml_data['hmlns:sample'][i]['@id']
            center_code = hml_data['hmlns:sample'][i]['@center-code']
            collection_method = hml_data['hmlns:sample'][i]['hmlns:collection-method']
            sample = Sample(center_code=center_code, id=sample_id,
                            collection_method=collection_method)
            typings = []
            for typing_data in hml_data['hmlns:sample'][i]['hmlns:typing']:
                typing_date = typing_data['@date']
                gene_family = typing_data['@gene-family']

                typing = Typing(date=typing_date, gene_family=gene_family)

                allele_assignments = []
                for assignment in typing_data['hmlns:allele-assignment']:
                    allele_db = assignment['@allele-db']
                    db = assignment['@allele-version']
                    type_date = assignment['@date']
                    gls = [gl.strip() for gl in assignment['hmlns:glstring']
                           if re.search("\*\d", gl)]
                    allele_assignment = AlleleAssignment(allele_db=allele_db,
                                                         allele_version=db,
                                                         date=type_date,
                                                         glstring=gls)
                    allele_assignments.append(allele_assignment)
                consensus_seqs = []
                if 'hmlns:consensus-sequence' in typing_data:
                    for consensus in typing_data['hmlns:consensus-sequence']:
                        blocks = []
                        for cbd in consensus['hmlns:consensus-sequence-block']:
                            consensus_seq = ''.join([c.strip() for c
                                                     in cbd['hmlns:sequence']
                                                     if re.search("\D", c)])
                            seq = Seq(consensus_seq, IUPAC.unambiguous_dna)
                            con_b = ConsensusSeqBlock(continuity=cbd['@continuity'],
                                                      description=cbd['@description'],
                                                      end=cbd['@end'],
                                                      expected_copy_number=cbd['@expected-copy-number'],
                                                      phase_set=cbd['@phase-set'],
                                                      reference_sequence_id=cbd['@reference-sequence-id'],
                                                      start=cbd['@start'],
                                                      strand=str(cbd['@strand']),
                                                      sequence=seq)
                            blocks.append(con_b)
                        ref_dbs = []
                        for ref_data in consensus['hmlns:reference-database']:
                            refseqs = []
                            for seq_data in ref_data['hmlns:reference-sequence']:
                                ref_seq = RefSequence(accession=seq_data['@accession'],
                                                      end=seq_data['@end'],
                                                      id=seq_data['@id'],
                                                      name=seq_data['@name'],
                                                      start=seq_data['@start'],
                                                      uri=seq_data['@uri'])
                                refseqs.append(ref_seq)
                            ref_db = RefDatabase(availability=ref_data['@availability'],
                                                 curated=ref_data['@curated'],
                                                 description=ref_data['@description'],
                                                 name=ref_data['@name'],
                                                 uri=ref_data['@uri'],
                                                 version=ref_data['@version'],
                                                 reference_sequence=refseqs)
                            ref_dbs.append(ref_db)
                        cons = Consensus(date=consensus['@date'],
                                         consensus_sequence_block=blocks,
                                         reference_database=ref_dbs)
                        consensus_seqs.append(cons)
                    typing.allele_assignment = allele_assignments
                    typing.consensus_sequence = consensus_seqs
                    typings.append(typing)
            sample.typing = typings
            sample.create_seqrecords()
            samples.append(sample)
        hml.sample = samples
        return hml

    def _get_version(self, hmlfile):
        """
        Sets the typing of this Sample.

        :param typing: The typing of this Sample.
        :type typing: List[Typing]
        """
        doc = ''
        with open(hmlfile) as fd:
            doc = xmltodict.parse(fd.read())
            fd.close()
        return doc['hml']['@version']


def tobiotype(hml, outdir, dtype='fasta', by='subject'):
    """
    Sets the typing of this Sample.

    :param typing: The typing of this Sample.
    :type typing: List[Typing]
    """
    if by == 'subject':
        for sample in hml.sample:
            sequence_records = []
            for loc in sample.seq_records:
                for seqrc in sample.seq_records[loc]:
                    sequence_records.append(seqrc)

            outfile = outdir + '/' + sample.id + '.' + dtype
            SeqIO.write(sequence_records, outfile, dtype)
    else:
        sequence_records = []
        for sample in hml.sample:
            for loc in sample.seq_records:
                for seqrc in sample.seq_records[loc]:
                    sequence_records.append(seqrc)

        outfile = outdir + '/' + hml.project_name + '.' + dtype
        SeqIO.write(sequence_records, outfile, dtype)


def toDF(hml):
    """
    Sets the typing of this Sample.

    :param typing: The typing of this Sample.
    :type typing: List[Typing]
    """
    data = []
    for sample in hml.sample:
        for loc in sample.seq_records:
            for seqrc in sample.seq_records[loc]:
                seq = str(seqrc.seq)
                gl = seqrc.name
                db = seqrc.description
                row = [sample.id, loc, gl, db, seq]
                data.append(row)
    return pd.DataFrame(data, columns=['ID', 'Locus', 'glstring', 'dbversion', 'sequence'])









