#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################
##  DendroPy Phylogenetic Computing Library.
##
##  Copyright 2010-2015 Jeet Sukumaran and Mark T. Holder.
##  All rights reserved.
##
##  See "LICENSE.rst" for terms and conditions of usage.
##
##  If you use this work or any portion thereof in published work,
##  please cite it as:
##
##     Sukumaran, J. and M. T. Holder. 2010. DendroPy: a Python library
##     for phylogenetic computing. Bioinformatics 26: 1569-1571.
##
##############################################################################

"""
Regression tests for ``dendropy.interop.genbank`` annotation export.

These use small, hand-built INSDSeq XML fragments (no network access) so
that they exercise ``GenBankAccessionFeature`` / ``GenBankAccessionRecord``
parsing and ``as_annotation()`` in isolation.
"""

import unittest
from xml.etree import ElementTree

from dendropy.interop import genbank


def _build_feature_xml():
    xml_str = """
    <INSDFeature>
      <INSDFeature_key>source</INSDFeature_key>
      <INSDFeature_location>1..100</INSDFeature_location>
      <INSDFeature_intervals>
        <INSDInterval>
          <INSDInterval_from>1</INSDInterval_from>
          <INSDInterval_to>100</INSDInterval_to>
          <INSDInterval_accession>TEST00001.1</INSDInterval_accession>
        </INSDInterval>
      </INSDFeature_intervals>
    </INSDFeature>
    """
    return ElementTree.fromstring(xml_str)


def _build_record_xml():
    # Deliberately omits several optional fields (e.g. INSDSeq_topology,
    # INSDSeq_strandedness, INSDSeq_division, dates, taxonomy) so that the
    # corresponding attributes on GenBankAccessionRecord remain None, as
    # is common for real (partial) GenBank records.
    xml_str = """
    <INSDSeq>
      <INSDSeq_locus>TESTLOCUS</INSDSeq_locus>
      <INSDSeq_length>100</INSDSeq_length>
      <INSDSeq_moltype>DNA</INSDSeq_moltype>
      <INSDSeq_definition>Test definition</INSDSeq_definition>
      <INSDSeq_primary-accession>TEST00001</INSDSeq_primary-accession>
      <INSDSeq_accession-version>TEST00001.1</INSDSeq_accession-version>
      <INSDSeq_other-seqids>
        <INSDSeqid>gi|12345</INSDSeqid>
      </INSDSeq_other-seqids>
      <INSDSeq_source>Test source organism</INSDSeq_source>
      <INSDSeq_organism>Testus organismus</INSDSeq_organism>
      <INSDSeq_feature-table>
        <INSDFeature>
          <INSDFeature_key>source</INSDFeature_key>
          <INSDFeature_location>1..100</INSDFeature_location>
          <INSDFeature_intervals>
            <INSDInterval>
              <INSDInterval_from>1</INSDInterval_from>
              <INSDInterval_to>100</INSDInterval_to>
              <INSDInterval_accession>TEST00001.1</INSDInterval_accession>
            </INSDInterval>
          </INSDFeature_intervals>
        </INSDFeature>
      </INSDSeq_feature-table>
    </INSDSeq>
    """
    return ElementTree.fromstring(xml_str)


class TestGenBankAccessionFeatureAsAnnotation(unittest.TestCase):

    def test_interval_parsing_scoped_to_intervals_element(self):
        # Regression test: parse_xml() used to search for "INSDInterval"
        # starting at the top-level feature element instead of the
        # "INSDFeature_intervals" sub-element, so intervals were silently
        # never parsed (ElementTree.findall() with a bare tag name only
        # matches direct children).
        feature = genbank.GenBankAccessionFeature(_build_feature_xml())
        self.assertEqual(len(feature.intervals), 1)
        interval = feature.intervals[0]
        self.assertEqual(interval.begin, "1")
        self.assertEqual(interval.end, "100")
        self.assertEqual(interval.accession, "TEST00001.1")

    def test_as_annotation_does_not_crash_and_uses_interval_values(self):
        feature = genbank.GenBankAccessionFeature(_build_feature_xml())
        # Regression test: as_annotation() used to call
        # getattr(self, item[0]) instead of getattr(interval, item[0]),
        # which raised AttributeError since GenBankAccessionFeature has no
        # "begin"/"end"/"accession" attributes (those belong to the
        # per-interval GenBankAccessionInterval objects).
        top = feature.as_annotation()
        intervals_annote = [
            a for a in top.annotations if a.name == "INSDSeq_intervals"
        ][0]
        interval_annote = [
            a for a in intervals_annote.annotations if a.name == "INSDInterval"
        ][0]
        sub_by_name = {a.name: a.value for a in interval_annote.annotations}
        self.assertEqual(sub_by_name.get("INSDInterval_from"), "1")
        self.assertEqual(sub_by_name.get("INSDInterval_to"), "100")
        self.assertEqual(sub_by_name.get("INSDInterval_accession"), "TEST00001.1")
        # Regression test: as_annotation() used to do
        # "interval_annote.annotations.add(interval_annote)" (adding the
        # annotation to itself) instead of "interval_annote.annotations.add(sub)".
        self.assertNotIn(interval_annote, list(interval_annote.annotations))


class TestGenBankAccessionRecordAsAnnotation(unittest.TestCase):

    def test_as_annotation_skips_none_fields_without_crash(self):
        record = genbank.GenBankAccessionRecord(db="nucleotide", xml=_build_record_xml())
        # sanity check: these optional fields really are absent/None
        self.assertIsNone(record.topology)
        self.assertIsNone(record.strandedness)
        self.assertIsNone(record.division)
        # Regression test: as_annotation() used to leave the local
        # variable "a" unbound whenever a field's value was None (neither
        # the hasattr(value, "as_annotation") nor the "value is not None"
        # branch fired), raising UnboundLocalError on top.annotations.add(a).
        top = record.as_annotation()
        names = [a.name for a in top.annotations]
        self.assertNotIn("INSDSeq_topology", names)
        self.assertNotIn("INSDSeq_strandedness", names)
        self.assertNotIn("INSDSeq_division", names)
        self.assertIn("INSDSeq_locus", names)
        self.assertIn("INSDSeq_organism", names)


if __name__ == "__main__":
    unittest.main()
