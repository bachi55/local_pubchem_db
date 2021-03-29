#####
# MIT License
#
# Copyright (c) 2017-2020 Eric Bach (eric.bach@aalto.fi)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#####

import os
import sqlite3
import unittest

from collections import OrderedDict

from pubchem2sqlite.utils import get_column_stmt, iter_sdf_file, extract_info_from_sdf, build_db


class TestParsingJSONDBSpecs(unittest.TestCase):
    def test_get_column_stmt(self):
        specs = OrderedDict([("MASS", {"DTYPE": "float", "NOT_NULL": False}),
                             ("INCHI", {"DTYPE": "string", "NOT_NULL": True}),
                             ("CID", {"DTYPE": "integer", "PRIMARY_KEY": True})])
        stmt = get_column_stmt(specs)
        self.assertEqual("MASS float,INCHI string not null,CID integer not null primary key", stmt)

        ######################################

        specs = OrderedDict([("MASS", {"DTYPE": "float", "NOT_NULL": False}),
                             ("INCHI", {"DTYPE": "string", "PRIMARY_KEY": True, "NOT_NULL": True}),
                             ("CID", {"DTYPE": "integer"})])
        stmt = get_column_stmt(specs)
        self.assertEqual("MASS float,INCHI string not null primary key,CID integer", stmt)

        ######################################

        specs = OrderedDict([("MASS", {"DTYPE": "float", "NOT_NULL": False}),
                             ("INCHI", {"DTYPE": "string", "PRIMARY_KEY": True}),
                             ("CID", {"DTYPE": "integer"})])
        stmt = get_column_stmt(specs)
        self.assertEqual("MASS float,INCHI string not null primary key,CID integer", stmt)

        ######################################

        # NOTE: The primary key label overwrites the not-null label. Primary keys are always not null.

        specs = OrderedDict([("MASS", {"DTYPE": "float", "NOT_NULL": False}),
                             ("INCHI", {"DTYPE": "string", "PRIMARY_KEY": True, "NOT_NULL": False}),
                             ("CID", {"DTYPE": "integer", "NOT_NULL": True})])
        stmt = get_column_stmt(specs)
        self.assertEqual("MASS float,INCHI string not null primary key,CID integer not null", stmt)


class TestSDFProcessing(unittest.TestCase):
    def setUp(self):
        self.base_dir = os.path.dirname(__file__)

    def test_sdf_molecule_iterator(self):
        cids = [31038, 31039, 31040]
        with open(os.path.join(self.base_dir, "sdf", "cmps_00_02.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                self.assertEqual(cids[idx], cid)

        cids = [34516, 34517, 34518]
        with open(os.path.join(self.base_dir, "sdf", "cmps_03_05.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                self.assertEqual(cids[idx], cid)

        cids = [46773, 46774]
        with open(os.path.join(self.base_dir, "sdf", "cmps_06_07.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                self.assertEqual(cids[idx], cid)

    def test_data_extraction(self):
        specs = {
            "columns": {
                "cid": {
                    "SD_TAG": ["PUBCHEM_COMPOUND_CID"],
                    "DTYPE": "integer",
                    "NOT_NULL": True,
                    "PRIMARY_KEY": True
                },
                "InChI": {
                    "SD_TAG": ["PUBCHEM_IUPAC_INCHI"],
                    "DTYPE": "varchar",
                    "NOT_NULL": True
                },
                "xlogp3": {
                    "SD_TAG": ["PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"],
                    "DTYPE": "real",
                    "NOT_NULL": False
                }
            }
        }

        inchis = [
            "InChI=1S/C18H31NO/c1-2-3-4-5-6-7-8-9-10-11-12-13-18-14-16-19(20)17-15-18/h14-17H,2-13H2,1H3",
            "InChI=1S/C11H18O2/c1-2-3-4-5-6-7-8-9-10-11(12)13/h1H,3-10H2,(H,12,13)",
            "InChI=1S/C5H6O5.2Na/c6-3(5(9)10)1-2-4(7)8;;/h1-2H2,(H,7,8)(H,9,10);;/q;2*+1/p-2"]
        xlogp3s = [6.6, 3.3, None]
        with open(os.path.join(self.base_dir, "sdf", "cmps_00_02.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                infos = extract_info_from_sdf(sdf, specs)

                self.assertEqual(cid, infos["cid"])
                self.assertEqual(inchis[idx], infos["InChI"])
                self.assertEqual(xlogp3s[idx], infos["xlogp3"])

        ############################################################################
        specs_no_AA_xlogp3 = specs
        specs_no_AA_xlogp3["columns"]["xlogp3"]["SD_TAG"] = ["PUBCHEM_XLOGP3"]

        inchis = [
            "InChI=1S/C18H31NO/c1-2-3-4-5-6-7-8-9-10-11-12-13-18-14-16-19(20)17-15-18/h14-17H,2-13H2,1H3",
            "InChI=1S/C11H18O2/c1-2-3-4-5-6-7-8-9-10-11(12)13/h1H,3-10H2,(H,12,13)",
            "InChI=1S/C5H6O5.2Na/c6-3(5(9)10)1-2-4(7)8;;/h1-2H2,(H,7,8)(H,9,10);;/q;2*+1/p-2"]
        xlogp3s = [None, 3.3, None]
        with open(os.path.join(self.base_dir, "sdf", "cmps_00_02.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                infos = extract_info_from_sdf(sdf, specs_no_AA_xlogp3)

                self.assertEqual(cid, infos["cid"])
                self.assertEqual(inchis[idx], infos["InChI"])
                self.assertEqual(xlogp3s[idx], infos["xlogp3"])

        ############################################################################
        specs_only_AA_xlogp3 = specs
        specs_only_AA_xlogp3["columns"]["xlogp3"]["SD_TAG"] = ["PUBCHEM_XLOGP3_AA"]

        inchis = [
            "InChI=1S/C18H31NO/c1-2-3-4-5-6-7-8-9-10-11-12-13-18-14-16-19(20)17-15-18/h14-17H,2-13H2,1H3",
            "InChI=1S/C11H18O2/c1-2-3-4-5-6-7-8-9-10-11(12)13/h1H,3-10H2,(H,12,13)",
            "InChI=1S/C5H6O5.2Na/c6-3(5(9)10)1-2-4(7)8;;/h1-2H2,(H,7,8)(H,9,10);;/q;2*+1/p-2"]
        xlogp3s = [6.6, None, None]
        with open(os.path.join(self.base_dir, "sdf", "cmps_00_02.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                infos = extract_info_from_sdf(sdf, specs_only_AA_xlogp3)

                self.assertEqual(cid, infos["cid"])
                self.assertEqual(inchis[idx], infos["InChI"])
                self.assertEqual(xlogp3s[idx], infos["xlogp3"])

    def test_data_transformation(self):
        specs = {
            "columns": {
                "cid": {
                    "SD_TAG": ["PUBCHEM_COMPOUND_CID"],
                    "DTYPE": "integer",
                    "NOT_NULL": True,
                    "PRIMARY_KEY": True,
                    "CREATE_LIKE": "lambda __x: 2 * __x"
                },
                "InChIKey": {
                    "SD_TAG": ["PUBCHEM_IUPAC_INCHIKEY"],
                    "DTYPE": "varchar",
                    "NOT_NULL": True
                },
                "InChIKey_1": {
                    "SD_TAG": ["PUBCHEM_IUPAC_INCHIKEY"],
                    "DTYPE": "varchar",
                    "NOT_NULL": True,
                    "CREATE_LIKE": "lambda __x: __x.split('-')[0]"
                },
                "xlogp3": {
                    "SD_TAG": ["PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"],
                    "DTYPE": "real",
                    "NOT_NULL": False,
                    "CREATE_LIKE": "lambda __x: round(__x)"
                }
            }
        }

        inchikeys = [
            "JGUZOCJCNMVJHU-UHFFFAOYSA-N",
            "OAOUTNMJEFWJPO-UHFFFAOYSA-N",
            "YBGBJYVHJTVUSL-UHFFFAOYSA-L"]
        xlogp3s = [6.6, 3.3, None]
        with open(os.path.join(self.base_dir, "sdf", "cmps_00_02.sdf"), "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                infos = extract_info_from_sdf(sdf, specs)

                self.assertEqual(2 * cid, infos["cid"])
                self.assertEqual(inchikeys[idx], infos["InChIKey"])
                self.assertEqual(inchikeys[idx].split("-")[0], infos["InChIKey_1"])

                if xlogp3s[idx] is None:
                    self.assertIsNone(infos["xlogp3"])
                else:
                    self.assertEqual(round(xlogp3s[idx]), infos["xlogp3"])


class TestDataImport(unittest.TestCase):
    def setUp(self):
        self.base_dir = os.path.dirname(__file__)
        os.makedirs(os.path.join(self.base_dir, "db"), exist_ok=True)
        self.db_fn = os.path.join(self.base_dir, "db", "pubchem.sqlite")

    def tearDown(self):
        """
        Close DB connection and remove DB file.
        """
        if self.conn:
            self.conn.close()

        if os.path.exists(self.db_fn):
            os.remove(self.db_fn)

    def test_db_import(self):
        specs = {
            "columns": {
                "cid": {
                    "SD_TAG": ["PUBCHEM_COMPOUND_CID"],
                    "DTYPE": "integer",
                    "NOT_NULL": True,
                    "PRIMARY_KEY": True
                },
                "inchikey": {
                    "SD_TAG": ["PUBCHEM_IUPAC_INCHIKEY"],
                    "DTYPE": "varchar",
                    "NOT_NULL": True
                },
                "InChI": {
                    "SD_TAG": ["PUBCHEM_IUPAC_INCHI"],
                    "DTYPE": "varchar",
                    "NOT_NULL": True
                },
                "xlogp3": {
                    "SD_TAG": ["PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"],
                    "DTYPE": "real",
                    "NOT_NULL": False
                }
            }
        }
        # Check that return value does not indicate any error
        self.assertEqual(0, build_db(self.base_dir, use_gzip=True, reset=True, db_specs=specs))

        self.conn = sqlite3.connect(self.db_fn)
        self.assertEqual(8,
                         self.conn.execute("SELECT count(*) FROM compounds").fetchall()[0][0])
        self.assertEqual("SISXGVIKZQKGLA-UHFFFAOYSA-N",
                         self.conn.execute("SELECT inchikey FROM compounds WHERE cid == 34516").fetchall()[0][0])
        self.assertEqual(6.6,
                         self.conn.execute("SELECT xlogp3 FROM compounds WHERE cid == 31038").fetchall()[0][0])
        self.assertEqual("InChI=1S/C5H6O5.2Na/c6-3(5(9)10)1-2-4(7)8;;/h1-2H2,(H,7,8)(H,9,10);;/q;2*+1/p-2",
                         self.conn.execute("SELECT InChI FROM compounds WHERE cid == 31040").fetchall()[0][0])

        ################################################################

        # NOTE: Each compound for which at least one "NOT_NULL" information is None (=NULL) will be not added to the DB.

        specs_not_null_xlogp = specs
        specs_not_null_xlogp["columns"]["xlogp3"]["NOT_NULL"] = True
        # Check that return value does not indicate any error
        self.assertEqual(0, build_db(self.base_dir, use_gzip=True, reset=True, db_specs=specs))

        self.conn = sqlite3.connect(self.db_fn)
        self.assertEqual(5,
                         self.conn.execute("SELECT count(*) FROM compounds").fetchall()[0][0])
        cids = [r[0] for r in self.conn.execute("SELECT cid FROM compounds")]
        self.assertNotIn(34516, cids)
        self.assertNotIn(31040, cids)
        self.assertNotIn(46774, cids)


if __name__ == '__main__':
    unittest.main()
