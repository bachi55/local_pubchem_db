import unittest

from misc import get_column_stmt, iter_sdf_file, extract_info_from_sdf


class TestParsingJSONDBSpecs(unittest.TestCase):
    def test_get_column_stmt(self):
        specs = {
            "columns": {
                "MASS": {
                    "DTYPE": "float", "NOT_NULL": False
                },
                "INCHI": {
                    "DTYPE": "string", "NOT_NULL": True
                },
                "CID": {
                    "DTYPE": "integer", "PRIMARY_KEY": True
                }
            }
        }

        stmt = get_column_stmt(specs["columns"])
        self.assertEqual("MASS float,INCHI string not null,CID integer not null primary key", stmt)


class TestSDFProcessing(unittest.TestCase):
    def test_sdf_molecule_iterator(self):
        cids = [31038, 31039, 31040]
        with open("cmps_00_02.sdf", "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                self.assertEqual(cids[idx], cid)

        cids = [34516, 34517, 34518]
        with open("cmps_03_05.sdf", "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                self.assertEqual(cids[idx], cid)

        cids = [46773, 46774]
        with open("cmps_06_07.sdf", "r") as sdf_file:
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
        with open("cmps_00_02.sdf", "r") as sdf_file:
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
        with open("cmps_00_02.sdf", "r") as sdf_file:
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
        with open("cmps_00_02.sdf", "r") as sdf_file:
            for idx, (cid, sdf) in enumerate(iter_sdf_file(sdf_file)):
                infos = extract_info_from_sdf(sdf, specs_only_AA_xlogp3)

                self.assertEqual(cid, infos["cid"])
                self.assertEqual(inchis[idx], infos["InChI"])
                self.assertEqual(xlogp3s[idx], infos["xlogp3"])


class TestDataImport(unittest.TestCase):
    def test_db_import(self):
        pass


if __name__ == '__main__':
    unittest.main()
