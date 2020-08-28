import glob
import sqlite3
import os
import sys
import gzip
import traceback
import argparse
import json

from datetime import datetime
from timeit import default_timer as timer

from misc import split_sdf_file, initialize_db, load_db_specifications

"""
This script can be used to parse a set of SDF files and add (some of) the
content to a SQLite database.

This code is based on Huibin's implementation.
"""


def get_sdf_files_not_in_db(db_connection, sdf_fn_in_folder):
    """
    Returns the sdf_file names (full path) which are not already in the DB.
    
    :param db_connection: sqlite3.Connection, database connection
    :param sdf_fn_in_folder: list of string, filenames of the sdf-files.
    """
    sdf_files_in_db = db_connection.execute("SELECT filename FROM sdf_file").fetchall()
    sdf_files_in_db = [str(x[0]) for x in sdf_files_in_db]
                             
    return list(filter(lambda x: os.path.basename(x) not in sdf_files_in_db, sdf_fn_in_folder))


def insert_info_from_sdf_strings(db_connection, l_cid_sdf):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table.

    :param db_connection: sqlite3.Connection, database connection
    :param l_cid_sdf: list of (cid, sdf-string)-tuples
    """
    start = timer()

    db_connection.execute("BEGIN")

    for cid, sdf in l_cid_sdf:
        mol_sdf = extract_info_from_sdf(sdf)

        assert (cid == mol_sdf["PUBCHEM_COMPOUND_CID"])
        if mol_sdf["PUBCHEM_COMPOUND_CID"] is None or mol_sdf["PUBCHEM_IUPAC_INCHI"] is None or \
                mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"] is None or mol_sdf["PUBCHEM_MOLECULAR_FORMULA"] is None:
            continue

        # TODO: Use execute many here. This should remove the need to explicitly have a begin statement here
        db_connection.execute("INSERT INTO info "
                              "(cid, iupac_name, iupac_inchi, iupac_inchikey, inchikey1, xlogp3,"
                              "monoisotopic_mass, molecular_formula, smiles, can_smiles, atom_def_stereo_count, "
                              "atom_udef_stereo_count, bond_def_stereo_count, bond_udef_stereo_count) "
                              "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
            cid,
            mol_sdf["PUBCHEM_IUPAC_NAME"],
            mol_sdf["PUBCHEM_IUPAC_INCHI"],
            mol_sdf["PUBCHEM_IUPAC_INCHIKEY"],
            mol_sdf["INCHIKEY1"],
            mol_sdf["PUBCHEM_XLOGP3*"],
            mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"],
            mol_sdf["PUBCHEM_MOLECULAR_FORMULA"],
            mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
            mol_sdf["PUBCHEM_OPENEYE_CAN_SMILES"],
            mol_sdf["PUBCHEM_ATOM_DEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_ATOM_UDEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_BOND_DEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_BOND_UDEF_STEREO_COUNT"])
        )

    db_connection.execute("COMMIT")

    end = timer()
    print("Extraction and insertion of the information took %.3fsec" % (end - start))


def extract_info_from_sdf(sdf):
    """
    Extract the information for a given molecules from its sdf-string. 
    
    :param sdf: string, containing the molecule's sdf
    
    return: dictionary, containing the extracted information.
    """
    
    info = {
        "PUBCHEM_COMPOUND_CID":           None,  # integer
        "PUBCHEM_IUPAC_NAME":             None,  # string
        "PUBCHEM_IUPAC_INCHI":            None,  # string
        "PUBCHEM_IUPAC_INCHIKEY":         None,  # string
        "INCHIKEY1":                      None,  # string
        "PUBCHEM_XLOGP3*":                None,  # float
        "PUBCHEM_MONOISOTOPIC_WEIGHT":    None,  # float
        "PUBCHEM_MOLECULAR_FORMULA":      None,  # string
        "PUBCHEM_OPENEYE_ISO_SMILES":     None,  # string
        "PUBCHEM_OPENEYE_CAN_SMILES":     None,  # string
        "PUBCHEM_ATOM_DEF_STEREO_COUNT":  None,  # integer
        "PUBCHEM_ATOM_UDEF_STEREO_COUNT": None,  # integer
        "PUBCHEM_BOND_DEF_STEREO_COUNT":  None,  # integer
        "PUBCHEM_BOND_UDEF_STEREO_COUNT": None   # integer
    }
    found = {key: False for key in info.keys()}

    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0

    while not all(found.values()) and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i = i + 1
            continue

        if not found["PUBCHEM_COMPOUND_CID"] and lines[i].find("PUBCHEM_COMPOUND_CID") != -1:
            found["PUBCHEM_COMPOUND_CID"] = True
            info["PUBCHEM_COMPOUND_CID"] = int(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_NAME"] and lines[i].find("PUBCHEM_IUPAC_NAME") != -1:
            found["PUBCHEM_IUPAC_NAME"] = True
            info["PUBCHEM_IUPAC_NAME"] = lines[i+1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_INCHI"] and lines[i].find("PUBCHEM_IUPAC_INCHI") != -1:
            found["PUBCHEM_IUPAC_INCHI"] = True
            info["PUBCHEM_IUPAC_INCHI"] = lines[i+1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_INCHIKEY"] and lines[i].find("PUBCHEM_IUPAC_INCHIKEY") != -1:
            found["PUBCHEM_IUPAC_INCHIKEY"] = True
            found["INCHIKEY1"] = True
            info["PUBCHEM_IUPAC_INCHIKEY"] = lines[i+1].strip()
            info["INCHIKEY1"] = info["PUBCHEM_IUPAC_INCHIKEY"].split("-")[0]
            i = i + 2
            continue

        # PubChem uses XLOGP3 and XLOGP3_AA without any distinction [1]. Therefore, here we search for both information
        # PUBCHEM_XLOGP3 and PUBCHEM_XLOGP3_AA. We take PUBCHEM_XLOGP3 if it was found earlier.
        #
        #  - XLOGP3: 'standard' model [2]
        #  - XLOGP3_AA: predicted logp using pure atom-additive model [2]
        #
        # [1]: http://www.sioc-ccbg.ac.cn/skins/ccbgwebsite/software/xlogp3/manual/XLOGP3_Manual.pdf
        # [2]: ftp://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf (page 9)
        if not found["PUBCHEM_XLOGP3*"] and lines[i].find("PUBCHEM_XLOGP3") != -1:
            found["PUBCHEM_XLOGP3*"] = True
            info["PUBCHEM_XLOGP3*"] = float(lines[i+1].strip())
            i = i + 2
            continue
        if not found["PUBCHEM_XLOGP3*"] and lines[i].find("PUBCHEM_XLOGP3_AA") != -1:
            found["PUBCHEM_XLOGP3*"] = True
            info["PUBCHEM_XLOGP3*"] = float(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_MONOISOTOPIC_WEIGHT"] and lines[i].find("PUBCHEM_MONOISOTOPIC_WEIGHT") != -1:
            found["PUBCHEM_MONOISOTOPIC_WEIGHT"] = True
            info["PUBCHEM_MONOISOTOPIC_WEIGHT"] = float(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_MOLECULAR_FORMULA"] and lines[i].find("PUBCHEM_MOLECULAR_FORMULA") != -1:
            found["PUBCHEM_MOLECULAR_FORMULA"] = True
            info["PUBCHEM_MOLECULAR_FORMULA"] = lines[i+1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_OPENEYE_ISO_SMILES"] and lines[i].find("PUBCHEM_OPENEYE_ISO_SMILES") != -1:
            found["PUBCHEM_OPENEYE_ISO_SMILES"] = True
            info["PUBCHEM_OPENEYE_ISO_SMILES"] = lines[i+1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_OPENEYE_CAN_SMILES"] and lines[i].find("PUBCHEM_OPENEYE_CAN_SMILES") != -1:
            found["PUBCHEM_OPENEYE_CAN_SMILES"] = True
            info["PUBCHEM_OPENEYE_CAN_SMILES"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_ATOM_DEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_ATOM_DEF_STEREO_COUNT") != -1:
            found["PUBCHEM_ATOM_DEF_STEREO_COUNT"] = True
            info["PUBCHEM_ATOM_DEF_STEREO_COUNT"] = int(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_ATOM_UDEF_STEREO_COUNT") != -1:
            found["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] = True
            info["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] = int(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_BOND_DEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_BOND_DEF_STEREO_COUNT") != -1:
            found["PUBCHEM_BOND_DEF_STEREO_COUNT"] = True
            info["PUBCHEM_BOND_DEF_STEREO_COUNT"] = int(lines[i+1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_BOND_UDEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_BOND_UDEF_STEREO_COUNT") != -1:
            found["PUBCHEM_BOND_UDEF_STEREO_COUNT"] = True
            info["PUBCHEM_BOND_UDEF_STEREO_COUNT"] = int(lines[i+1].strip())
            i = i + 2
            continue

        i = i + 1

    return info


def main():
    arg_parser = argparse.ArgumentParser()

    # Required arguments:
    arg_parser.add_argument("base_dir", type=str,
                            help="Base-directory containing containing the 'db/' and 'sdf/' folders.")

    # Optional arguments
    arg_parser.add_argument("--gzip", action="store_true",
                            help="If true, sdf-files are assumed to be compressed using gzip and do have file extension"
                                 "'.gz'.")
    arg_parser.add_argument("--reset", action="store_true",
                            help="If true, all existing tables will be deleted and the DB will be re-build.")
    arg_parser.add_argument("--db_layout_fn", type=str, default="./default_db_layout.json",
                            help="JSON-file specifying the database layout.")

    # Parse arguments:
    args = arg_parser.parse_args()

    # Which SDF reader to use
    if args.gzip:
        sdf_opener = lambda fn: gzip.open(fn, "rt")
    else:
        sdf_opener = lambda fn: open(fn, "r")

    # Load the DB layout
    db_specs = load_db_specifications(args.db_layout_fn)

    # Directory paths to the SDF and output DB file
    sdf_dir = os.path.join(args.base_dir, "sdf")
    db_dir = os.path.join(args.base_dir, "db")

    # Connect to the 'pubchem' database.
    conn = sqlite3.connect(os.path.join(db_dir, "pubchem_" + str(datetime.now().date()) + ".db"))

    try:
        # Initialize the database
        with conn:
            initialize_db(conn, db_specs, reset=args.reset)

        # Get a list of all sdf-files available and reduce it to the ones still
        # needed to be processed.
        fn_patter = "*.sdf.gz" if args.gzip else "*.sdf"
        sdf_files = glob.glob(os.path.join(sdf_dir, fn_patter))
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (before filtering): %d" % n_sdf_files)

        with conn:
            sdf_files = get_sdf_files_not_in_db(conn, sdf_files)
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (after filtering): %d" % n_sdf_files)

        if n_sdf_files > 0:
            # Iterate over the sdf-files and add them one by one
            for ii, sdf_fn in enumerate(sdf_files):
                print("Process sdf-file: %s (%d/%d)" % (os.path.basename(sdf_fn), ii + 1, n_sdf_files))

                # parse and insert current sdf-file
                with sdf_opener(sdf_fn) as sdf_file, conn:
                    insert_info_from_sdf_strings(conn, split_sdf_file(sdf_file))

                    # add current sdf-file to the list of completed sdf-files
                    conn.execute("INSERT INTO sdf_file (filename, lowest_cid, highest_cid) VALUES(?,?,?)", (
                        os.path.basename(sdf_fn),
                        os.path.basename(sdf_fn).split(".")[0].split("_")[1],
                        os.path.basename(sdf_fn).split(".")[0].split("_")[2]))

        # Make the column 'monoisotopic_mass' an index column. Later, when we query by
        # monoisotopic mass with ppm-window we _hugely_ speed up the query.
        with conn:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_monoisotopic_mass ON info(monoisotopic_mass)")
        print("Create index on monoisotopic mass.")

        # Create an index on the inchikey1 to query stereo-isomers from the database given
        # their 2D structure only.
        with conn:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_inchikey1 ON info(inchikey1)")
        print("Create index on inchikey1.")

        # Create an index on the inchi's
        with conn:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_inchi ON info(iupac_inchi)")
        print("Create index on inchi.")

        # Create an index on the molecular formula
        with conn:
            conn.execute("CREATE INDEX IF NOT EXISTS idx_molecular_formula ON info(molecular_formula)")
        print("Create index on molecular formula.")

        return_code = 0

    except sqlite3.ProgrammingError as err:
        print("Programming error: '" + err.args[0] + "'.")
        traceback.print_exc()
        return_code = 1
    except sqlite3.DatabaseError as err:
        print("Database error: '" + err.args[0] + "'.")
        traceback.print_exc()
        return_code = 1
    except IOError as err:
        print("An IOError occurred: '" + os.strerror(err.args[0]) + "'.")
        traceback.print_exc()
        return_code = 1
    except Exception as err:
        print(err.args[0])
        traceback.print_exc()
        return_code = 1

    finally:
        conn.close()

    return return_code


if __name__ == "__main__":
    sys.exit(main())
