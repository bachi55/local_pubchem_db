import glob
import sqlite3
import os
import sys
from timeit import default_timer as timer

from misc import split_sdf_file

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


def initialize_db(db_connection, reset=False):
    """
    Initialization of the DB: 
        - Remove existing tables if desired.
        - Create not present tables. 
        
    :param db_connection: sqlite3.Connection, connection to the database to modify.
    :param reset: boolean, If 'True' tables in the database are dropped and recreated.
    """
    if not reset:
        return

    db_connection.execute("DROP TABLE IF EXISTS sdf_file")
    db_connection.execute("CREATE TABLE sdf_file("
                          "filename     varchar primary key,"
                          "lowest_cid   integer,"
                          "highest_cid  integer);")

    db_connection.execute("DROP TABLE IF EXISTS info")
    db_connection.execute("CREATE TABLE info("
                          "cid                     integer primary key not null,"
                          "iupac_name              varchar                     ,"
                          "iupac_inchi             varchar             not null,"
                          "iupac_inchikey          varchar                     ,"
                          "inchikey1               varchar                     ,"
                          "xlogp3                  real                        ,"
                          "monoisotopic_mass       real                not null,"
                          "molecular_formula       varchar             not null,"
                          "smiles                  varchar                     ,"
                          "can_smiles              varchar,"
                          "atom_def_stereo_count   integer,"
                          "atom_udef_stereo_count  integer,"
                          "bond_def_stereo_count   integer,"
                          "bond_udef_stereo_count  integer);")


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Usage: python %s <PUBCHEM_DB_BASEDIR> <RESET_DB>" % os.path.basename(sys.argv[0]))
        print("Example:")
        print("\t python %s /run/media/bach/Intenso/data/pubchem/ True" % os.path.basename(sys.argv[0]))
        exit(1)

    basedir = sys.argv[1]
    reset_database = eval(sys.argv[2])

    sdf_dir = basedir + "/sdf/"
    db_dir = basedir + "/db/"

    # Connect to the 'pubchem' database.
    conn = sqlite3.connect(db_dir + "/pubchem.db", isolation_level=None)
    try:
        # Initialize the database
        with conn:
            initialize_db(conn, reset=reset_database)

        # Get a list of all sdf-files available and reduce it to the ones still
        # needed to be processed.
        sdf_files = glob.glob(sdf_dir + "*.sdf")
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (before filtering): %d" % n_sdf_files)

        with conn:
            sdf_files = get_sdf_files_not_in_db(conn, sdf_files)
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (after filtering): %d" % n_sdf_files)

        if n_sdf_files > 0:
            # Iterate over the sdf-files and add them one by one
            for ii, sdf_fn in enumerate(sdf_files):
                print("Process sdf-file: %d/%d" % (ii + 1, n_sdf_files))

                # parse and insert current sdf-file
                with open(sdf_fn, "r") as sdf_file, conn:
                    insert_info_from_sdf_strings(conn, split_sdf_file(sdf_file))

                    # add current sdf-file to the list of completed sdf-files
                    conn.execute("INSERT INTO sdf_file (filename, lowest_cid, highest_cid) VALUES(?,?,?)", (
                        os.path.basename(sdf_fn),
                        os.path.basename(sdf_fn).split(".")[0].split("_")[1],
                        os.path.basename(sdf_fn).split(".")[0].split("_")[2]))

            # Make the column 'monoisotopic_mass' an index column. Later, when we query by
            # monoisotopic mass with ppm-window we _hugely_ speed up the query.
            with conn:
                conn.execute("CREATE INDEX idx_monoisotopic_mass ON info(monoisotopic_mass)")

            # Create an index on the inchikey1 to query stereo-isomers from the database given
            # their 2D structure only.
            with conn:
                conn.execute("CREATE INDEX idx_inchikey1 ON info(inchikey1)")

            # Create an index on the inchi's
            with conn:
                conn.execute("CREATE INDEX idx_inchi on info(iupac_inchi)")

    except sqlite3.ProgrammingError as err:
        print("Programming error: '" + err.args[0] + "'.")
    except sqlite3.DatabaseError as err:
        print("Database error: '" + err.args[0] + "'.")
    except IOError as err:
        print("An IOError occurred: '" + os.strerror(err.args[0]) + "'.")
    except Exception as err:
        print(err.args[0])
    finally:
        conn.close()
