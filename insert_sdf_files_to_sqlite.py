import queue
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

Database layout: 
    DB: pubchem --> TABLE-1: 'sdf', contains the sdf files for each compound-id
                --> TABLE-2: 'info', contains several information for each 
                                     compound-id, e.g. exact mass, mol-formula, 
                                     ...
                --> TABLE-3: 'fp', contains the fingerprint-vectors for each 
                                   compound.
                --> TABLE-4: 'sdf_file', contains the filenames of the 
                                         sdf-files successfully added to the 
                                         DB.
"""

"""
Table creation: 'sdf_file'
    create table sdf_file(filename                 varchar primary key,
                          last_modification        real, 
                          last_modification_ctime  varchar, 
                          lowest_cid               integer, 
                          highest_cid              integer
                         );
Table creation: 'sdf'
    create table sdf(cid integer primary key, sdf_string varchar);
Table creation: 'info'
    create table info(cid               integer primary key, 
                      name              varchar,
                      exact_mass        real, 
                      inchi             varchar,
                      smiles_canonical  varchar,
                      smiles            varchar,
                      molecular_formula varchar
                      );
Table creation: 'fp'
    create table fp(cid               integer primary key,
                    cdk_substructure  varchar,
                    pubchem           varchar,
                    klekota_roth      varchar,
                    fp3               varchar,
                    maccs             varchar,
                    fp4               varchar,
                   );
"""


def insert_sdf(db_connection, sdf_directory, max_num_attempts=5):
    """
    Function to insert a set of sdf files into a database. The directory 
    containing the sdf files is screened and all contained sdf files are 
    processed.
    
    :param db_connection: sqlite3.Connection, to the pubchem database
    :param sdf_directory: string, Path to the sdf-directory containing the pubchem sdf files
    :param max_num_attempts: integer, number of tries to insert a single sdf file into the DB in case of an error
    """
   
    # Get a list of all sdf-files available and reduce it to the ones still 
    # needed to be processed.
    sdf_filenames = glob.glob(sdf_directory + "*.sdf")
    sdf_filenames = get_sdf_files_not_in_db(db_connection, sdf_filenames)
    
    # Add the sdf-files into a Queue to process them. Each sdf-file is stored
    # together with an integer indicating how often this file was already tried 
    # to be added to the DB. It can happen that there is an IO-Error and we need 
    # to retry to add a sdf-file. However, we want to restrict the maximum 
    # number of tries.
    que_sdf_files = queue.Queue()
    for sdf_filename in sdf_filenames:
        que_sdf_files.put((sdf_filename, 0))
        
    que_sdf_files_failed = queue.Queue()
    
    n_files = que_sdf_files.qsize()
    print("Sdf-files to process: %d" % n_files)

    while not que_sdf_files.empty():                
        start = timer()

        # Open the sdf-file
        sdf_filename, num_attempts = que_sdf_files.get()
        try:
            with open(sdf_filename) as sdf_file_stream:
                molecules = split_sdf_file(sdf_file_stream)
        except IOError as e:
            print("An IOError occured: '" + os.strerror(e.args[0]) + "' when processing " +
                  os.path.basename(sdf_filename) + ".")
            
            # Increase the number of attempts due to an unsuccessful try.
            num_attempts = num_attempts + 1
            if num_attempts < max_num_attempts:
                # Enqueue the sdf-file again for a later try.
                que_sdf_files.put((sdf_filename, num_attempts))
            else: 
                que_sdf_files_failed.put(sdf_filename)
            
            continue
        
        # Insert the sdfs into the database
        try:
            # Commit all molecules in a block to decrease the number of 
            # transactions.
            db_connection.execute("BEGIN")
            # This version takes about 3s for 25000 molecules.
            db_connection.executemany("INSERT INTO sdf VALUES(?,?)", molecules)
            
            db_connection.execute("INSERT INTO sdf_file (filename, lowest_cid, highest_cid) VALUES(?,?,?)", (
                os.path.basename(sdf_filename),
                os.path.basename(sdf_filename).split(".")[0].split("_")[1],
                os.path.basename(sdf_filename).split(".")[0].split("_")[2]))
            
            db_connection.execute("COMMIT")
        except sqlite3.ProgrammingError as error:
            # Error caused by the programmer, e.g. incorrect query. This error 
            # indicates a logical error, which needs to be fixed. 
            print("A programming error occurred: '" + error.args[0] + "', when processing " +
                  os.path.basename(sdf_filename) + ".")
            db_connection.execute("ROLLBACK")
            raise
        except (sqlite3.OperationalError, sqlite3.InternalError) as error:
            # Error caused by the database, e.g. lost connection, memory 
            # allocation error, DB out of sync. We should re-try to add the 
            # current sdf-file.
            print("An operational / internal error occurred: '" + error.args[0] + "', when processing " +
                  os.path.basename(sdf_filename) + ".")
            db_connection.execute("ROLLBACK")
            
            # Increase the number of attempts due to an unsuccessful try.
            num_attempts = num_attempts + 1
            if num_attempts < max_num_attempts:
                # Enqueue the sdf-file again for a later try.
                que_sdf_files.put((sdf_filename, num_attempts))
            else: 
                que_sdf_files_failed.put(sdf_filename)
            
            continue 
                
        end = timer()
        print("Processed: " + os.path.basename(sdf_filename) + " - " + str(round(end - start, 2)) + "sec (" +
              str(que_sdf_files.qsize()) + " remain)")


def get_sdf_files_not_in_db(db_connection, sdf_fn_in_folder):
    """
    Returns the sdf_file names (full path) which are not already in the DB.
    
    :param db_connection: sqlite3.Connection, database connection
    :param sdf_fn_in_folder: list of string, filenames of the sdf-files.
    """
    sdf_files_in_db = db_connection.execute("SELECT filename FROM sdf_file").fetchall()
    sdf_files_in_db = [str(x[0]) for x in sdf_files_in_db]
                             
    return list(filter(lambda x: os.path.basename(x) not in sdf_files_in_db, sdf_fn_in_folder))


def insert_info(db_connection, chunk_size=10000):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table. 
    
    :param db_connection: sqlite3.Connection, database connection
    :param chunk_size: integer, number of CIDs to insert into the DB with a COMMIT.
    """
    start = timer()
    
    try:
        inserted_cids = 0

        db_connection.execute("BEGIN")

        for cid, sdf in db_connection.execute("SELECT cid, sdf_string FROM sdf"):

            # Parse the sdf-string of the current molecule.        
            # This loop takes about 8s for 25000 molecules.
            mol_sdf = extract_info_from_sdf(sdf)

            assert (cid == mol_sdf["PUBCHEM_COMPOUND_CID"])

            db_connection.execute("INSERT INTO info \
                          (cid, iupac_name, iupac_inchi, iupac_inchikey, xlogp3, \
                           monoisotopic_mass, molecular_formula, smiles) \
                          VALUES(?,?,?,?,?,?,?,?)", (
                cid,
                mol_sdf["PUBCHEM_IUPAC_NAME"],
                mol_sdf["PUBCHEM_IUPAC_INCHI"],
                mol_sdf["PUBCHEM_IUPAC_INCHIKEY"],
                mol_sdf["PUBCHEM_XLOGP3*"],
                mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"],
                mol_sdf["PUBCHEM_MOLECULAR_FORMULA"],
                mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"]))
            
            # One transaction encompasses 'chunk_size' molecules.
            if inserted_cids == chunk_size:
                inserted_cids = 0
                print("Current CID: %d" % cid)

                db_connection.execute("COMMIT")
                db_connection.execute("BEGIN")
            else:
                inserted_cids = inserted_cids + 1
                
        db_connection.execute("COMMIT")
    except sqlite3.Error as error:
        print("An error occured: " + error.args[0])
        db_connection.execute("ROLLBACK")
            
    end = timer()
    print("Extraction and insertation of the information took %.3fsec" % (end - start))


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

        db_connection.execute("INSERT INTO info \
                      (cid, iupac_name, iupac_inchi, iupac_inchikey, inchikey1, xlogp3, \
                       monoisotopic_mass, molecular_formula, smiles, atom_def_stereo_count, atom_udef_stereo_count, \
                       bond_def_stereo_count, bond_udef_stereo_count) \
                      VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", (
            cid,
            mol_sdf["PUBCHEM_IUPAC_NAME"],
            mol_sdf["PUBCHEM_IUPAC_INCHI"],
            mol_sdf["PUBCHEM_IUPAC_INCHIKEY"],
            mol_sdf["INCHIKEY1"],
            mol_sdf["PUBCHEM_XLOGP3*"],
            mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"],
            mol_sdf["PUBCHEM_MOLECULAR_FORMULA"],
            mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
            mol_sdf["PUBCHEM_ATOM_DEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_ATOM_UDEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_BOND_DEF_STEREO_COUNT"],
            mol_sdf["PUBCHEM_BOND_UDEF_STEREO_COUNT"]
        ))

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


def initialize_db(db_connection, add_sdf=True, add_info=True, add_sdf_file=True, reset=False):
    """
    Initialization of the DB: 
        - Remove existing tables if desired.
        - Create not present tables. 
        
    :param db_connection: sqlite3.Connection, connection to the database to modify.
    :param add_sdf: boolean, should the sdf table be added containing the fill sdf-string for each CID
    :param add_info: boolean, should the info table be added containing information from the parsed sdf-strings for each
        CID
    :param add_sdf_file: boolean, should the sdf_file table be added containing the filenames of all processed sdf-files
    :param reset: boolean, If 'True' tables in the database are dropped and recreated.
    """

    # TODO: "CREATE TABLE IF NOT EXISTS ..." circumvents explicit determination of table existence.

    # Check whether 'sdf' exists:
    if add_sdf:
        cur = db_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sdf'")
        len_cur = len(cur.fetchall())

        if reset & (len_cur > 0):
            db_connection.execute("DROP TABLE sdf")
            len_cur = 0

        if len_cur == 0:
            db_connection.execute("CREATE TABLE sdf(cid integer primary key, sdf_string varchar)")

    # Check whether 'info' exists:
    if add_info:
        cur = db_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='info'")
        len_cur = len(cur.fetchall())

        if reset & (len_cur > 0):
            db_connection.execute("DROP TABLE info")
            len_cur = 0

        if len_cur == 0:
            db_connection.execute("CREATE TABLE info(cid         integer primary key not null, \
                                             iupac_name          varchar                     , \
                                             iupac_inchi         varchar             not null, \
                                             iupac_inchikey      varchar                     , \
                                             inchikey1           varchar                     , \
                                             xlogp3              real                        , \
                                             monoisotopic_mass   real                not null, \
                                             molecular_formula   varchar             not null, \
                                             smiles              varchar                     , \
                                             atom_def_stereo_count integer, \
                                             atom_udef_stereo_count integer, \
                                             bond_def_stereo_count integer, \
                                             bond_udef_stereo_count integer);")


    # Check whether 'sdf_file' exsits:
    if add_sdf_file:
        cur = db_connection.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='sdf_file'")
        len_cur = len(cur.fetchall())

        if reset & (len_cur > 0):
            db_connection.execute("DROP TABLE sdf_file")
            len_cur = 0

        if len_cur == 0:
            db_connection.execute("CREATE TABLE sdf_file(filename        varchar primary key, \
                                                 lowest_cid              integer,             \
                                                 highest_cid             integer              \
                                                );")


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
            initialize_db(conn, add_sdf=False, reset=reset_database)

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
