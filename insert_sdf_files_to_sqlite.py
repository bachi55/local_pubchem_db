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

# import pybel
import queue
import glob
import sqlite3
import re # Regular expression
import os 
import time
import numpy as np
from timeit import default_timer as timer


def insert_sdf (conn, sdf_dir, max_num_attempts = 5):
    """
    Function to insert a set of sdf files into a database. The directory 
    containing the sdf files is screened and all contained sdf files are 
    processed.
    
    :param conn: sqlite3.Connection, to the pubchem database
    :param sdf_dir: string, Path to the sdf-directory containing the pubchem sdf files
    """
   
    # Get a list of all sdf-files available and reduce it to the ones still 
    # needed to be processed.
    sdf_files = glob.glob (sdf_dir + "*.sdf")
    sdf_files = get_sdf_files_not_in_DB (conn, sdf_files)
    
    # Add the sdf-files into a Queue to process them. Each sdf-file is stored
    # together with an integer indicating how often this file was already tried 
    # to be added to the DB. It can happen that there is an IO-Error and we need 
    # to retry to add a sdf-file. However, we want to restrict the maximum 
    # number of tries.
    que_sdf_files = queue.Queue()
    for sdf_file in sdf_files:
        que_sdf_files.put ((sdf_file, 0))
        
    que_sdf_files_failed = queue.Queue()
    
    n_files = que_sdf_files.qsize()
    print ("Sdf-files to process: %d" % n_files)

    while not que_sdf_files.empty():                
        start = timer()

        # Open the sdf-file
        sdf_fn, num_attempts = que_sdf_files.get()
        try:
            with open (sdf_fn) as sdf_file:
                molecules = split_sdf_file (sdf_file)
        except IOError as e:
            print ("An IOError occured: '" + os.strerror (e.args[0]) + "' when processing " \
                   + os.path.basename (sdf_fn) + ".")
            
            # Increase the number of attempts due to an unsuccessful try.
            num_attempts = num_attempts + 1
            if num_attempts < max_num_attempts:
                # Enqueue the sdf-file again for a later try.
                que_sdf_files.put ((sdf_fn, num_attempts))
            else: 
                que_sdf_files_failed.put (sdf_fn)
            
            continue
        
        # Insert the sdfs into the database
        try:
            # Commit all molecules in a block to decrease the number of 
            # transactions.
            conn.execute ("BEGIN")    
            # This version takes about 3s for 25000 molecules.
            conn.executemany ("INSERT INTO sdf VALUES(?,?)", molecules)
            
            conn.execute ("INSERT INTO sdf_file                \
                           (filename, lowest_cid, highest_cid) \
                           VALUES(?,?,?)",
                          (os.path.basename (sdf_fn),
                           os.path.basename (sdf_fn).split(".")[0].split("_")[1],
                           os.path.basename (sdf_fn).split(".")[0].split("_")[2]))
            
            conn.execute ("COMMIT")
        except sqlite3.ProgrammingError as e: 
            # Error caused by the programmer, e.g. incorrect query. This error 
            # indicates a logical error, which needs to be fixed. 
            print ("A programming error occured: '" + e.args[0] + "', when \
                    processing " + os.path.basename (sdf_fn) + ".")
            conn.execute ("ROLLBACK")
            raise
        except (sqlite3.OperationalError, sqlite3.InternalError) as e: 
            # Error caused by the database, e.g. lost connection, memory 
            # allocation error, DB out of sync. We should re-try to add the 
            # current sdf-file.
            print ("An operational / internal error occured: '" + e.args[0] + "', when \
                    processing " + os.path.basename (sdf_fn) + ".")
            conn.execute ("ROLLBACK")
            
            # Increase the number of attempts due to an unsuccessful try.
            num_attempts = num_attempts + 1
            if num_attempts < max_num_attempts:
                # Enqueue the sdf-file again for a later try.
                que_sdf_files.put ((sdf_fn, num_attempts))
            else: 
                que_sdf_files_failed.put (sdf_fn)
            
            continue 
                
        end = timer()
        print ("Processed: " + os.path.basename (sdf_fn) + \
               " - " + str (round (end - start, 2)) + "sec (" + \
               str (que_sdf_files.qsize()) + " remain)")

def get_sdf_files_not_in_DB (conn, sdf_fn_in_folder):
    """
    Returns the sdf_file names (full path) which are not already in the DB.
    
    :param conn: sqlite3.Connection, database connection
    :param sdf_fn_in_folder: list of string, filenames of the sdf-files.
    """
    sdf_files_in_DB = conn.execute ("SELECT filename FROM sdf_file").fetchall()
    sdf_files_in_DB = [str (x[0]) for x in sdf_files_in_DB]
                             
    return list (filter (lambda x: os.path.basename (x) not in sdf_files_in_DB, sdf_fn_in_folder))
        
def insert_info (conn, chunk_size = 10000):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table. 
    
    :param conn: sqlite3.Connection, database connection
    :param chunk_size: integer, number of CIDs to insert into the DB with a COMMIT.
    """
    start = timer()
    
    try:
        inserted_cids = 0

        conn.execute ("BEGIN")

        for cid, sdf in conn.execute ("SELECT cid, sdf_string FROM sdf"):

            # Parse the sdf-string of the current molecule.        
            # This loop takes about 8s for 25000 molecules.
            mol_sdf = extract_info_from_sdf (sdf)

            assert (cid == mol_sdf["PUBCHEM_COMPOUND_CID"])

            conn.execute ("INSERT INTO info \
                          (cid, iupac_name, iupac_inchi, iupac_inchikey, xlogp3, \
                           monoisotopic_mass, molecular_formula, molecular_weight, smiles, isotopic_atom_count) \
                          VALUES(?,?,?,?,?,?,?,?,?)",
                          (cid,
                           mol_sdf["PUBCHEM_IUPAC_NAME"],
                           mol_sdf["PUBCHEM_IUPAC_INCHI"],
                           mol_sdf["PUBCHEM_IUPAC_INCHIKEY"],
                           mol_sdf["PUBCHEM_XLOGP3_AA"],
                           mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"],
                           mol_sdf["PUBCHEM_MOLECULAR_FORMULA"],
                           mol_sdf["PUBCHEM_MOLECULAR_WEIGHT"],
                           mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
                           mol_sdf["PUBCHEM_ISOTOPIC_ATOM_COUNT"]))
            
            # One transaction encompasses 'chunk_size' molecules.
            if inserted_cids == chunk_size:
                inserted_cids = 0
                print("Current CID: %d" % cid)

                conn.execute ("COMMIT")
                conn.execute ("BEGIN")
            else:
                inserted_cids = inserted_cids + 1
                
        conn.execute ("COMMIT")
    except sqlite3.Error as e:
        print ("An error occured: " + e.args[0])
        conn.execute ("ROLLBACK")
            
    end = timer()
    print ("Extraction and insertation of the information took %.3fsec" % (end - start))

def insert_info_from_sdf_strings (conn, l_cid_sdf):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table.

    :param conn: sqlite3.Connection, database connection
    :param l_cid_sdf: list of (cid, sdf-string)-tuples
    """
    start = timer()

    conn.execute ("BEGIN")

    for cid, sdf in l_cid_sdf:
        mol_sdf = extract_info_from_sdf (sdf)

        assert (cid == mol_sdf["PUBCHEM_COMPOUND_CID"])

        conn.execute ("INSERT INTO info \
                      (cid, iupac_name, iupac_inchi, iupac_inchikey, xlogp3, \
                       monoisotopic_mass, molecular_formula, molecular_weight, smiles, isotopic_atom_count) \
                      VALUES(?,?,?,?,?,?,?,?,?,?)",
                      (cid,
                       mol_sdf["PUBCHEM_IUPAC_NAME"],
                       mol_sdf["PUBCHEM_IUPAC_INCHI"],
                       mol_sdf["PUBCHEM_IUPAC_INCHIKEY"],
                       mol_sdf["PUBCHEM_XLOGP3_AA"],
                       mol_sdf["PUBCHEM_MONOISOTOPIC_WEIGHT"],
                       mol_sdf["PUBCHEM_MOLECULAR_FORMULA"],
                       mol_sdf["PUBCHEM_MOLECULAR_WEIGHT"],
                       mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
                       mol_sdf["PUBCHEM_ISOTOPIC_ATOM_COUNT"]))

    conn.execute ("COMMIT")

    end = timer()
    print ("Extraction and insertion of the information took %.3fsec" % (end - start))

def extract_info_from_sdf (sdf):
    """
    Extract the information for a given molecules from its sdf-string. 
    
    :param sdf: string, containing the molecule's sdf
    
    return: dictionary, containing the extracted information.
    """
    
    info = {
        "PUBCHEM_COMPOUND_CID"        : np.nan, # integer
        "PUBCHEM_IUPAC_NAME"          : "NULL", # string
        "PUBCHEM_IUPAC_INCHI"         : "NULL", # string
        "PUBCHEM_IUPAC_INCHIKEY"      : "NULL", # string
        "PUBCHEM_XLOGP3_AA"           : np.nan, # float
        "PUBCHEM_MONOISOTOPIC_WEIGHT" : np.nan, # float
        "PUBCHEM_MOLECULAR_FORMULA"   : "NULL", # string
        "PUBCHEM_MOLECULAR_WEIGHT"    : np.nan, # float
        "PUBCHEM_OPENEYE_ISO_SMILES"  : "NULL", # string
        "PUBCHEM_ISOTOPIC_ATOM_COUNT" : np.nan  # integer
    }
    found = {key: False for key in info.keys()}

    lines  = sdf.split ("\n")
    n_line = len (lines)
    i      = 0
    
    while not all (found.values()) and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i = i+1
            continue

        if not found["PUBCHEM_COMPOUND_CID"] and lines[i].find("PUBCHEM_COMPOUND_CID") != -1:
            found["PUBCHEM_COMPOUND_CID"] = True
            info["PUBCHEM_COMPOUND_CID"] = int (lines[i+1].strip())
            i = i+2
            continue

        if not found["PUBCHEM_IUPAC_NAME"] and lines[i].find("PUBCHEM_IUPAC_NAME") != -1:
            found["PUBCHEM_IUPAC_NAME"] = True
            info["PUBCHEM_IUPAC_NAME"] = lines[i+1].strip()
            i = i+2
            continue

        if not found["PUBCHEM_IUPAC_INCHI"] and lines[i].find("PUBCHEM_IUPAC_INCHI") != -1:
            found["PUBCHEM_IUPAC_INCHI"] = True
            info["PUBCHEM_IUPAC_INCHI"] = lines[i+1].strip()
            i = i+2
            continue

        if not found["PUBCHEM_IUPAC_INCHIKEY"] and lines[i].find("PUBCHEM_IUPAC_INCHIKEY") != -1:
            found["PUBCHEM_IUPAC_INCHIKEY"] = True
            info["PUBCHEM_IUPAC_INCHIKEY"] = lines[i+1].strip()
            i = i+2
            continue

        if not found["PUBCHEM_XLOGP3_AA"] and lines[i].find("PUBCHEM_XLOGP3_AA") != -1:
            found["PUBCHEM_XLOGP3_AA"] = True
            info["PUBCHEM_XLOGP3_AA"] = float (lines[i+1].strip())
            i = i+2
            continue

        if not found["PUBCHEM_MONOISOTOPIC_WEIGHT"] and lines[i].find("PUBCHEM_MONOISOTOPIC_WEIGHT") != -1:
            found["PUBCHEM_MONOISOTOPIC_WEIGHT"] = True
            info["PUBCHEM_MONOISOTOPIC_WEIGHT"] = float (lines[i+1].strip())
            i = i+2
            continue

        if not found["PUBCHEM_MOLECULAR_FORMULA"] and lines[i].find("PUBCHEM_MOLECULAR_FORMULA") != -1:
            found["PUBCHEM_MOLECULAR_FORMULA"] = True
            info["PUBCHEM_MOLECULAR_FORMULA"] = lines[i+1].strip()
            i = i+2
            continue

        if not found["PUBCHEM_MOLECULAR_WEIGHT"] and lines[i].find("PUBCHEM_MOLECULAR_WEIGHT") != -1:
            found["PUBCHEM_MOLECULAR_WEIGHT"] = True
            info["PUBCHEM_MOLECULAR_WEIGHT"] = float (lines[i+1].strip())
            i = i+2
            continue

        if not found["PUBCHEM_OPENEYE_ISO_SMILES"] and lines[i].find("PUBCHEM_OPENEYE_ISO_SMILES") != -1:
            found["PUBCHEM_OPENEYE_ISO_SMILES"] = True
            info["PUBCHEM_OPENEYE_ISO_SMILES"] = lines[i+1].strip()
            i = i+2
            continue

        if not found["PUBCHEM_ISOTOPIC_ATOM_COUNT"] and lines[i].find("PUBCHEM_ISOTOPIC_ATOM_COUNT") != -1:
            found["PUBCHEM_ISOTOPIC_ATOM_COUNT"] = True
            info["PUBCHEM_ISOTOPIC_ATOM_COUNT"] = int (lines[i+1].strip())
            i = i+2
            continue

        i = i+1
        
    return info
        
def split_sdf_file (sdf_file):
    """
    Split a sdf-file containing the information for several molecules into a 
    list of (cid, sdf-string).
    
    :param sdf_file: file stream, pointing to the sdf-file. Result of 'open'

    :return: list of (cid, sdf-string)-tuples
    """
    sdf_lines  = sdf_file.read()
    start     = 0
    sdfs = []

    while True:
        end_pos = sdf_lines.find ("$$$$", start)
        if end_pos == -1:
            break

        # Extract the sdf-string and the cid
        sdf_str = sdf_lines[start:end_pos-1]
        sdf_str = sdf_str.replace ("'", "")
        cid     = int (re.findall ("<PUBCHEM_COMPOUND_CID>\n([0-9]+)", sdf_str)[0])

        sdfs.append ((cid, sdf_str))

        start = end_pos + 5

    return sdfs
         
def initialize_DB (conn, add_sdf = True, add_info = True, add_sdf_file = True, reset = False):
    """
    Initialization of the DB: 
        - Remove existing tables if desired.
        - Create not present tables. 
        
    :param conn: sqlite3.Connection, connection to the database to modify.
    :param reset: boolean, If 'True' tables in the database are dropped and recreated.
    """

    # Check whether 'sdf' exists:
    if (add_sdf):
        cur = conn.execute ("SELECT name FROM sqlite_master WHERE type='table' AND name='sdf'")
        len_cur = len (cur.fetchall())

        if reset & (len_cur > 0):
            conn.execute ("DROP TABLE sdf")
            len_cur = 0

        if len_cur == 0:
            conn.execute ("CREATE TABLE sdf(cid integer primary key, sdf_string varchar)")

    # Check whether 'info' exists:
    if (add_info):
        cur = conn.execute ("SELECT name FROM sqlite_master WHERE type='table' AND name='info'")
        len_cur = len (cur.fetchall())

        if reset & (len_cur > 0):
            conn.execute ("DROP TABLE info")
            len_cur = 0

        if len_cur == 0:
            conn.execute ("CREATE TABLE info(cid                 integer primary key, \
                                             iupac_name          varchar,             \
                                             iupac_inchi         varchar,             \
                                             iupac_inchikey      varchar,             \
                                             xlogp3              real,                \
                                             monoisotopic_mass   real,                \
                                             molecular_formula   varchar,             \
                                             molecular_weight    real,                \
                                             smiles              varchar,             \
                                             isotopic_atom_count integer              \
                                             );")

    # Check whether 'sdf_file' exsits:
    if (add_sdf_file):
        cur = conn.execute ("SELECT name FROM sqlite_master WHERE type='table' AND name='sdf_file'")
        len_cur = len (cur.fetchall())

        if reset & (len_cur > 0):
            conn.execute ("DROP TABLE sdf_file")
            len_cur = 0

        if len_cur == 0:
            conn.execute ("CREATE TABLE sdf_file(filename                varchar primary key, \
                                                 lowest_cid              integer,             \
                                                 highest_cid             integer              \
                                                );")

if __name__ == "__main__":

    sdf_dir = "/m/cs/project/kepaco/pubchem_local/compounds_sdf/"
    db_dir  = "/m/cs/project/kepaco/pubchem_local/db/"

    sdf_dir_sandbox = "/media/intenso/data/sandbox/"
    db_dir_sandbox  = "/media/intenso/data/sandbox/"
    reset_database = True

    # Connect to the 'pubchem' database.
    conn = sqlite3.connect (db_dir_sandbox + "/pubchem.db", isolation_level = None)
    try:
        # Initialize the database
        with conn:
             initialize_DB (conn, add_sdf = False, reset = reset_database)

        # Get a list of all sdf-files available and reduce it to the ones still
        # needed to be processed.
        sdf_files = glob.glob (sdf_dir_sandbox + "*.sdf")
        print ("Sdf-files to process (before filtering): %d" % len (sdf_files))

        with conn:
            sdf_files = get_sdf_files_not_in_DB (conn, sdf_files)

        print ("Sdf-files to process: %d" % len (sdf_files))

        # Iterate over the sdf-fles and add them one by one
        for ii, sdf_fn in enumerate (sdf_files):
            print ("Process sdf-file: %d/%d" % (ii + 1, len (sdf_files)))

            # parse and insert current sdf-file
            with open (sdf_fn) as sdf_file, conn:
                insert_info_from_sdf_strings (conn, split_sdf_file (sdf_file))

                # add current sdf-file to the list of completed sdf-files
                conn.execute ("INSERT INTO sdf_file               \
                              (filename, lowest_cid, highest_cid) \
                              VALUES(?,?,?)",
                              (os.path.basename (sdf_fn),
                               os.path.basename (sdf_fn).split(".")[0].split("_")[1],
                               os.path.basename (sdf_fn).split(".")[0].split("_")[2]))

    except sqlite3.ProgrammingError as e:
        print ("Programming error: '" + e.args[0] + "'.")
    except sqlite3.DatabaseError as e:
        print ("Database error: '" + e.args[0] + "'.")
    except IOError as e:
        print ("An IOError occured: '" + os.strerror (e.args[0]) + "'.")
    except Exception as e:
        print (e.args[0])
    finally:
        conn.close()


    # Insert the sdf string from the sdf-files
    # insert_sdf (conn, sdf_dir_sandbox)
    # insert_info (conn, 25000)





#insert_sdf (conn, sdf_dir)


#cur.execute ("INSERT INTO sdf VALUES(?,?)", 
#             (molob.data["PUBCHEM_COMPOUND_CID"],
#              molob.write ("sdf")))
#conn.commit()

            # This version takes about 35s for 25000 molecules. 
#            for mol in pybel.readfile ("sdf", sdf_fn):
#                print ("Pubchem-ID: " + mol.data["PUBCHEM_COMPOUND_CID"])
#                cur.execute ("INSERT INTO sdf VALUES(?,?)", 
#                             (mol.data["PUBCHEM_COMPOUND_CID"],
#                              mol.write ("sdf")))

            # This loop takes about 23s for 25000 molecules.
#            mol_sdf = pybel.readstring ("sdf", str (row[1]))
#            conn.execute ("INSERT INTO info VALUES(?,?,?,?,?,?)",
#                          (mol_cid, 
#                           mol_sdf.data["PUBCHEM_EXACT_MASS"],
#                           mol_sdf.data["PUBCHEM_IUPAC_INCHI"],
#                           mol_sdf.data["PUBCHEM_OPENEYE_CAN_SMILES"],
#                           mol_sdf.data["PUBCHEM_OPENEYE_ISO_SMILES"],
#                           mol_sdf.data["PUBCHEM_MOLECULAR_FORMULA"]))
# 
