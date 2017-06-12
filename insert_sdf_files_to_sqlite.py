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
"""

"""
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
Table creation: 'fps'
    create table fps (cid               integer primary key,
                      cdk_substructure  varchar,
                      pubchem           varchar,
                      klekota_roth      varchar,
                      fp3               varchar,
                      maccs             varchar,
                      fp4               varchar,
                      
                     );
"""

import pybel
import glob
import sqlite3
import re # Regular expression
import os 
from timeit import default_timer as timer

def insert_sdf (conn, sdf_dir):
    """
    Function to insert a set of sdf files into a database. The directory 
    containing the sdf files is screened and all contained sdf files are 
    processed.
    
    conn: Database connection
    sdf_dir: Path to the sdf-directory. 
    """
    
    for sdf_fn in glob.glob (sdf_dir + "*.sdf"):
        start = timer()
        
        # Open the sdf-file
        try:
            sdf_file  = open (sdf_fn)
            molecules = split_sdf_file (sdf_file)
            sdf_file.close()
        except IOError:
            print ("IOError when processing: " + os.path.basename (sdf_fn) + ".")
            continue
        
        # Insert the sdfs into the database
        try:
            # Commit all molecules in a block to decrease the number of 
            # transactions.
            conn.execute ("BEGIN")    
            # This version takes about 3s for 25000 molecules.
            conn.executemany ("INSERT INTO sdf VALUES(?,?)", molecules)
            conn.execute ("COMMIT")
        except sqlite3.Error as e:
            print ("An error occured: '" + e.args[0] +  "', when processing " + \
                   os.path.basename (sdf_fn) + ".")
            conn.execute ("ROLLBACK")
            continue
                
        end = timer()
        print ("Processed: " + sdf_fn + " - " + str (end - start) + "sec")
        
def insert_info (conn):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table. 
    
    conn: database connection
    """
    start = timer()
    
    try:
        conn.execute ("BEGIN")
             
        for row in conn.execute ("SELECT cid, sdf_string FROM sdf"):
            mol_cid = row[0]
            
            # Parse the sdf-string of the current molecule.        
            # This loop takes about 8s for 25000 molecules.
            mol_sdf = extract_info_from_sdf (str (row[1]))
           
            conn.execute ("INSERT INTO info \
                          (cid, name, exact_mass, inchi, smiles_canonical, \
                           smiles, molecular_formula) \
                          VALUES(?,?,?,?,?,?,?)",
                          (mol_cid, 
                           mol_sdf["PUBCHEM_IUPAC_NAME"],
                           mol_sdf["PUBCHEM_EXACT_MASS"],
                           mol_sdf["PUBCHEM_IUPAC_INCHI"],
                           mol_sdf["PUBCHEM_OPENEYE_CAN_SMILES"],
                           mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
                           mol_sdf["PUBCHEM_MOLECULAR_FORMULA"]))
            
            # One transaction encompasses 25000 molecules.
            if mol_cid % 25000 == 0:    
                conn.execute ("COMMIT")
                conn.execute ("BEGIN")
                
        conn.execute ("COMMIT")
    except sqlite3.Error as e:
        print ("An error occured: " + e.args[0])
        conn.execute ("ROLLBACK")
            
    end = timer()
    print ("Extraction and insertation of the information - " \
           + str (end - start) + "sec")
    
def insert_fps (conn):
    """
    Insert the fingerprint-vectors of different fingerprint definitions:
        [CSI:FingerID]
        - CDK Substructure fingerprints             - CDK package   (Java)
        - PubChem (CACTVS) fingerprints (SUBSKEYS)  - pubchem       (sdf-file)
        - Klekota–Roth fingerprints                 - CDK package   (Java)
        - FP3 fingerprints                          - Pybel package (python)
        - MACCS fingerprints                        - Pybel package (python)
        
        - Circular fingerprints                     - CDK package   (Java)
        
    Fingerprint storage strategy in the database:
        1) "nbits;one-bit_1,one-bit_2,..."
        2) String of zeros and ones of length nbits. 
    """    
    
    start = timer()
    
    try:
        conn.execute ("BEGIN")
             
        for row in conn.execute ("SELECT cid, inchi FROM info"):
            mol_cid = row[0]
            
            # Parse the sdf-string of the current molecule.        
            # This loop takes about 8s for 25000 molecules.
            mol_sdf = calculate_fingerprint_babel (str (row[1]))
           
            conn.execute ("INSERT INTO info \
                          (cid, name, exact_mass, inchi, smiles_canonical, \
                           smiles, molecular_formula) \
                          VALUES(?,?,?,?,?,?,?)",
                          (mol_cid, 
                           mol_sdf["PUBCHEM_IUPAC_NAME"],
                           mol_sdf["PUBCHEM_EXACT_MASS"],
                           mol_sdf["PUBCHEM_IUPAC_INCHI"],
                           mol_sdf["PUBCHEM_OPENEYE_CAN_SMILES"],
                           mol_sdf["PUBCHEM_OPENEYE_ISO_SMILES"],
                           mol_sdf["PUBCHEM_MOLECULAR_FORMULA"]))
            
            # One transaction encompasses 25000 molecules.
            if mol_cid % 25000 == 0:    
                conn.execute ("COMMIT")
                conn.execute ("BEGIN")
                
        conn.execute ("COMMIT")
    except sqlite3.Error as e:
        print ("An error occured: " + e.args[0])
        conn.execute ("ROLLBACK")
            
    end = timer()
    print ("Calculation and insertation of the fingerprints - " \
           + str (end - start) + "sec")
    
def extract_info_from_sdf (sdf_string):
    """
    Extract the information for a given molecules from its sdf-string. 
    
    sdf_string: string containing the molecule's sdf
    
    return: 
        dictionary 
    """
    
    info = {"PUBCHEM_IUPAC_NAME"         : "NULL",
            "PUBCHEM_EXACT_MASS"         : -1.0,
            "PUBCHEM_IUPAC_INCHI"        : "NULL",
            "PUBCHEM_OPENEYE_CAN_SMILES" : "NULL",
            "PUBCHEM_OPENEYE_ISO_SMILES" : "NULL",
            "PUBCHEM_MOLECULAR_FORMULA"  : "NULL"}
    
    lines  = sdf_string.split ("\n")
    n_line = len (lines)
    i      = 0
    
    while i < n_line:
        if lines[i].find("PUBCHEM_IUPAC_NAME") != -1:
            info["PUBCHEM_IUPAC_NAME"] = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_EXACT_MASS") != -1:
            info["PUBCHEM_EXACT_MASS"] = float (lines[i+1].strip())
            i = i+1
        if lines[i].find("PUBCHEM_IUPAC_INCHI") != -1:
            info["PUBCHEM_IUPAC_INCHI"] = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_OPENEYE_CAN_SMILES") != -1:
            info["PUBCHEM_OPENEYE_CAN_SMILES"] = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_OPENEYE_ISO_SMILES") != -1:
            info["PUBCHEM_OPENEYE_ISO_SMILES"] = lines[i+1].strip()
            i = i+1
        if lines[i].find("PUBCHEM_MOLECULAR_FORMULA") != -1:
            info["PUBCHEM_MOLECULAR_FORMULA"] = lines[i+1].strip()
            i = i+1
        i = i+1
        
    return (info)
        
def split_sdf_file (sdf_file):
    """
    Split a sdf-file containing the information for several molecules into a 
    list of (compound-id, sdf-molecule). 
    
    sdf_file: File-object pointing to the sdf-file.
    """
    sdf_data  = sdf_file.read()
    start     = 0
    molecules = []
    
    while True:
        end_pos = sdf_data.find ("$$$$", start)
        if end_pos == -1:
            break
        
        # Extract the sdf-string and the cid
        sdf_str = sdf_data[start:end_pos-1]
        sdf_str = sdf_str.replace ("'", "")
        cid     = int (re.findall ("<PUBCHEM_COMPOUND_CID>\n([0-9]+)", sdf_str)[0])
        
        molecules.append ((cid, sdf_str))
        
        start = end_pos + 5
        
    return (molecules)
            
sdf_dir = "/home/bach/mnt/on_triton/project/pubchem_local/sandbox/"  
db_dir  = "/home/bach/mnt/on_triton/project/pubchem_local/sandbox/"

# Connect to the 'pubchem' database and get a cursor.
conn = sqlite3.connect (db_dir + "pubchem")
conn.isolation_level = None

#insert_sdf (conn, sdf_dir)
insert_info (conn)

conn.close()


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