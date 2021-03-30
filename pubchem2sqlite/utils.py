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

import re
import os
import json
import gzip
import glob
import sqlite3
import traceback

from collections import OrderedDict
from timeit import default_timer as timer


def _as_dtype(val, dtype):
    """
    Convert a value represented using a string variable to the desired output data-type.

    :param val: string, representation of the value to convert

    :param dtype: string, desired output data-type

    :return: input values converted to the output data-type
    """
    if dtype in ["integer", "int"]:
        oval = int(val)
    elif dtype in ["real", "float", "double"]:
        oval = float(val)
    elif dtype in ["varchar", "character", "text"]:
        oval = val
    else:
        raise ValueError("Invalid dtype: %s." % dtype)

    return oval


def extract_info_from_sdf(sdf, db_specs):
    """
    Extract the information for a given molecules from its sdf-string.

    :param sdf: string, containing the molecule's sdf

    :param db_specs: dictionary, containing the database specifications (see 'default_db_layout.json')

    :return: OrderedDict, containing the extracted information.
    """
    infos = OrderedDict([(k, None) for k in db_specs["columns"]])
    missing_infos = set(infos.keys())

    #####################################################
    # TODO: This can be moves somewhere else. Its a global setup for all sdf-files.
    # type of the information associated with the SD-tag
    dtypes = {}
    create_likes = {}
    sdtags2info = {}
    for k, specs in db_specs["columns"].items():
        dtypes[k] = specs["DTYPE"].lower()

        # Check for value transformation function
        if "CREATE_LIKE" in specs:
            create_likes[k] = eval(specs["CREATE_LIKE"])

        for sdtag in ["> <%s>" % v for v in specs["SD_TAG"]]:
            try:
                sdtags2info[sdtag].append(k)
            except KeyError:
                sdtags2info[sdtag] = [k]
    #####################################################

    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0

    while missing_infos and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i += 1
            continue

        if lines[i].startswith("> ") and (lines[i] in sdtags2info):
            for info in sdtags2info[lines[i]]:
                val = _as_dtype(lines[i + 1], dtypes[info])

                # Apply value transformation if provided
                if info in create_likes:
                    val = create_likes[info](val)

                infos[info] = val

                missing_infos.remove(info)

            i += 2

        i += 1

    return infos


def insert_info_from_sdf_strings(db_connection, db_specs, iter_cid_sdf):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table.

    :param db_connection: sqlite3.Connection, database connection

    :param db_specs: dictionary, containing the database specifications (see 'default_db_layout.json')

    :param iter_cid_sdf: list of (cid, sdf-string)-tuples

    :return: scalar, number of inserted rows (compounds)
    """
    start = timer()

    col_names = ",".join(db_specs["columns"].keys())
    placeholders = ",".join(["?"] * len(db_specs["columns"]))

    # Get all columns that are requested to be NOT NULL
    not_null_cols = set(col for col, specs in db_specs["columns"].items() if specs.get("NOT_NULL", False))

    n_inserted = 0

    for _, sdf in iter_cid_sdf:
        # Extract information
        mol_sdf = extract_info_from_sdf(sdf, db_specs)

        # Skip PubChem entries that do not provide all information requested to be NOT NULL
        skip_sdf = False
        for col in not_null_cols:
            if mol_sdf[col] is None:
                skip_sdf = True
                break
        if skip_sdf:
            continue

        # Insert extracted information to the database
        db_connection.execute("INSERT INTO compounds (%s) VALUES (%s)" % (col_names, placeholders),
                              tuple(v for v in mol_sdf.values()))
        n_inserted += 1

    end = timer()
    print("Extraction and insertion of the information took %.3fsec" % (end - start))

    return n_inserted


def load_db_specifications(fn):
    """
    Loads a database layout specified in a json-file into a dictionary.

    :param fn: string, path to the json-file specifying the database layout. (see 'default_db_layout.json')

    :return: dictionary
    """
    with open(fn, "r") as json_file:
        specs = json.loads(json_file.read(), object_pairs_hook=OrderedDict)
    return specs


def get_column_stmt(column_specs):
    stmt_columns = []

    has_primary_key = False  # Primary keys must (if) be defined on a single column.

    for name, spec in column_specs.items():
        new_col = [name, spec["DTYPE"]]

        if spec.get("NOT_NULL", False) or spec.get("PRIMARY_KEY", False):
            new_col.append("not null")

        if spec.get("PRIMARY_KEY", False):
            if has_primary_key:
                raise ValueError("Primary keys must be defined on a single column.")

            new_col.append("primary key")
            has_primary_key = True

        stmt_columns.append(" ".join(new_col))

    return ",".join(stmt_columns)


def initialize_db(db_connection, specs, reset=False):
    """
    Initialization of the DB:
        - Remove existing tables if desired.
        - Create not present tables.

    :param db_connection: sqlite3.Connection, connection to the database to modify.

    :param specs: dictionary, specifying the database layout, e.g. columns, primary key and indices.
        See: 'load_db_specifications' for details.

    :param reset: boolean, If 'True' tables in the database are dropped and recreated.
    """
    if reset:
        db_connection.execute("DROP TABLE IF EXISTS sdf_file")
        db_connection.execute("DROP TABLE IF EXISTS compounds")

    # Table to keep track about the sdf-chucks (files) included
    db_connection.execute("CREATE TABLE IF NOT EXISTS sdf_file("
                          "filename          VARCHAR NOT NULL PRIMARY KEY,"
                          "lowest_cid        INTEGER NOT NULL,"
                          "highest_cid       INTEGER NOT NULL,"
                          "date_added        VARCHAR NOT NULL,"
                          "n_compounds       INTEGER NOT NULL)")

    # Table containing the compounds
    db_connection.execute("CREATE TABLE IF NOT EXISTS compounds(%s)" % get_column_stmt(specs["columns"]))


def split_sdf_file(sdf_file_stream):
    """
    Split a sdf-file containing the information for several molecules into a
    list of (cid, sdf-string).

    :param sdf_file_stream: file stream, pointing to the sdf-file. Result of 'open'

    :return: list of (cid, sdf-string)-tuples
    """
    return list(iter_sdf_file(sdf_file_stream))


def iter_sdf_file(sdf_file_stream):
    """
    Split a sdf-file containing the information for several molecules. Returns an iterator over the list
    of (cid, sdf-string)-tuples.

    :param sdf_file_stream: file stream, pointing to the sdf-file. Result of 'open'

    :yields: (cid, sdf-string)-tuple
    """
    sdf_lines = sdf_file_stream.read()
    start = 0

    while True:
        end_pos = sdf_lines.find("$$$$", start)
        if end_pos == -1:
            break

        # Extract the sdf-string and the cid
        sdf_str = sdf_lines[start:end_pos - 1]
        sdf_str = sdf_str.replace("'", "")
        cid = int(re.findall("<PUBCHEM_COMPOUND_CID>\n([0-9]+)", sdf_str)[0])

        start = end_pos + 5

        yield cid, sdf_str


def get_sdf_files_not_in_db(db_connection, sdf_fn_in_folder):
    """
    Returns the sdf_file names (full path) which are not already in the DB.

    :param db_connection: sqlite3.Connection, database connection
    :param sdf_fn_in_folder: list of string, filenames of the sdf-files.
    """
    sdf_files_in_db = db_connection.execute("SELECT filename FROM sdf_file").fetchall()
    sdf_files_in_db = [str(x[0]) for x in sdf_files_in_db]

    return sorted(list(filter(lambda x: os.path.basename(x) not in sdf_files_in_db, sdf_fn_in_folder)))


def opensdf(fn, use_gzip):
    if use_gzip:
        return gzip.open(fn, "rt")
    else:
        return open(fn, "r")


def build_db(base_dir, use_gzip, reset, db_specs):
    # Directory paths to the SDF and output DB file
    sdf_dir = os.path.join(base_dir, "sdf")
    db_dir = os.path.join(base_dir, "db")

    # Connect to the 'pubchem' database.
    conn = sqlite3.connect(os.path.join(db_dir, "pubchem.sqlite"))

    try:
        # Initialize the database
        with conn:
            initialize_db(conn, db_specs, reset=reset)

        # Get a list of all sdf-files available and reduce it to the ones still
        # needed to be processed.
        fn_patter = "*.sdf.gz" if use_gzip else "*.sdf"
        sdf_files = glob.glob(os.path.join(sdf_dir, fn_patter))
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (before filtering): %d" % n_sdf_files)

        sdf_files = get_sdf_files_not_in_db(conn, sdf_files)
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process (after filtering): %d" % n_sdf_files)

        if n_sdf_files > 0:
            # Iterate over the sdf-files and add them one by one
            for ii, sdf_fn in enumerate(sdf_files):
                print("Process sdf-file: %s (%d/%d)" % (os.path.basename(sdf_fn), ii + 1, n_sdf_files))

                # parse and insert current sdf-file
                with opensdf(sdf_fn, use_gzip) as sdf_file, conn:
                    n_inserted = insert_info_from_sdf_strings(conn, db_specs, iter_sdf_file(sdf_file))
                    print("Number of inserted compounds: %d" % n_inserted)

                    # add current sdf-file to the list of completed sdf-files
                    conn.execute("INSERT INTO sdf_file (filename, lowest_cid, highest_cid, date_added, n_compounds) "
                                 "  VALUES(?, ?, ?, DATE('now'), ?)",
                                 (os.path.basename(sdf_fn),
                                  os.path.basename(sdf_fn).split(".")[0].split("_")[1],
                                  os.path.basename(sdf_fn).split(".")[0].split("_")[2],
                                  n_inserted))

        # Create indices after all sdf-files have been parsed and imported to the database.
        for colname, specs in db_specs["columns"].items():
            if specs.get("WITH_INDEX", False):
                idx_name = "idx_%s" % colname
                with conn:
                    conn.execute("DROP INDEX IF EXISTS %s" % idx_name)
                    conn.execute("CREATE INDEX %s ON compounds(%s)" % (idx_name, colname))
                print("Create index on '%s'." % colname)

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
