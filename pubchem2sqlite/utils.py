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
import json
import os

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

    # type of the information associated with the SD-tag
    dtypes = {k: specs["DTYPE"].lower() for k, specs in db_specs["columns"].items()}
    sdtags = {k: ["> <%s>" % v for v in specs["SD_TAG"]] for k, specs in db_specs["columns"].items()}

    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0

    while missing_infos and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i += 1
            continue

        # Check whether the current line contains any required information
        found_info_in_line = False
        for info in missing_infos:
            for sd_tag in sdtags[info]:
                if sd_tag == lines[i]:
                    infos[info] = _as_dtype(lines[i + 1], dtypes[info])

                    found_info_in_line = True
                    missing_infos.remove(info)

                    i += 2
                    break

            if found_info_in_line:
                break

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
