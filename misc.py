import re
import json


def load_db_specifications(fn):
    """
    Loads a database layout specified in a json-file into a dictionary.

    :param fn: string, path to the json-file specifying the database layout. See 'default_db_layout.json' as an example.

    :return: dictionary
    """
    with open(fn, "r") as json_file:
        specs = json.loads(json_file.read())
    return specs


def get_column_stmt(column_specs):
    stmt_columns = []

    has_primary_key = False  # Primary keys must (if) be defined on a single column.

    for name, spec in column_specs.items():
        new_col = [name, spec["DTYPE"]]

        if spec.get("NOT_NULL", True):
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
    if not reset:
        return

    # Table to keep track about the sdf-chucks (files) included
    db_connection.execute("DROP TABLE IF EXISTS sdf_file")
    db_connection.execute("CREATE TABLE sdf_file("
                          "filename          VARCHAR PRIMARY KEY,"
                          "lowest_cid        INTEGER,"
                          "highest_cid       INTEGER,"
                          "last_modification REAL)")

    # Table containing the compounds
    db_connection.execute("DROP TABLE IF EXISTS compounds")
    db_connection.execute("CREATE TABLE compounds(%s)" % get_column_stmt(specs["columns"]))


def split_sdf_file(sdf_file_stream):
    """
    Split a sdf-file containing the information for several molecules into a
    list of (cid, sdf-string).

    :param sdf_file_stream: file stream, pointing to the sdf-file. Result of 'open'

    :return: list of (cid, sdf-string)-tuples
    """
    sdf_lines = sdf_file_stream.read()
    start = 0
    sdfs = []

    while True:
        end_pos = sdf_lines.find("$$$$", start)
        if end_pos == -1:
            break

        # Extract the sdf-string and the cid
        sdf_str = sdf_lines[start:end_pos - 1]
        sdf_str = sdf_str.replace("'", "")
        cid = int(re.findall("<PUBCHEM_COMPOUND_CID>\n([0-9]+)", sdf_str)[0])

        sdfs.append((cid, sdf_str))

        start = end_pos + 5

    return sdfs


def extract_single_info_from_sdf(sdf, field_identifier, dtype, return_cid=True):
    """
    Extract the information for a given molecules from its sdf-string.

    :param sdf: string, containing the molecule's sdf
    :param field_identifier: string, which field of the sdf should be extracted
    :param dtype: string, of which type is the desired information. Can be: "integer", "real" or "varchar"
    :param return_cid: boolean, indicating whether the pubchem id should be returned as well, e.g. for consistency
        checks.

    return: dictionary, containing the extracted information: {field_identifier: value}
    """
    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0

    info = {field_identifier: None}
    if return_cid:
        info["PUBCHEM_COMPOUND_CID"] = None

    found = {key: False for key in info.keys()}

    while not all(found.values()) and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i += 1
            continue

        if return_cid and not found["PUBCHEM_COMPOUND_CID"] and lines[i].find("PUBCHEM_COMPOUND_CID") != -1:
            found["PUBCHEM_COMPOUND_CID"] = True
            info["PUBCHEM_COMPOUND_CID"] = int(lines[i + 1].strip())
            i += 3
            continue

        if not found[field_identifier] and lines[i].find(field_identifier) != -1:
            found[field_identifier] = True
            if dtype == "integer":
                info[field_identifier] = int(lines[i + 1].strip())
            elif dtype == "real":
                info[field_identifier] = float(lines[i + 1].strip())
            elif dtype == "varchar":
                info[field_identifier] = lines[i + 1].strip()
            else:
                raise ValueError("Invalid dtype: %s." % dtype)

            i += 3
            continue

        i += 1

    return info
