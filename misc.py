import re
import json

from collections import OrderedDict
from timeit import default_timer as timer


def _as_dtype(val, dtype):
    """
    Convert a value represented using a string variable to the desired output data-type.

    :param val:
    :param dtype:
    :return:
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

    lines = sdf.split("\n")
    n_line = len(lines)
    i = 0

    while any([v is None for v in infos.values()]) and i < (n_line - 1):
        # skip empty lines
        if not len(lines[i]):
            i += 1
            continue

        # Check whether the current line contains any required information
        found_info_in_line = False
        for info in infos:
            if infos[info] is None:
                dtype = db_specs["columns"][info]["DTYPE"].lower()  # type of the information associated with the SD-tag

                for sd_tag in db_specs["columns"][info]["SD_TAG"]:
                    if lines[i].find(sd_tag) != -1:
                        found_info_in_line = True
                        infos[info] = _as_dtype(lines[i + 1].strip(), dtype)

                        i += 2
                        break

            if found_info_in_line:
                break

        i += 1

    return infos


def extract_info_from_sdf_OLD(sdf):
    """
    Extract the information for a given molecules from its sdf-string.

    :param sdf: string, containing the molecule's sdf

    :return: dictionary, containing the extracted information.
    """

    info = {
        "PUBCHEM_COMPOUND_CID": None,  # integer
        "PUBCHEM_IUPAC_NAME": None,  # string
        "PUBCHEM_IUPAC_INCHI": None,  # string
        "PUBCHEM_IUPAC_INCHIKEY": None,  # string
        "INCHIKEY1": None,  # string
        "PUBCHEM_XLOGP3*": None,  # float
        "PUBCHEM_MONOISOTOPIC_WEIGHT": None,  # float
        "PUBCHEM_MOLECULAR_FORMULA": None,  # string
        "PUBCHEM_OPENEYE_ISO_SMILES": None,  # string
        "PUBCHEM_OPENEYE_CAN_SMILES": None,  # string
        "PUBCHEM_ATOM_DEF_STEREO_COUNT": None,  # integer
        "PUBCHEM_ATOM_UDEF_STEREO_COUNT": None,  # integer
        "PUBCHEM_BOND_DEF_STEREO_COUNT": None,  # integer
        "PUBCHEM_BOND_UDEF_STEREO_COUNT": None  # integer
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
            info["PUBCHEM_COMPOUND_CID"] = int(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_NAME"] and lines[i].find("PUBCHEM_IUPAC_NAME") != -1:
            found["PUBCHEM_IUPAC_NAME"] = True
            info["PUBCHEM_IUPAC_NAME"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_INCHI"] and lines[i].find("PUBCHEM_IUPAC_INCHI") != -1:
            found["PUBCHEM_IUPAC_INCHI"] = True
            info["PUBCHEM_IUPAC_INCHI"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_IUPAC_INCHIKEY"] and lines[i].find("PUBCHEM_IUPAC_INCHIKEY") != -1:
            found["PUBCHEM_IUPAC_INCHIKEY"] = True
            found["INCHIKEY1"] = True
            info["PUBCHEM_IUPAC_INCHIKEY"] = lines[i + 1].strip()
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
            info["PUBCHEM_XLOGP3*"] = float(lines[i + 1].strip())
            i = i + 2
            continue
        if not found["PUBCHEM_XLOGP3*"] and lines[i].find("PUBCHEM_XLOGP3_AA") != -1:
            found["PUBCHEM_XLOGP3*"] = True
            info["PUBCHEM_XLOGP3*"] = float(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_MONOISOTOPIC_WEIGHT"] and lines[i].find("PUBCHEM_MONOISOTOPIC_WEIGHT") != -1:
            found["PUBCHEM_MONOISOTOPIC_WEIGHT"] = True
            info["PUBCHEM_MONOISOTOPIC_WEIGHT"] = float(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_MOLECULAR_FORMULA"] and lines[i].find("PUBCHEM_MOLECULAR_FORMULA") != -1:
            found["PUBCHEM_MOLECULAR_FORMULA"] = True
            info["PUBCHEM_MOLECULAR_FORMULA"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_OPENEYE_ISO_SMILES"] and lines[i].find("PUBCHEM_OPENEYE_ISO_SMILES") != -1:
            found["PUBCHEM_OPENEYE_ISO_SMILES"] = True
            info["PUBCHEM_OPENEYE_ISO_SMILES"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_OPENEYE_CAN_SMILES"] and lines[i].find("PUBCHEM_OPENEYE_CAN_SMILES") != -1:
            found["PUBCHEM_OPENEYE_CAN_SMILES"] = True
            info["PUBCHEM_OPENEYE_CAN_SMILES"] = lines[i + 1].strip()
            i = i + 2
            continue

        if not found["PUBCHEM_ATOM_DEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_ATOM_DEF_STEREO_COUNT") != -1:
            found["PUBCHEM_ATOM_DEF_STEREO_COUNT"] = True
            info["PUBCHEM_ATOM_DEF_STEREO_COUNT"] = int(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_ATOM_UDEF_STEREO_COUNT") != -1:
            found["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] = True
            info["PUBCHEM_ATOM_UDEF_STEREO_COUNT"] = int(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_BOND_DEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_BOND_DEF_STEREO_COUNT") != -1:
            found["PUBCHEM_BOND_DEF_STEREO_COUNT"] = True
            info["PUBCHEM_BOND_DEF_STEREO_COUNT"] = int(lines[i + 1].strip())
            i = i + 2
            continue

        if not found["PUBCHEM_BOND_UDEF_STEREO_COUNT"] and lines[i].find("PUBCHEM_BOND_UDEF_STEREO_COUNT") != -1:
            found["PUBCHEM_BOND_UDEF_STEREO_COUNT"] = True
            info["PUBCHEM_BOND_UDEF_STEREO_COUNT"] = int(lines[i + 1].strip())
            i = i + 2
            continue

        i = i + 1

    return info


def insert_info_from_sdf_strings(db_connection, db_specs, iter_cid_sdf):
    """
    Insert the basic information of each molecule in the DB (represented using
    its sdf-string) and insert those information into the 'info' table.

    :param db_connection: sqlite3.Connection, database connection

    :param db_specs: dictionary, containing the database specifications (see 'default_db_layout.json')

    :param iter_cid_sdf: list of (cid, sdf-string)-tuples
    """
    start = timer()

    n_cols = len(db_specs["columns"])

    new_rows = []

    for cid, sdf in iter_cid_sdf:
        mol_sdf = extract_info_from_sdf(sdf, db_specs)

        # We skip PubChem entries that do not provide all information requested to be NOT NULL
        for k, v in db_specs["columns"].items():
            if v.get("NOT_NULL", False) and mol_sdf[k] is None:
                continue

        # db_connection.execute("INSERT INTO compounds (%s) VALUES (%s)"
        #                       % (",".join(db_specs["columns"].keys()), ",".join(["?"] * n_cols)),
        #                       tuple(v for v in mol_sdf.values()))

        new_rows.append(tuple(v for v in mol_sdf.values()))

    n_cols = len(db_specs["columns"])
    db_connection.executemany("INSERT INTO compounds (%s) VALUES (%s)"
                              % (",".join(db_specs["columns"].keys()), ",".join(["?"] * n_cols)),
                              new_rows)

    end = timer()
    print("Extraction and insertion of the information took %.3fsec" % (end - start))


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
                          "highest_cid       INTEGER)")

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
