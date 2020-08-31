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

import glob
import sqlite3
import os
import sys
import gzip
import traceback
import argparse

from pubchem2sqlite.utils import iter_sdf_file, initialize_db, load_db_specifications, insert_info_from_sdf_strings
from pubchem2sqlite.utils import get_sdf_files_not_in_db


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


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()

    # Required arguments
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

    # Parse arguments
    args = arg_parser.parse_args()

    # Load the DB layout
    db_specs = load_db_specifications(args.db_layout_fn)

    # Build the database
    sys.exit(build_db(args.base_dir, args.gzip, args.reset, db_specs))
