import glob
import sqlite3
import os
import sys
import gzip
import traceback
import argparse

from datetime import datetime

from misc import iter_sdf_file, initialize_db, load_db_specifications, insert_info_from_sdf_strings
from misc import get_sdf_files_not_in_db

"""
This script can be used to parse a set of SDF files and add (some of) the
content to a SQLite database.

This code is based on Huibin's implementation.
"""


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
                    insert_info_from_sdf_strings(conn, db_specs, iter_sdf_file(sdf_file))

                    # add current sdf-file to the list of completed sdf-files
                    conn.execute("INSERT INTO sdf_file (filename, lowest_cid, highest_cid) "
                                 "  VALUES(?,?,?)",
                                 (os.path.basename(sdf_fn),
                                  os.path.basename(sdf_fn).split(".")[0].split("_")[1],
                                  os.path.basename(sdf_fn).split(".")[0].split("_")[2]))

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
    sys.exit(main())
