import sqlite3
import os
import sys
import glob

from timeit import default_timer as timer
from misc import split_sdf_file, extract_single_info_from_sdf

if __name__ == "__main__":

    if len(sys.argv) < 5:
        print("Usage: python %s <PUBCHEM_DB_BASEDIR> <NEW_COLUMN_NAME> <NEW_COLUMN_TYPE> <SDF_IDENTIFIER>" %
              os.path.basename(sys.argv[0]))
        print("Example:")
        print("\t python %s /run/media/bach/Intenso/data/pubchem/ can_smiles varchar PUBCHEM_OPENEYE_CAN_SMILES" %
              os.path.basename(sys.argv[0]))
        exit(1)

    basedir = sys.argv[1]
    new_column_name = sys.argv[2]
    new_column_type = sys.argv[3]
    sdf_indentifier = sys.argv[4]

    if new_column_type in ["integer", "real"]:
        raise NotImplemented("Currently only 'varchar' columns can be added.")

    sdf_dir = basedir + "/sdf/"
    db_dir = basedir + "/db/"

    # Connect to the 'pubchem' database.
    conn = sqlite3.connect(db_dir + "/pubchem.db", isolation_level=None)
    try:
        # Add new column
        with conn:
            # FIXME: Its not really clear, why we cannot use the '?' placeholders here.
            conn.execute("ALTER TABLE info ADD COLUMN %s %s" % (new_column_name, new_column_type))

        # Get a list of all sdf-files available and reduce it to the ones still
        # needed to be processed.
        sdf_files = glob.glob(sdf_dir + "*.sdf")
        n_sdf_files = len(sdf_files)
        print("Sdf-files to process: %d" % n_sdf_files)

        if n_sdf_files > 0:
            # Iterate over the sdf-files and add them one by one
            for ii, sdf_fn in enumerate(sdf_files):
                print("Process sdf-file: %d/%d" % (ii + 1, n_sdf_files))

                # parse current sdf-file
                with open(sdf_fn, "r") as sdf_file:
                    l_cid_sdf = split_sdf_file(sdf_file)

                # insert information
                with conn:
                    start = timer()

                    conn.execute("BEGIN")

                    for cid, sdf in l_cid_sdf:
                        # extract desired information from sdf
                        info = extract_single_info_from_sdf(sdf, field_identifier=sdf_indentifier,
                                                            dtype=new_column_type, return_cid=True)

                        assert (cid == info["PUBCHEM_COMPOUND_CID"])

                        # alter entry
                        conn.execute("UPDATE info SET %s = '%s' WHERE cid is %d" % (new_column_name,
                                                                                    info[sdf_indentifier], cid))

                    conn.execute("COMMIT")

                    end = timer()
                    print("Extraction and update of the information took %.3fsec" % (end - start))

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
