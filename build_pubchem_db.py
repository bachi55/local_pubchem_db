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

import sys
import argparse

from pubchem2sqlite import load_db_specifications, build_db


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
