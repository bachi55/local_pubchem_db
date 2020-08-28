import unittest

from misc import get_column_stmt


class ParsingJSONDBSpecs(unittest.TestCase):
    def test_get_column_stmt(self):
        specs = {
            "columns": {
                "MASS": {
                    "DTYPE": "float", "NOT_NULL": False
                },
                "INCHI": {
                    "DTYPE": "string", "NOT_NULL": True
                },
                "CID": {
                    "DTYPE": "integer", "PRIMARY_KEY": True
                }
            }
        }

        stmt = get_column_stmt(specs["columns"])
        self.assertEqual("MASS float,INCHI string not null,CID integer not null primary key", stmt)


if __name__ == '__main__':
    unittest.main()
