# Build an SQLite Database from Pubchem

This library can be used to build an [SQLite](https://www.sqlite.org/index.html) database (DB) containing all compounds in [PubChem](https://pubchem.ncbi.nlm.nih.gov/). A simple [json-file](default_db_layout.json) is used to specify the layout of the SQLite DB allowing to controll which PubChem [SD-tags](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf) should be included.

## Features: 

- Construct a local copy of the PubChem *compound* database from publicly available *sdf*-files.
- Local database is SQLite database
- Database layout and the PubChem information to extract can be specified. 

## Installation 
 
The library libraries only requirement is Python (tested on 3.5 - 3.9), no external libraries need to be installed. Simply run:
```bash
git clone https://github.com/bachi55/local_pubchem_db
cd local_pubchem_db
pip install .
```

## Usage

### (1) Download PubChem SDF-files

PubChem allows to [download a snapshot](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/) of their current database (~75GB). You can use the following command to retrieve the sdf-files: 
```bash
wget -e robots=off --recursive --no-parent https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/
mkdir -p pubchem/sdf
mv ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/*.sdf.gz pubchem/sdf
rm -r ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF
```

### (2) Specifiy the DB Layout

You can specify the DB layout using a simple json-file (see also this [example](default_db_layout.json)):
```json
{
    "columns": {
        "InChIKey": {
            "SD_TAG": ["PUBCHEM_IUPAC_INCHIKEY"],
            "DTYPE": "varchar",
            "NOT_NULL": true,
            "PRIMARY_KEY": true
        },
        "InChI": {
            "SD_TAG": ["PUBCHEM_IUPAC_INCHI"],
            "DTYPE": "varchar",
            "NOT_NULL": true
        },
        "xlogp3": {
            "SD_TAG": ["PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"],
            "DTYPE": "real",
            "NOT_NULL": false
        },
        "cid": {
            "SD_TAG": ["PUBCHEM_COMPOUND_CID"],
            "DTYPE": "integer",
            "NOT_NULL": true,
        }
    }
}
```
The json-file contains consists of nested dictionaries. The one containing the DB layout is accessed via ```columns```. This dictionary contains a key for each column name in the [*compounds* table](https://github.com/bachi55/local_pubchem_db/blob/bd339a19ffd8b442eda54f0b8684270eabf4c357/pubchem2sqlite/utils.py#L162) and its properites:

| Key | Description |
| --- | --- |
| [SD_TAG](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf) | List of SD-Tags correponding to the column. For example ```[PUBCHEM_IUPAC_INCHI]``` refering to the InChI of a compound, or, ```["PUBCHEM_XLOGP3", "PUBCHEM_XLOGP3_AA"]``` to its predicted XLogP3 value (here more than one tag can occure in the sdf-files). | 
| [DTYPE](https://github.com/bachi55/local_pubchem_db/blob/bd339a19ffd8b442eda54f0b8684270eabf4c357/pubchem2sqlite/utils.py#L9) | Type of the information behind the sd-tag. Must be a valid SQLite type. |
| NOT_NULL | If true, the column cannot contain null entries. If a compound does not have the requested property and the corresponding column should not be null, *it will not be added* to tha DB. | 
| PRIMARY_KEY | If true, the column is used as primary key for the *compounds* table. Currently, only one column can be specified as primary key. | 
| WITH_INDEX | If true, an index is created for the column. This allows faster queries with constraints on this column. | 

### (3) Build the Database

The DB is build using: 
```bash
cd pubchem
mkdir db
# Create your db-layout
touch db/db_layout.json
# Edit it ... 
# vim db/db_layout.json
python /path/to/local_pubchem_db/build_pubchem_db.py pubchem --gzip --db_layout_fn=db/db_layout.json
```
Please note: PubChem contains *many* compound and the resulting SQLite file can be very large. Especially when the indices are created sqlite might require twice the memory of the DB for some time. Ensure your drive has enough memory. If you face the problem, that your systems temp-directory runs out of memory, when creating the indices, have a look [here](https://sqlite.org/tempfiles.html#temporary_file_storage_locations). You can temporarily change sqlite's temp-directory by setting ```SQLITE_TMPDIR```:
```bash
SQLITE_TMPDIR=/my/large/disk/temp python /path/to/local_pubchem_db/build_pubchem_db.py pubchem --gzip --db_layout_fn=db/db_layout.json
```
## Version History

#### 0.2:
- Add flexibility for the DB-layout using a json-file
- Add tests
- Publish code as package
- Change license to MIT

#### 0.1:
- Initial Version
