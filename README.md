# Build an SQLite Database from Pubchem

This library can be used to build an [SQLite](https://www.sqlite.org/index.html) database (DB) containing all compounds in [PubChem](https://pubchem.ncbi.nlm.nih.gov/). A simple [json-file](default_db_layout.json) is used to specify the layout of the SQLite DB allowing to controll which PubChem [SD-tags](https://ftp.ncbi.nlm.nih.gov/pubchem/specifications/pubchem_sdtags.pdf) should be included.

## Features: 

- Construct a local copy of the PubChem *compound* database from publicly available *sdf*-files.
- Local database is SQLite database
- Database layout and the PubChem information to extract can be specified. 

## Installation 
 
The library libraries only requirement is Python (TODO: add version), no external libraries need to be installed. Simply run:
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
