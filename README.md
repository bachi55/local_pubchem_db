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
```
