from setuptools import setup, find_packages

setup(
    name="pubchem2sqlite",
    version="0.2.0",
    license="MIT",
    packages=find_packages(exclude=["tests", "build_pubchem_db.py"]),

    # Minimum requirements the package was tested with
    install_requires=[
        "setuptools>=46.1"
    ],

    # Metadata
    author="Eric Bach",
    author_email="eric.bach@aalto.fi",
    description="Simple library to construct a local copy of the PubChem compound database using SQLite.",
    url="https://github.com/bachi55/local_pubchem_db",
)