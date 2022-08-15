"""Access to AlphaFold predictions."""

import argparse
import contextlib
import ftplib
import json
import gzip
import os
import shutil
import tarfile
import tempfile
from unicodedata import decomposition

from urllib.request import urlopen
from urllib.request import urlretrieve
from urllib.request import urlcleanup


class AFDB:
    """Docstring still needs to be completed."""

    def __init__(
        self,
        server="ftp://ftp.ebi.ac.uk",
        version="latest",
        afdb=None,
        file_format="mmcif",
        verbose=True,
    ):
        """
        Docstring needs to be completed.

        A class that provides access to both bulk downloads and individual downloads of
        AF2 predictions.
        """
        # Set server and version. For now there is only one server but this leaves the possibility
        # of regionally hosted servers in the future. I don't know why someone would want to access
        # older versionfs of the database, but version option allows for that.
        self.afdb_server = server
        self.version = version
        self.file_format = file_format

        if afdb:
            self.local_afdb = afdb
        else:
            self.local_afdb = os.getcwd()

        # Enable or disable verbose
        self._verbose = verbose

    def get_metadata(self):
        """
        Docstring needs to be completed.

        Downloads and parses the AFDB metadata JSON file.
        """
        url = self.afdb_server + "/pub/databases/alphafold/download_metadata.json"
        if self._verbose:
            print("Retrieving metadata...")
        with contextlib.closing(urlopen(url)) as response:
            decoded = response.read().decode("UTF-8")
            metadata = json.loads(decoded)
        return metadata

    def retrieve_species(
        self, species, pdir=None, file_format=None, overwrite=False, decompress=True
    ):
        """
        Docstring needs to be completed.

        Provides an interface to download AF2 predictions in bulk for species that
        have tarballs provided on the FTP server.
        """
        # Set file format
        # Need to sort out keeping both
        file_format = self.file_format.lower()

        # Retrieve metadata
        afdb_metadata = self.get_metadata()

        # Retrieve list of species archives
        species_archive = {
            archive["species"]: archive["archive_name"]
            for archive in afdb_metadata
            if "species" in archive.keys()
        }

        # Check if bulk download exists for species
        if species not in species_archive.keys():
            raise KeyError(
                f"Specified species '{species}' is not available for bulk "
                "download at this time."
            )

        archive_file = species_archive[f"{species}"]

        # Where do the files get saved
        if pdir is None:
            path = self.local_afdb
        else:
            path = pdir
        species_dir = species.lower().replace(" ", "_")
        final_directory = os.path.join(path, species_dir)

        # Skip download if the directory already exists
        if not overwrite:
            if os.path.exists(final_directory):
                if self._verbose:
                    print(
                        f"AlphaFold2 bulk download directory for {species} "
                        f"already exists at: {final_directory}."
                    )
                return final_directory

        # Make directory if it doesn't exist
        if not os.path.exists(final_directory):
            os.mkdir(final_directory)

        # Retrieve the tarball
        if self._verbose:
            idx = self.afdb_server.find("://")
        if idx >= 0:
            ftp = ftplib.FTP(self.afdb_server[idx + 3 :])
        else:
            ftp = ftplib.FTP(self.afdb_server)
        ftp.login()
        ftp.cwd("/pub/databases/alphafold/latest/")

        with tempfile.NamedTemporaryFile() as tmpfile:
            if self._verbose:
                print(
                    f"Bulk downloading AlphaFold2 predicted structures for {species} in {file_format} format..."
                )
            ftp.retrbinary("RETR " + archive_file, tmpfile.write)
            tar = tarfile.open(os.path.abspath(tmpfile.name))
            tar.extractall(final_directory)

        ftp.quit()

        if decompress:
            for model in os.listdir(final_directory):
                file = os.path.join(final_directory, model)
                with gzip.open(file, "rb") as gz:
                    with open(file[:-3], "wb") as final_file:
                        shutil.copyfileobj(gz, final_file)
                    os.remove(file)

        if file_format == "mmcif":
            for model in os.listdir(final_directory):
                file = os.path.join(final_directory, model)
                if not file.endswith(".cif"):
                    os.remove(file)
        elif file_format == "pdb":
            for model in os.listdir(final_directory):
                file = os.path.join(final_directory, model)
                if not file.endswith(".pdb"):
                    os.remove(file)

        return final_directory

    def retrieve_swissprot(
        self, pdir=None, file_format=None, overwrite=False, decompress=True
    ):
        """
        Docstring needs to be completed.

        Provides an interface to download AF2 predictions in bulk for SwissProt.
        """
        # Set file format
        file_format = self.file_format.lower()

        archive = {"mmcif": "swissprot_pdb_v3.tar", "pdb": "swissprot_cif_v3.tar"}

        # Where do the files get saved
        if pdir is None:
            path = self.local_afdb
        else:
            path = pdir
        final_directory = os.path.join(path, f"swissprot_{file_format}")

        # Skip download if the directory already exists
        if not overwrite:
            if os.path.exists(final_directory):
                if self._verbose:
                    print(
                        f"AlphaFold2 bulk download directory for Swiss-Prot "
                        f"already exists at: {final_directory}."
                    )
                return final_directory

        # Make directory if it doesn't exist
        if not os.path.exists(final_directory):
            os.mkdir(final_directory)

        # Retrieve the tarball
        if self._verbose:
            idx = self.afdb_server.find("://")
        if idx >= 0:
            ftp = ftplib.FTP(self.afdb_server[idx + 3 :])
        else:
            ftp = ftplib.FTP(self.afdb_server)
        ftp.login()
        ftp.cwd("/pub/databases/alphafold/latest/")

        with tempfile.NamedTemporaryFile() as tmpfile:
            if self._verbose:
                print(
                    f"Bulk downloading AlphaFold2 predicted Swiss-Prot structures in {file_format} format..."
                )
                print(
                    "Swiss-Prot tarballs are >20 GB in size, this may take a while...\n"
                )
            ftp.retrbinary("RETR " + archive[file_format], tmpfile.write)
            tar = tarfile.open(os.path.abspath(tmpfile.name))
            tar.extractall(final_directory)

        ftp.quit()

        if decompress:
            for model in os.listdir(final_directory):
                file = os.path.join(final_directory, model)
                with gzip.open(file, "rb") as gz:
                    with open(file[:-3], "wb") as final_file:
                        shutil.copyfileobj(gz, final_file)
                    os.remove(file)

        return final_directory

    def retrieve_mane(self, pdir=None, overwrite=False, decompress=True):
        """
        Docstring needs to be completed.

        Provides an interface to download MANE Select dataset.
        """
        # Where do the files get saved
        if pdir is None:
            path = self.local_afdb
        else:
            path = pdir
        final_directory = os.path.join(path, "mane_select")

        # Skip download if the directory already exists
        if not overwrite:
            if os.path.exists(final_directory):
                if self._verbose:
                    print(
                        f"AlphaFold2 bulk download directory for MANE Select "
                        f"already exists at: {final_directory}."
                    )
                return final_directory

        # Make directory if it doesn't exist
        if not os.path.exists(final_directory):
            os.mkdir(final_directory)

        # Retrieve the tarball
        if self._verbose:
            idx = self.afdb_server.find("://")
        if idx >= 0:
            ftp = ftplib.FTP(self.afdb_server[idx + 3 :])
        else:
            ftp = ftplib.FTP(self.afdb_server)
        ftp.login()
        ftp.cwd("/pub/databases/alphafold/latest/")

        with tempfile.NamedTemporaryFile() as tmpfile:
            if self._verbose:
                print("Bulk downloading AlphaFold2 predicted MANE Select structures...")
            ftp.retrbinary("RETR mane_overlap_v3.tar", tmpfile.write)
            tar = tarfile.open(os.path.abspath(tmpfile.name))
            tar.extractall(final_directory)

        ftp.quit()

        if decompress:
            for model in os.listdir(final_directory):
                file = os.path.join(final_directory, model)
                with gzip.open(file, "rb") as gz:
                    with open(file[:-3], "wb") as final_file:
                        shutil.copyfileobj(gz, final_file)
                    os.remove(file)
        return final_directory

    def retrieve_prediction(
        self, uniprot_id, pdir=None, file_format=None, overwrite=False
    ):
        """
        Docstring needs to be completed.

        Provides an interface to download individual AF2 predictions by
        UniProt ID.
        """
        # Set file format
        file_format = self.file_format.lower()

        # Ensure ID is uppercase
        uniprot_id = uniprot_id.upper()

        archive = {
            "mmcif": f"AF-{uniprot_id}-F1-model_v2.cif",
            "pdb": f"AF-{uniprot_id}-F1-model_v2.pdb",
        }

        # Where do the files get saved
        if pdir is None:
            path = self.local_afdb
        else:
            path = pdir
        final_directory = os.path.join(path, "individual_files")
        if not os.path.isdir(final_directory):
            os.mkdir(final_directory)
        final_file = os.path.join(final_directory, archive[file_format])

        # Skip download if the file already exists
        if not overwrite:
            if os.path.exists(final_file):
                if self._verbose:
                    print(f"AlphaFold2 prediction exists: {final_file}")
                return final_file

        url = f"https://alphafold.ebi.ac.uk/files/{archive[file_format]}"

        if self._verbose:
            print(f"Downloading AlphaFold2 prediction for UniProt ID {uniprot_id}...")
        try:
            urlcleanup()
            urlretrieve(url, final_file)
        except OSError:
            print(f"Download failed! Is {uniprot_id} the correct UniProt ID?")

        return final_file


if __name__ == "__main__":

    # Instantiate argument parser/subparser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description="Placeholder"
    )

    # Set AFDB location
    parser.add_argument(
        "afdb_path",
        nargs="?",
        default=os.getcwd(),
        help="local path to use, defaults to the current working directory",
    )

    # Add mutually exclusive group for modes
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-i",
        "--identifiers",
        metavar="uniprot id(s)",
        help="retrieve predictions(s) - accepts a single uniprot id or a comma separated list of uniprot ids",
    )
    mode.add_argument(
        "-m",
        "--mane",
        action="store_true",
        help="bulk download AlphaFold2 predictions for the MANE Select dataset",
    )
    mode.add_argument(
        "-s",
        "--species",
        metavar='"Species name"',
        nargs=1,
        help="bulk download AlphaFold2 predictions for provided species name (in quotation marks)",
    )
    mode.add_argument(
        "-w",
        "--swissprot",
        metavar="format",
        nargs=1,
        help="bulk download AlphaFold2 predictions for SwissProt in a specified format",
    )

    # Add optional format argument that select choice
    options = parser.add_argument_group()
    options.add_argument(
        "-d",
        "--decompress",
        choices=(True, False),
        default=True,
        help="decompress bulk downloaded predictions, defaults to True",
    )
    options.add_argument(
        "-f",
        "--format",
        choices=("mmcif", "pdb"),
        default="mmcif",
        help="format in which to download predictions, defaults to mmcif",
    )
    options.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="overwrite existing prediction files",
    )

    # Parse arguments
    args = parser.parse_args()

    # Set parameters
    afdb_path = args.afdb_path
    decompress = args.decompress
    file_format = args.format
    overwrite = args.overwrite

    # Instantiate AFDB object
    af = AFDB(afdb=afdb_path)

    # Run
    if args.identifiers:
        uniprot_ids = args.identifiers.split(",")
        for id in uniprot_ids:
            af.retrieve_prediction(
                id, pdir=afdb_path, file_format=file_format, overwrite=overwrite
            )
    elif args.mane:
        af.retrieve_mane(pdir=afdb_path, overwrite=overwrite, decompress=decompress)
    elif args.species:
        species = args.species[0]
        af.retrieve_species(
            species,
            pdir=afdb_path,
            file_format=file_format,
            overwrite=overwrite,
            decompress=decompress,
        )
    elif args.swissprot:
        af.retrieve_swissprot(
            pdir=afdb_path,
            file_format=file_format,
            overwrite=overwrite,
            decompress=decompress,
        )
