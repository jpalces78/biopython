# Copyright 2016 by Jacek Smietanski.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Testing download access to AlphaFold Protein Structure Database"""

import contextlib
import os
import shutil
import tempfile
import unittest

from Bio.PDB.AFDB import AFDB

class TestAFDB(unittest.TestCase):
    """Test methods of the AFDB class."""

    @contextlib.contextmanager
    def make_temp_directory(self, directory):
        temp_dir = tempfile.mkdtemp(dir=directory)
        try:
            yield temp_dir
        finally:
            shutil.rmtree(temp_dir)

    def check_retrieve_species(self, species, file_format, expected_num_files):
        with self.make_temp_directory(os.getcwd()) as tmp:
            af = AFDB()
            pdir = os.path.join(tmp, "test_afdb")
            os.mkdir(pdir)
            af.retrieve_species(species, pdir=pdir, file_format=file_format)
            num_files = sum(len(files) for _, _, files in os.walk(pdir))
            self.assertEqual(num_files, expected_num_files)

    def test_retrieve_species(self):
        """Test retrieving AF2 predictions for a species."""
        species = "Oryza sativa"
        self.check_retrieve_species(species, "mmcif", 10)  # Change the expected number of files accordingly

    # Add more test methods for other functionalities of the AFDB class


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
