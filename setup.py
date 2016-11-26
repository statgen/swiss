from distutils.core import setup
from glob import glob
import os

def subfiles(dirpath):
  for d, ds, fs in os.walk(dirpath):
    for f in fs:
      if f.endswith(".vcf.gz") and not f.startswith("trim"):
        continue

      if "fusion" in f:
        continue

      yield os.path.join(d,f).replace("swiss/","")

setup(
  name = "swiss",
  version = "1.0.0",
  author = "Ryan Welch",
  author_email = "welchr@umich.edu",

  # Packages
  packages = ["swiss"],

  package_data = {
    "swiss": [f for f in subfiles("swiss/data")]
  },

  # Details
  url = "https://github.com/welchr/Swiss",
  license = "LICENSE.txt",
  description = "Software to help identify overlap between association scan results and GWAS hit catalogs.",

  long_description=open("README.md").read(),

  install_requires = [
    "termcolor",
    "pandas",
    "pysam",
    "bx-python"
  ],

  scripts = [
    "bin/swiss",
  ],

  classifiers = [
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 2.7',
    #'Programming Language :: Python :: 3.5',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ]
)

