from distutils.core import setup
from glob import glob
import os

def datafiles(dirpath):
  for d, ds, fs in os.walk(dirpath):
    if "genome" in d:
      continue

    for f in fs:
      if f.endswith(".py") or f.endswith(".pyc"):
        continue

      if ".vcf.gz" in f and not f.startswith("trim"):
        continue

      if "fusion" in f:
        continue

      if "got2d" in f.lower():
        continue

      yield os.path.join(d,f).replace("swiss/","")

setup(
  name = "swiss",
  version = "1.0.0",
  author = "Ryan Welch",
  author_email = "welchr@umich.edu",

  # Packages
  packages = ["swiss","swiss/conf"],

  package_data = {
    "swiss": [f for f in datafiles("swiss/data")],
    "swiss/conf": ["default.yaml"]
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
    "bx-python",
    "pyyaml",
    "appdirs"
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
    'Operating System :: POSIX :: Linux',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
  ]
)

