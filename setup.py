from setuptools import setup

setup(
  name = "swiss",
  version = "1.1.1",
  author = "Ryan Welch",
  author_email = "welchr@umich.edu",

  # Packages
  packages = ["swiss","swiss/conf"],

  package_data = {
    "swiss/conf": ["default.yaml"]
  },

  # Details
  url = "https://github.com/welchr/Swiss",
  license = "LICENSE.txt",
  description = "Software to help identify overlap between association scan results and GWAS hit catalogs.",

  long_description = open("README.md").read(),

  install_requires = open("requirements.txt").read().split(),

  scripts = [
    "bin/swiss",
    "bin/swiss-create-data"
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

