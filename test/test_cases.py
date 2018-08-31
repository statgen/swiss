# -*- coding: utf-8 -*-
"""
This file contains basic test cases for swiss.
"""

from __future__ import print_function
import inspect
import csv
import os
from os import path
from math import isnan
import pytest

try:
  from swiss.main import main as swiss_main
except:
  # We're running tests for the binary/compiled version of swiss.
  from subprocess import check_call
  def swiss_main(arg_string):
    print(arg_string)
    check_call(arg_string, shell=True)

  os.environ["PATH"] = os.environ["SWISS_BINARY_PATH"] + os.pathsep + os.environ["PATH"]

def assert_files_exist(*files):
  """
  Assert that all files passed to this function exist.
  """

  for filepath in files:
    assert filepath.exists()

def test_ffa(tmpdir):
  frame = inspect.currentframe()
  func = frame.f_code.co_name
  prefix = tmpdir.join(func)
  whereami = path.join(path.dirname(__file__))

  data = path.join(whereami, "data/test_hg19.gz")
  gwascat = path.join(whereami, "data/gwascat_ebi_GRCh37p13.tab")
  args = "swiss --assoc {data} --multi-assoc --trait SM --build hg19 " \
         "--ld-clump-source 1000G_2014-11_EUR --ld-gwas-source 1000G_2014-11_EUR --gwas-cat {gwascat} " \
         "--ld-clump --clump-p 1e-10 --out {prefix}".format(data=data, gwascat=gwascat, prefix=prefix)
  swiss_main(args)

  file_clump = prefix + ".SM.clump"
  file_neargwas = prefix + ".SM.near-gwas.tab"
  file_ldgwas = prefix + ".SM.ld-gwas.tab"

  assert_files_exist(file_clump, file_neargwas, file_ldgwas)

  clump_data = dict()
  with open(str(file_clump)) as infile:
    reader = csv.DictReader(infile, delimiter="\t")
    for line in reader:
      vid = line["MARKER_ID"]
      clump_data[vid] = line

  assert clump_data["4:73303394_G/A"]["ld_with"] == ""
  pytest.approx(float(clump_data["4:73303394_G/A"]["PVALUE"]), 6.04e-11)

  for v in clump_data["4:74033564_G/C"]["ld_with_values"].split(","):
    pytest.approx(float(v), 0.50)

  expected_ld = [1.00, 0.63, 0.65, 0.64, 0.64]
  for exp, obs in zip(expected_ld, clump_data["19:45413233_G/T"]["ld_with_values"].split(",")):
    pytest.approx(exp, float(obs))

def test_hg38(tmpdir):
  """
  Run an hg38 result file through swiss and check for expected output.
  """

  frame = inspect.currentframe()
  func = frame.f_code.co_name
  prefix = tmpdir.join(func)
  whereami = path.join(path.dirname(__file__))

  data = path.join(whereami, "data/test_hg38.gz")
  gwascat = path.join(whereami, "data/gwascat_ebi_GRCh38p7.tab")
  args = "swiss --assoc {data} --gwas-cat {gwascat} --variant-col VARIANT " \
         "--chrom-col CHROM --pos-col POS --trait BMI --build hg38 " \
         "--ld-clump-source 1000G_2014-11_EUR --ld-gwas-source 1000G_2014-11_EUR " \
         "--ld-clump --clump-p 5e-08 --out {prefix}".format(data=data, prefix=prefix, gwascat=gwascat)
  swiss_main(args)

  file_clump = prefix + ".clump"
  file_neargwas = prefix + ".near-gwas.tab"
  file_ldgwas = prefix + ".ld-gwas.tab"

  assert_files_exist(file_clump, file_neargwas, file_ldgwas)

  with open(str(file_clump)) as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    i = None
    for i, row in enumerate(reader):
      # There should only be one row in this result set.
      assert i == 0

      # Check that p-value matches expected value.
      pvalue = float(row["PVALUE"])
      pytest.approx(pvalue, 8.14562e-09)

      ld_values = [float(x) for x in row["ld_with_values"].split(",")]
      expected_ld_values = [0.95, 1.00, 1.00, 1.00, 1.00, 1.00, 0.96, 0.83]
      for pair in zip(ld_values, expected_ld_values):
        pytest.approx(*pair)

def test_distclump(tmpdir):
  """
  Run a test case using distance based clumping.
  """

  frame = inspect.currentframe()
  func = frame.f_code.co_name
  prefix = tmpdir.join(func)
  whereami = path.join(path.dirname(__file__))

  data = path.join(whereami, "data/test_hg19.gz")
  gwascat = path.join(whereami, "data/gwascat_ebi_GRCh37p13.tab")
  args = "swiss --assoc {data} --multi-assoc --trait SM --build hg19 " \
         "--dist-clump --ld-gwas-source 1000G_2014-11_EUR --gwas-cat {gwascat} " \
         "--clump-p 1e-10 --out {prefix}".format(data=data, prefix=prefix, gwascat=gwascat)
  swiss_main(args)

  file_clump = prefix + ".SM.clump"
  file_neargwas = prefix + ".SM.near-gwas.tab"
  file_ldgwas = prefix + ".SM.ld-gwas.tab"

  assert_files_exist(file_clump, file_neargwas, file_ldgwas)

  with open(str(file_clump)) as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    i = None
    variants = []
    pvalues = []
    for i, row in enumerate(reader):
      variants.append(row["variant"])
      pvalues.append(float(row["PVALUE"]))

    expected_variants = "4:73491622_G/A 4:73768622_G/T 4:74033564_G/C 4:74512280_G/A 19:45413233_G/T".split()
    expected_pvalues = [float(x) for x in "2.0899999999999998e-12 1.01e-13 4.22e-14 2.13e-12 2.6800000000000002e-17".split()]

    for pair in zip(variants, expected_variants):
      assert pair[0] == pair[1]

    for pair in zip(pvalues, expected_pvalues):
      pytest.approx(*pair)

  with open(str(file_neargwas)) as infile:
    reader = csv.DictReader(infile, delimiter="\t")

    i = None
    gwas_variants = []
    distances = []
    for i, row in enumerate(reader):
      if i > 5:
        break

      gwas_variants.append(row["GWAS_EPACTS"])
      distances.append(int(row["ASSOC_GWAS_DIST"]))

    expected_gwas_variants = "4:73322565_T/C 4:73645351_G/A 4:73515313_T/C 4:73696709_A/T 4:73515825_C/T".split()
    expected_distances = [int(x) for x in "169057 153729 23691 205087 24203".split()]

    for pair in zip(gwas_variants, expected_gwas_variants):
      assert pair[0] == pair[1]

    for pair in zip(distances, expected_distances):
      assert pair[0] == pair[1]

def test_known_gwas_is_tophit(tmpdir):
  frame = inspect.currentframe()
  func = frame.f_code.co_name
  prefix = tmpdir.join(func)
  whereami = path.join(path.dirname(__file__))

  data = path.join(whereami, "data/top_hit_is_gwas.tab")
  gwascat = path.join(whereami, "data/gwascat_ebi_GRCh37p13.tab")
  args = "swiss --assoc {data} --build hg19 --ld-clump-source 1000G_2014-11_ALL " \
         "--ld-gwas-source 1000G_2014-11_ALL --gwas-cat {gwascat} " \
         "--variant-col EPACTS --pval-col PVAL --ld-clump --clump-p 5e-08 --out {prefix}".format(data=data, prefix=prefix, gwascat=gwascat)
  swiss_main(args)

  file_clump = prefix + ".clump"
  file_neargwas = prefix + ".near-gwas.tab"
  file_ldgwas = prefix + ".ld-gwas.tab"

  assert_files_exist(file_clump, file_neargwas, file_ldgwas)

  gwas_ld_variants = []
  with open(str(file_ldgwas)) as infile:
    reader = csv.DictReader(infile, delimiter="\t")
    for line in reader:
      vid = line["ASSOC_EPACTS"]
      gwas_ld_variants.append(vid)

  assert "10:114754088_T/C" in gwas_ld_variants

