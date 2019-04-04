#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2014 Ryan Welch, The University of Michigan
#
# Swiss is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Swiss is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import gzip, os, subprocess, sys, hashlib, re, time, traceback
import numpy as np
import pandas as pd
import pysam
from hashlib import md5
from LDRegionCache import *
from textwrap import fill
from utils import *
from itertools import imap, chain, izip, count
from collections import defaultdict

import os, sys, tempfile, pysam, json, optparse, gzip
import pandas as pd
from subprocess import Popen, PIPE

def vcf_get_header(vcf_file):
  if vcf_file.endswith(".gz"):
    f = gzip.open(vcf_file);
  else:
    f = open(vcf_file);

  with f:
    for line in f:
      if line.startswith("#CHROM"):
        return line

def vcf_get_line(vcf_file,chrom=None,pos=None,epacts_id=None,snps_only=False):
  """
  Given a VCF file and a variant specified as either chrom/pos or an EPACTS ID, return the rows from the file
  matching the variant.

  If only chrom/pos is given, this could in theory overlap multiple variants (since SNPs and indels can have
  the same starting position.) All rows overlapping the chrom/pos will be returned.

  If an EPACTS ID is specified, this function will only return the row that matches not just the chrom/pos,
  but also the alleles specified in the EPACTS ID. This should result in a unique match.

  :param vcf_file:
  :param chrom:
  :param pos:
  :param epacts_id:
  :param snps_only:
  :return:
  """

  if epacts_id is not None:
    chrom, pos, ref, alt = parse_epacts(epacts_id)[0:4]
    alleles = (ref,alt)

  if chrom is None or pos is None:
    raise Exception, "Must provide chrom and pos, or epacts_id"

  tfile = pysam.Tabixfile(vcf_file)

  rows = []
  for row in tfile.fetch("%s:%s-%s" % (chrom,pos,pos)):
    variant = row.split("\t")
    vchrom, vpos, vid, vref, valt = variant[0:5]

    if snps_only and (len(vref) > 1 or len(valt) > 1):
      # It's not a SNP.
      continue

    if epacts_id is not None and (vref not in alleles) and (valt not in alleles):
      # If we were given an EPACTS ID to find, the variant at this position
      # must also match the alleles in the EPACTS ID. They don't, so skip it.
      continue

    rows.append(row)

  return rows

class PlinkLDSettings:
  def __init__(self,vcf_path,tabix_path,plink_path):
    # Check each path. If it doesn't exist, try to find it relative
    # to the m2zfast root directory.
    for arg,value in locals().items():
      if arg == 'self':
        continue

      path = find_systematic(value)
      if path is None or not os.path.exists(path):
        if arg == "tabix_path":
          error("cannot find tabix - please set the path in the configuration file, or make it available on your PATH.")
        else:
          error("path either does not exist or insufficient permissions to access it: %s" % str(value))
      else:
        exec "%s = \"%s\"" % (arg,path)

    if '.json' in vcf_path:
      import json
      with open(vcf_path) as jsin:
        self.vcf_path = json.load(jsin)
    else:
      self.vcf_path = vcf_path

    self.tabix_path = tabix_path
    self.plink_path = plink_path

  def create_ld_cache_key(self):
    key_string = self.vcf_path
    key = hashlib.sha512(key_string).hexdigest()
    return key

class VariantHash:
  def __init__(self):
    self.counter = count()
    self.h_hash = dict()
    self.h_variant = dict()

  def hash(self,variant):
    if variant in self.h_variant:
      return self.h_variant.get(variant)
    else:
      h = "id{}".format(next(self.counter))
      self.h_hash[h] = variant
      self.h_variant[variant] = h
      return h

  def variant(self,hashh):
    v = self.h_hash.get(hashh)
    if v is None:
      raise ValueError("Tried to unhash a variant that hadn't been hashed before: " + v)
    return v

class PlinkLDFinder():
  def __init__(self, plink_settings,cache=None,cleanup=True,verbose=False,plink_args=None):
    if not isinstance(plink_settings, PlinkLDSettings):
      raise ValueError

    self.data = {}
    self.variant = None
    self.settings = plink_settings
    self.debug = False
    self.start = None
    self.stop = None
    self.chr = None
    self.cache = cache
    self.cleanup = cleanup
    self.verbose = verbose
    self.plink_args = plink_args

    self.calc_ok = True

    self.min_r2 = None

  def write(self,filename):
    try:
      f = open(filename,'w')
      print >> f, "snp1 snp2 dprime rsquare"

      if len(self.data) == 0:
        return False

      for snp in self.data:
        print >> f, "%s %s %s %s" % (
          snp,
          self.variant,
          str(self.data.get(snp)[0]),
          str(self.data.get(snp)[1])
        )

      f.close()
    except:
      print >> sys.stderr, "Error: could not write computed LD to disk, permissions?"
      return False

    return True

  def _check_geno_paths(self):
    files = []
    if type(self.settings.vcf_path) is dict:
      map(files.append,self.settings.vcf_path.itervalues())
    else:
      files.append(self.settings.vcf_path)

    for file in files:
      if not os.path.exists(file):
        msg = fill("Error: could not find required file to generate LD "
                  "%s. Check your conf file to make sure paths are "
                  "correct. " % file)
        die(msg)

      if not os.path.exists(file + ".tbi"):
        msg = fill("Error: could not find tabix file {}, VCF must be tabix indexed and bgzipped".format(file + ".tbi"))
        die(msg)

  def compute(self, variant, chr, start, stop, min_r2=0):
    self.variant = variant
    self.start = max(0,start)
    self.stop = stop
    self.chr = chr
    self.data = None
    self.min_r2 = min_r2
    self.calc_ok = True

    # If the cache has data for this SNP and region, use it.
    # Otherwise, compute it.
    if self.cache:
      if self.cache.hasRegion(variant, start, stop):
        self.data = self.cache.getAllLD(variant)
      else:
        self._check_geno_paths()
        self.data = self._run_ld()
        self.cache.updateLD(variant, start, stop, self.data)
    else:
      self._check_geno_paths()
      self.data = self._run_ld()

    # Complete successfully?
    return self.calc_ok

  def _run_ld(self):
    if type(self.settings.vcf_path) is dict:
      vcf = self.settings.vcf_path.get(self.chr)
      if vcf is None:
        self.calc_ok = False
        print >> sys.stderr, "Error: no VCF file available for chromosome '%s' - maybe a X vs 23 issue?" % self.chr
        return {}
    else:
      vcf = self.settings.vcf_path

    region = "%s:%s-%s" % (self.chr,self.start,self.stop)

    data = {}
    try:
      data = self._ld_refvar_region(vcf,self.variant,region,min_r2=self.min_r2,ignore_filter=True)
    except Exception as e:
      if SWISS_DEBUG: raise

      print >> sys.stderr, e.message
      self.calc_ok = False

    return data

  def _ld_refvar_region(self,vcf_file,variant,region,min_r2=0,ignore_filter=False,ignore_indels=False):
    """
    Given a VCF file, variant (EPACTS ID), and a region, compute LD between the variant and all other variants
    in the region. This function uses PLINK 1.9 to do the LD calculations.

    Returns a data frame with the following columns:
    CHROM, POS - obvious
    SNP - EPACTS ID
    MAF - obvious
    R2 - obvious
    DP - D'

    Args:
      vcf_file: path to VCF file
      variant: variant in EPACTS format (e.g. chr:pos_ref/alt)
      region: region within to compute, specify like "chr:start-end"
      min_r2: only return LD pairs with r2 > min value
      ignore_filter: should we ignore the filter column in the VCF?
      ignore_indels: should we skip over indel variants in the VCF?
    """

    # Use tabix to pull out the region of interest from the VCF file.
    # plink1.9 does not have the ability to do this on its own.
    tabix_cmd = "{tabix} -h {0} {1}".format(vcf_file,region,tabix=self.settings.tabix_path)
    proc_tabix = Popen(tabix_cmd,shell=True,stdout=PIPE,stderr=PIPE)

    # We need to hash variant names for now (PLINK limits length of string)
    hasher = VariantHash()
    variant_hash = hasher.hash(variant)

    # Use plink1.9 to calculate LD - we'll feed it VCF lines directly to its STDIN.
    tmpout = tempfile.mktemp(dir=os.getcwd())
    extra_args = self.plink_args if self.plink_args is not None else ""
    plink_cmd = "{plink} --vcf /dev/fd/0 --r2 gz dprime with-freqs yes-really --ld-snp {0} --ld-window-kb 99999 --ld-window 99999 --threads 1 " \
      "--ld-window-r2 {min_r2} {extra_args} --out {1}".format(variant_hash,tmpout,plink=self.settings.plink_path,min_r2=min_r2,extra_args=extra_args)
    proc_ld = Popen(
      plink_cmd,
      shell=True,
      stdin=PIPE,
      stdout=PIPE,
      stderr=PIPE
    )

    def cleanup():
      for ext in [".log",".nosex",".ld.gz","-temporary.bed","-temporary.bim","-temporary.fam"]:
        delfile = tmpout + ext
        try:
          os.remove(delfile)
        except:
          pass

    # Loop over VCF lines, changing the ID to be a more descriptive EPACTS ID (chr:pos_ref/alt_id).
    for line in proc_tabix.stdout:
      if line.startswith("#"):
        proc_ld.stdin.write(line)
        continue

      ls = line.split("\t")
      ls[-1] = ls[-1].rstrip()

      chrom, pos, vid, ref, alt, qual, filt = ls[0:7]

      # If this variant isn't marked as PASS, we shouldn't consider it.
      if not ignore_filter and "PASS" not in filt:
        continue

      # Skip indel?
      if (len(ref) > 1 or len(alt) > 1) and ignore_indels:
        continue

      # If this variant already has an ID, we'll tack it on to the end of the EPACTS ID.
      # Skipping this - the "extra" part of the EPACTS ID was causing mismatches in other places...
      # if vid != ".":
      #   new_id = "%s:%s_%s/%s_%s" % (chrom,pos,ref,alt,vid)
      # else:
      #   new_id = "%s:%s_%s/%s" % (chrom,pos,ref,alt)

      # Make an EPACTS ID for this variant
      new_id = "%s:%s_%s/%s" % (chrom,pos,ref,alt)

      # We have to hash the variant name going to PLINK (it can't handle very long strings > a few thousand characters)
      id_hash = hasher.hash(new_id)

      # Write the variant out to plink1.9's STDIN.
      print >> proc_ld.stdin, "\t".join([chrom,pos,id_hash,ref,alt] + ls[5:])

    # Wait for plink to finish calculating LD.
    ld_stdout, ld_stderr = proc_ld.communicate()

    if ld_stderr != '' or 'Error' in ld_stdout or proc_ld.returncode != 0:
      if not SWISS_DEBUG:
        cleanup()

      raise Exception, "\n" + ld_stdout + "\n\n" + ld_stderr

    # Return a data frame of LD statistics.
    df = pd.read_table(tmpout + ".ld.gz",sep="\s+",compression="gzip")

    # Drop the row that corresponds to LD with itself
    df = df[df.SNP_B != df.SNP_A]

    # Convert from variant hashes back to real variant IDs
    df.loc[:,"SNP_A"] = df.loc[:,"SNP_A"].map(hasher.variant)
    df.loc[:,"SNP_B"] = df.loc[:,"SNP_B"].map(hasher.variant)

    # Small changes to names.
    df.rename(columns = lambda x: x.replace("CHR","CHROM"),inplace=True)
    df.rename(columns = lambda x: x.replace("BP","POS"),inplace=True)
    df.rename(columns = lambda x: x.replace("SNP","VARIANT"),inplace=True)

    # Run cleanup
    if not SWISS_DEBUG:
      cleanup()

    # Slight modification for swiss: it expects return to be dictionary of variant --> (dprime,rsq)
    ld_data = dict(zip(df.VARIANT_B,df.apply(lambda x: tuple([x["DP"],x["R2"]]),axis=1)))

    return ld_data

# def ld_snp_pair(vcf_file,variant1,variant2,min_r2=0,ignore_filter=False,ignore_indels=False):
#   """
#   Given a VCF file and two variants, compute LD between the variants.
#   This function uses PLINK 1.9 to do the LD calculations.
#
#   Requires:
#
#   TABIX_PATH - global constant set to the tabix path, or just "tabix" if it is on the user's path
#   PLINK_PATH - same as above, except for PLINK (must be 1.9, NOT Sean Purcell's older version)
#
#   Returns a data frame with the following columns:
#   CHR, BP - obvious
#   SNP - EPACTS ID
#   MAF - obvious
#   R2 - obvious
#   DP - D'
#   """
#
#   parsed_v1 = parse_epacts(variant1)
#   parsed_v2 = parse_epacts(variant2)
#
#   # The two variants should be on the same chromosome.
#   if parsed_v1[0] != parsed_v2[0]:
#     error("reference variant is not on the same chromosome as other variant")
#
#   # Use plink1.9 to calculate LD - we'll feed it VCF lines directly to its STDIN.
#   tmpout = tempfile.mktemp(dir=os.getcwd())
#   proc_ld = Popen(
#     "{plink} --vcf /dev/fd/0 --r2 gz dprime with-freqs yes-really --ld-snps {0} --ld-window-kb 2000000000 --ld-window 2000000000 --threads 1 "
#     "--ld-window-r2 {min_r2} --out {1}".format(variant1,tmpout,plink=PLINK_PATH,min_r2=min_r2),
#     shell=True,
#     stdin=PIPE,
#     stdout=PIPE,
#     stderr=PIPE
#   )
#
#   # Send the header to PLINK
#   print >> proc_ld.stdin, vcf_get_header(vcf_file)
#
#   # Grab just the two lines we need (for both variants) from the VCF and send them to plink directly
#   # Note: they need to be sorted, for some odd reason.
#   lines = []
#
#   # If it's the same variant twice, only give PLINK the variant once. It doesn't like it twice for some reason.
#   if parsed_v1[0:4] == parsed_v2[0:4]:
#     lines.extend(vcf_get_line(vcf_file,epacts_id=variant1))
#   else:
#     # The two variants aren't the same.
#     # But PLINK does expect them to come in sorted order.
#     if parsed_v1[1] <= parsed_v2[1]:
#       lines.extend(vcf_get_line(vcf_file,epacts_id=variant1))
#       lines.extend(vcf_get_line(vcf_file,epacts_id=variant2))
#     else:
#       lines.extend(vcf_get_line(vcf_file,epacts_id=variant2))
#       lines.extend(vcf_get_line(vcf_file,epacts_id=variant1))
#
#   # Loop over VCF lines, changing the ID to be a more descriptive EPACTS ID (chr:pos_ref/alt_id).
#   for line in lines:
#     ls = line.split("\t")
#     ls[-1] = ls[-1].rstrip()
#
#     chrom, pos, vid, ref, alt, qual, filt = ls[0:7]
#
#     # If this variant isn't marked as PASS, we shouldn't consider it.
#     if not ignore_filter and filt != "PASS":
#       continue
#
#     if (len(ref) > 1 or len(alt) > 1) and ignore_indels:
#       continue
#
#     # If this variant already has an ID, we'll tack it on to the end of the EPACTS ID.
#     # Skipping this - the "extra" part of the EPACTS ID was causing mismatches in other places...
#     # if vid != ".":
#     #   new_id = "%s:%s_%s/%s_%s" % (chrom,pos,ref,alt,vid)
#     # else:
#     #   new_id = "%s:%s_%s/%s" % (chrom,pos,ref,alt)
#
#     # Make an EPACTS ID for this variant, and set it to the ID column in the row.
#     new_id = "%s:%s_%s/%s" % (chrom,pos,ref,alt)
#
#     if (len(ref) + len(alt)) > 10000:
#       trunc_id = "%s:%s_%s.../%s..." % (chrom,pos,ref[0:10],alt[0:10])
#       warn("variant alleles are too long for PLINK 1.9, skipping (truncated alleles to first 10 bp): %s" % trunc_id)
#       continue
#
#     # Write the variant out to plink1.9's STDIN.
#     print >> proc_ld.stdin, "\t".join([chrom,pos,new_id,ref,alt] + ls[5:])
#
#   # Wait for plink to finish calculating LD.
#   ld_stdout, ld_stderr = proc_ld.communicate()
#
#   if ld_stderr != '' or 'Error' in ld_stdout or proc_ld.returncode != 0:
#     raise Exception, "\n" + ld_stdout + "\n\n" + ld_stderr
#
#   # Return a data frame of LD statistics.
#   df = pd.read_table(tmpout + ".ld.gz",sep="\s+",compression="gzip")
#
#     # If var1 and var2 aren't the same, drop rows where the two are equivalent.
#   if parsed_v1[0:4] != parsed_v2[0:4]:
#     df = df[df.SNP_B != df.SNP_A]
#
#   # Small changes to names.
#   df.rename(columns = lambda x: x.replace("CHR","CHROM"),inplace=True)
#   df.rename(columns = lambda x: x.replace("BP","POS"),inplace=True)
#   df.rename(columns = lambda x: x.replace("SNP","VARIANT"),inplace=True)
#
#   # This just returns the first row as a data frame. If you do 0,: then it just returns a series of the first row.
#   df = df.iloc[0:1,:]
#
#   # Cleanup temporary files.
#   for ext in [".log",".nosex",".ld.gz","-temporary.bed","-temporary.bim","-temporary.fam"]:
#     delfile = tmpout + ext
#     try:
#       os.remove(delfile)
#     except:
#       pass
#
#   return df
