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

import os, sys
import pandas as pd
from pysam import VariantFile
from termcolor import *
from Variant import *
from subprocess import Popen,PIPE
from utils import *
from copy import deepcopy
from multiprocessing import *
from itertools import *

# Worker function for calculating LD buddies and looking for GWAS catalog overlap.
# This function is executed by multiprocessing pool, and needs to be importable at the top
# level of the module. Weird python thing.
def worker_ld_multi(args):
  v = args[0]
  conf = args[1]

  finder = conf['finder']
  ld_dist = conf['dist']
  ld_thresh = conf['ld_thresh']
  trait = conf['trait']
  gwascat = conf['gwascat']

  print "Working on finding GWAS catalog overlap for variant %s (%s)" % (v.vid,v.epacts)

  ld_ok = finder.compute(v.epacts,v.chrom,v.pos - ld_dist,v.pos + ld_dist,ld_thresh)

  if ld_ok:
    # In this case we also want to include the reference variant.
    finder.data[v.epacts] = (1,1) # dprime,r2

    ld_snps = {j for j in finder.data.iterkeys()}
    cat_rows = gwascat[gwascat['EPACTS'].isin(ld_snps)]

    cat_rows['ASSOC_VARIANT'] = v.vid
    cat_rows['ASSOC_EPACTS'] = v.epacts
    cat_rows['ASSOC_CHRPOS'] = v.chrpos

    vtraits = pd.Series(v.traits).dropna().unique()
    if len(vtraits) > 0:
      trait = ",".join(sorted(vtraits))
    else:
      if trait is None:
        trait = "NA"

      # Otherwise, we just use "trait" from conf as above.

    cat_rows["ASSOC_TRAIT"] = trait

    cat_rows['ASSOC_GWAS_LD'] = [finder.data.get(x)[1] for x in cat_rows['EPACTS']]

  else:
    cat_rows = None
    warning("could not calculate LD for variant %s (%s)" % (v.vid,v.chrpos))

  return cat_rows

class GWASCatalog:
  def __init__(self,filepath):
    self.filepath = filepath
    if not os.path.isfile(self.filepath):
      raise IOError, "Cannot find file: %s" % self.filepath

    self.col_vid = "VARIANT"
    self.col_epacts = "EPACTS"
    self.col_chr = "CHR"
    self.col_pos = "POS"
    self.col_trait = "PHENO"
    self.col_trait_group = "GROUP"
    #self.col_pvalue = "P_VALUE"
    self.col_logpvalue = "LOG_PVAL"
    self.col_ref = "REF"
    self.col_alt = "ALT"

    self.data = None
    self.all_cols = None

    self._load()

  def _load(self):
    """
    Load the GWAS catalog specified in self.filepath.

    Assumptions: tab-delimited file
    Post-processing:
      removes null chromosome and position rows
      adds a chrpos column e.g. 23:393919

    Returns:

    """

    self.data = pd.read_table(self.filepath,sep="\t")
    
    # Remove any entries missing CHR and POS. 
    self.data = self.data[self.data[self.col_pos].notnull()]
    self.data = self.data[self.data[self.col_chr].notnull()]

    # Remove any variants with missing identifier
    self.data = self.data[self.data[self.col_epacts].notnull()]

    # Save a list of the original columns in this catalog. 
    self.all_cols = [x for x in self.data.columns]

  def get_traits(self):
    return set(self.data[self.col_trait])

  def get_trait_groups(self):
    return set(self.data[self.col_trait_group])

  def get_trait_group_pairs(self):
    return [(x[0], sorted(x[1].PHENO.unique())) for x in self.data.groupby("GROUP")]

  def num_variants(self):
    return self.data[self.col_epacts].unique().shape[0]

  def num_rows(self):
    return self.data.shape[0]

  # Return all rows for a particular variant.
  def at_variant(self,variant):
    return self.data[self.data[self.col_epacts] == variant.epacts]

  def get_variants(self):
    def row_to_snp(row):
      s = Variant()

      s.vid = row[self.col_vid]
      s.chrom = row[self.col_chr]
      s.pos = int(row[self.col_pos])
      s.chrpos = "chr%s:%s" % (s.chrom,s.pos)
      s.ref = row[self.col_ref]
      s.alt = row[self.col_alt]
      s.epacts = "{}:{}_{}/{}".format(s.chrom,s.pos,s.ref,s.alt)

      return s

    as_variants = self.data.apply(row_to_snp,axis=1)

    return as_variants

  def variants_missing_vcf(self,vcf_file):
    cat_chroms = set(self.data[self.col_chr].unique())
    cat_variants = set(self.data[self.col_epacts].unique())

    vcf_variants = set()
    for cat_chrom in cat_chroms:
      print >> sys.stderr, "Checking chromosome %s..." % str(cat_chrom)

      if '.json' in vcf_file:
        import json
        with open(vcf_file) as jsin:
          vcf_dict = json.load(jsin)

        vcf = vcf_dict.get(cat_chrom)
        if vcf is None:
          warning("GWAS catalog has variants on chromosome %s, but could not find this chromosome in your VCF (or JSON) file: %s" % (cat_chrom,vcf_file))
          continue
      else:
        vcf = vcf_file

      vcf_pysam = VariantFile(vcf)

      # Subset catalog to chromosome
      df_cat_for_chrom = self.data.query("{} == '{}'".format(self.col_chr,cat_chrom))

      # Catalog has repeated rows for variants depending on the number of traits * citations
      # But we just need each variant once
      df_cat_for_chrom = df_cat_for_chrom.drop_duplicates(self.col_epacts)

      # Loop over subsetted catalog, check if variant is in VCF
      for idx, row in df_cat_for_chrom.iterrows():
        chrom, pos = row[self.col_chr], row[self.col_pos]

        for rec in vcf_pysam.fetch(chrom,pos,pos):
          epacts = "{}:{}_{}/{}".format(rec.chrom,rec.pos,rec.ref,rec.alt)
          vcf_variants.add(epacts)

    missing_variants = cat_variants.difference(vcf_variants)
    missing_rows = self.data[self.data[self.col_epacts].isin(missing_variants)]

    return missing_rows

  # Parallel version of variants_in_ld
  def variants_in_ld_multi(self,assoc,finder,ld_thresh=0.1,dist=1e6,num_threads=2):
    dist = int(dist)
    variants = assoc.get_variants()
    trait = assoc.trait

    # Pack up all of the settings needed for the worker process
    def get_local_settings():
      while 1:
        settings = {
          'finder' : deepcopy(finder),
          'ld_thresh' : ld_thresh,
          'dist' : dist,
          'trait' : trait,
          'gwascat' : self.data
        }

        yield settings

    # These piece of code exists to prevent creating multiprocesses
    # if thread count is only 1. Why? Because if we're already running multiprocessed,
    # like in the case of --multi-assoc, trying to spawn another multiprocess causes a crash.
    if num_threads > 1:
      # Pool of workers
      pool = Pool(num_threads)

      # Block until all of the GWAS catalog lookups per variant have completed
      results = pool.map(worker_ld_multi,izip(variants,get_local_settings()))

    else:
      results = map(worker_ld_multi,izip(variants,get_local_settings()))

    # Which variants didn't complete OK?
    failed_ld_variants = []
    for i in xrange(len(results)):
      if results[i] is None:
        failed_ld_variants.append(variants[i])

    try:
      ld_catalog = pd.concat(results)
    except:
      ld_catalog = None

    # Change column names to be more descriptive.
    if ld_catalog is not None:
      ld_catalog.rename(
        columns = dict(zip(self.all_cols,map(lambda x: "GWAS_" + x,self.all_cols))),
        inplace = True
      )

      # TODO: clean this up, string literals instead of asking the assoc object what the columns are!
      # Re-order columns.
      lead_cols = ['ASSOC_VARIANT','ASSOC_EPACTS','ASSOC_CHRPOS','ASSOC_TRAIT','GWAS_VARIANT','GWAS_EPACTS','GWAS_CHRPOS','ASSOC_GWAS_LD']
      all_cols = ld_catalog.columns.tolist()
      other_cols = filter(lambda x: x not in lead_cols,all_cols)
      col_order = lead_cols + other_cols
      ld_catalog = ld_catalog[col_order]

      # Remove unnecessary columns.
      ld_catalog = ld_catalog.drop(['GWAS_CHR','GWAS_POS'],axis=1)

    return ld_catalog, failed_ld_variants

  def variants_nearby(self,assoc,dist=1e5):
    dist = int(dist)
    variants = assoc.get_variants()

    dist_catalog = None
    for v in variants:
      # Could do better with interval tree or indexing, but performance is probably fine here (only 10K rows.)
      is_near = (self.data[self.col_chr].map(str) == str(v.chrom)) & (self.data[self.col_pos] > v.pos - dist) & (self.data[self.col_pos] < v.pos + dist)
      cat_rows = self.data[is_near]

      cat_rows['ASSOC_VARIANT'] = v.vid
      cat_rows['ASSOC_EPACTS'] = v.epacts
      cat_rows['ASSOC_CHRPOS'] = v.chrpos

      vtraits = pd.Series(v.traits).dropna().unique()
      if len(vtraits) > 0:
        trait = ",".join(sorted(vtraits))
      elif assoc.trait is not None:
        trait = assoc.trait
      else:
        trait = "NA"

      cat_rows["ASSOC_TRAIT"] = trait

      cat_rows['ASSOC_GWAS_DIST'] = abs(cat_rows[self.col_pos] - v.pos)

      dist_catalog = pd.concat([dist_catalog,cat_rows])

    # Change column names to be more descriptive.
    if dist_catalog is not None:
      dist_catalog.rename(
        columns = dict(zip(self.all_cols,map(lambda x: "GWAS_" + x,self.all_cols))),
        inplace = True
      )

      # Re-order columns.
      lead_cols = ['ASSOC_VARIANT','ASSOC_EPACTS','ASSOC_CHRPOS','ASSOC_TRAIT','GWAS_VARIANT','GWAS_EPACTS','GWAS_CHRPOS','ASSOC_GWAS_DIST']
      all_cols = dist_catalog.columns.tolist()
      other_cols = filter(lambda x: x not in lead_cols,all_cols)
      col_order = lead_cols + other_cols
      dist_catalog = dist_catalog[col_order]

      # Remove unnecessary columns.
      dist_catalog = dist_catalog.drop(['GWAS_CHR','GWAS_POS'],axis=1)

    return dist_catalog
