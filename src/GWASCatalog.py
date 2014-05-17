#!/usr/bin/env python
import pandas as pd
import os
from termcolor import *
from Variant import *
from subprocess import Popen,PIPE
from utils import *
from copy import deepcopy
from multiprocessing import *
from itertools import *
import sys

# Worker function for calculating LD buddies and looking for GWAS catalog overlap.
# This function is executed by multiprocessing pool, and needs to be importable at the top
# level of the module. Weird python thing.
def worker_ld_multi(args):
  v = args[0];
  conf = args[1];

  finder = conf['finder'];
  dist = conf['dist'];
  ld_thresh = conf['ld_thresh'];
  trait = conf['trait'];
  gwascat = conf['gwascat'];

  print "Working on finding GWAS catalog overlap for variant %s (%s)" % (v.name,v.chrpos);

  ld_ok = finder.compute(v.chrpos,v.chrom,v.pos - dist, v.pos + dist,ld_thresh);
  if ld_ok:
    ld_snps = {j for j in finder.data.iterkeys()};
    cat_rows = gwascat[gwascat['CHRPOS'].isin(ld_snps)];

    cat_rows['ASSOC_MARKER'] = v.name;
    cat_rows['ASSOC_CHRPOS'] = v.chrpos;

    if trait is not None:
      cat_rows['ASSOC_TRAIT'] = trait;
    else:
      cat_rows['ASSOC_TRAIT'] = "NA";

    cat_rows['ASSOC_GWAS_LD'] = [finder.data.get(x)[1] for x in cat_rows['CHRPOS']];

  else:
    cat_rows = None;
    warning("could not calculate LD for variant %s (%s) - you should try a different source of LD information to properly clump these variants." % (v.name,v.chrpos));

  return cat_rows;

class GWASCatalog:
  def __init__(self,filepath):
    self.filepath = filepath;
    if not os.path.isfile(self.filepath):
      raise IOError, "Cannot find file: %s" % self.filepath;

    self.col_snp = "SNP";
    self.col_chr = "CHR";
    self.col_pos = "POS";
    self.col_trait = "PHENO";
    self.col_trait_group = "Group";
    self.col_pvalue = "P_VALUE";
    self.build = "hg18";

    self.data = None;
    self.all_cols = None;
    self.pos_index = {};

    self._load();

  # Load the GWAS catalog specified in self.filepath. 
  # Assumptions: tab-delimited file
  # Post-processing: 
  #   removes null chromosome and position rows
  #   adds a chrpos column e.g. 23:393919
  def _load(self):
    self.data = pd.read_table(self.filepath,sep="\t");
    
    # Remove any entries missing CHR and POS. 
    self.data = self.data[self.data[self.col_pos].notnull()];
    self.data = self.data[self.data[self.col_chr].notnull()];

    # Create a chrpos variant column e.g. X:29291
    self.data['CHRPOS'] = self.data[self.col_chr].map(str) + ":" + self.data[self.col_pos].map(lambda x: str(int(x)));

    # Save a list of the original columns in this catalog. 
    self.all_cols = [x for x in self.data.columns];

    # Create position index. 
    self._reindex();

  def _reindex(self):
    self.pos_index = {};

    # Create position index. 
    for row_index, row in self.data.iterrows():
      row_chrom = row[self.col_chr];
      row_pos = row[self.col_pos];
      self.pos_index.setdefault((row_chrom,row_pos),[]).append(row_index);

  def get_traits(self):
    return set(self.data[self.col_trait]);

  def get_trait_groups(self):
    return set(self.data[self.col_trait_group]);

  def get_trait_group_pairs(self):
    return [(x[0], sorted(x[1].PHENO.unique())) for x in self.data.groupby("Group")];

  def num_variants(self):
    return self.data.shape[0];

  # Return all rows at a given position. 
  def at_position(self,chrom,pos):
    indices = self.pos_index.get((chrom,pos));
    if indices is not None:
      return self.data.ix[indices];
    else:
      return None;

  # Return all rows for a particular variant's chrom/pos. 
  def at_variant(self,variant):
    return self.at_position(variant.chrom,variant.pos);

  def get_variants(self):
    def row_to_snp(row):
      s = Variant();
      s.name = row[self.col_snp];
      s.chrom = row[self.col_chr];
      s.pos = int(row[self.col_pos]);
      s.chrpos = "chr%s:%s" % (s.chrom,s.pos);
      return s;

    as_snps = self.data.apply(row_to_snp,axis=1);

    return as_snps;

  def variants_missing_vcf(self,vcf_file):
    chrpos_catalog = set(self.data['CHRPOS']);
    
    cat_chroms = set([i.split(":")[0] for i in chrpos_catalog]);
    
    vcf_regions = set();
    for cat_chrom in cat_chroms:
      tabix_regions = set();
      for v in chrpos_catalog:
        (chrom,pos) = v.split(":");
        tabix_str = "%s:%s-%s" % (chrom,pos,pos);
        tabix_regions.add(tabix_str);

      if '.json' in vcf_file:
        import json
        with open(vcf_file) as jsin:
          vcf_dict = json.load(jsin);
        
        vcf = vcf_dict.get(cat_chrom);
        if vcf is None:
          warning("Could not find chromosome %s according to VCF specification (JSON) file: %s" % (cat_chrom,vcf_file));
          continue;
      else:
        vcf = vcf_file;

      cmd = "tabix {vcf} {regions} | cut -f 1-3".format(
        vcf = vcf,
        regions = " ".join(tabix_regions)
      );

      proc = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE);
      (stdout,stderr) = proc.communicate();

      if stderr != '':
        print colored("Error: ",'red') + "could not execute tabix on VCF file %s, error was:\n %s" % (vcf,stderr);
        raise Exception;

      for line in stdout.split("\n"):
        if line.strip() == '' or line is None:
          continue;

        lsplit = line.split();
        vcf_regions.add("{chrom}:{pos}".format(chrom=lsplit[0],pos=lsplit[1]));

    missing_chrpos = chrpos_catalog.difference(vcf_regions);
    missing_rows = self.data[self.data['CHRPOS'].isin(missing_chrpos)];

    return missing_rows;

  # # Given a variant, find other variants in this catalog in LD with it.
  # # Uses whichever finder/LD source you choose to provide to calculate LD.
  # def variants_in_ld(self,assoc,finder,ld_thresh=0.1,dist=1e6):
  #   dist = int(dist);
  #   variants = assoc.get_snps();
  #   trait = assoc.trait;
  #
  #   ld_catalog = None;
  #   failed_ld_variants = [];
  #   for v in variants:
  #     print "Working on variant %s (%s)" % (v.name,v.chrpos);
  #
  #     ld_ok = finder.compute(v.chrpos,v.chrom,v.pos - dist, v.pos + dist,ld_thresh);
  #     if ld_ok:
  #       ld_snps = {v for v in finder.data.iterkeys()};
  #       #cat_rows = self.data[self.data['CHRPOS'].map(lambda x: x in ld_snps)];
  #       cat_rows = self.data[self.data['CHRPOS'].isin(ld_snps)];
  #
  #       cat_rows['ASSOC_MARKER'] = v.name;
  #       cat_rows['ASSOC_CHRPOS'] = v.chrpos;
  #
  #       if trait is not None:
  #         cat_rows['ASSOC_TRAIT'] = trait;
  #       else:
  #         cat_rows['ASSOC_TRAIT'] = "NA";
  #
  #       cat_rows['ASSOC_GWAS_LD'] = [finder.data.get(x)[1] for x in cat_rows['CHRPOS']];
  #
  #       ld_catalog = pd.concat([ld_catalog,cat_rows]);
  #     else:
  #       failed_ld_variants.append(v);
  #
  #       warning("could not calculate LD for variant %s (%s) - you should try a different source of LD information to properly clump these variants." % (v.name,v.chrpos));
  #
  #   # Change column names to be more descriptive.
  #   if ld_catalog is not None:
  #     ld_catalog.rename(
  #       columns = dict(zip(self.all_cols,map(lambda x: "GWAS_" + x,self.all_cols))),
  #       inplace = True
  #     );
  #
  # #    ld_catalog.rename(columns = {
  # #      'SNP' : 'GWAS_SNP',
  # #      'CHR' : 'GWAS_CHR',
  # #      'POS' : 'GWAS_POS',
  # #      'CHRPOS' : "GWAS_CHRPOS"
  # #    },inplace=True);
  #
  #     # TODO: clean this up, string literals instead of asking the assoc object what the columns are!
  #     # Re-order columns.
  #     lead_cols = ['ASSOC_MARKER','ASSOC_CHRPOS','ASSOC_TRAIT','GWAS_SNP','GWAS_CHRPOS','ASSOC_GWAS_LD'];
  #     all_cols = ld_catalog.columns.tolist();
  #     other_cols = filter(lambda x: x not in lead_cols,all_cols);
  #     col_order = lead_cols + other_cols;
  #     ld_catalog = ld_catalog[col_order];
  #
  #     # Remove unnecessary columns.
  #     ld_catalog = ld_catalog.drop(['GWAS_CHR','GWAS_POS'],axis=1);
  #
  #   return ld_catalog, failed_ld_variants;

  # Parallel version of variants_in_ld
  def variants_in_ld_multi(self,assoc,finder,ld_thresh=0.1,dist=1e6,num_threads=2):
    dist = int(dist);
    variants = assoc.get_snps();
    trait = assoc.trait;

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

        yield settings;

    # These piece of code exists to prevent creating multiprocesses
    # if thread count is only 1. Why? Because if we're already running multiprocessed,
    # like in the case of --multi-assoc, trying to spawn another multiprocess causes a crash.
    if num_threads > 1:
      # Pool of workers
      pool = Pool(num_threads);

      # Block until all of the GWAS catalog lookups per variant have completed
      results = pool.map(worker_ld_multi,izip(variants,get_local_settings()));

    else:
      results = map(worker_ld_multi,izip(variants,get_local_settings()));

    # Which variants didn't complete OK?
    failed_ld_variants = [];
    for i in xrange(len(results)):
      if results[i] is None:
        failed_ld_variants.append(variants[i]);

    try:
      ld_catalog = pd.concat(results);
    except:
      ld_catalog = None;

    # Change column names to be more descriptive.
    if ld_catalog is not None:
      ld_catalog.rename(
        columns = dict(zip(self.all_cols,map(lambda x: "GWAS_" + x,self.all_cols))),
        inplace = True
      );

      # TODO: clean this up, string literals instead of asking the assoc object what the columns are!
      # Re-order columns.
      lead_cols = ['ASSOC_MARKER','ASSOC_CHRPOS','ASSOC_TRAIT','GWAS_SNP','GWAS_CHRPOS','ASSOC_GWAS_LD'];
      all_cols = ld_catalog.columns.tolist();
      other_cols = filter(lambda x: x not in lead_cols,all_cols);
      col_order = lead_cols + other_cols;
      ld_catalog = ld_catalog[col_order];

      # Remove unnecessary columns.
      ld_catalog = ld_catalog.drop(['GWAS_CHR','GWAS_POS'],axis=1);

    return ld_catalog, failed_ld_variants;

  def variants_nearby(self,assoc,dist=1e5):
    dist = int(dist);
    variants = assoc.get_snps();
    trait = assoc.trait;

    dist_catalog = None;
    for v in variants:
      # Could do better with interval tree or indexing, but performance is probably fine here (only 10K rows.)
      is_near = (self.data['CHR'].map(str) == str(v.chrom)) & (self.data['POS'] > v.pos - dist) & (self.data['POS'] < v.pos + dist);
      cat_rows = self.data[is_near];

      cat_rows['ASSOC_MARKER'] = v.name;
      cat_rows['ASSOC_CHRPOS'] = v.chrpos;

      if trait is not None:
        cat_rows['ASSOC_TRAIT'] = trait;
      else:
        cat_rows['ASSOC_TRAIT'] = "NA";

      cat_rows['ASSOC_GWAS_DIST'] = abs(cat_rows['POS'] - v.pos);

      dist_catalog = pd.concat([dist_catalog,cat_rows]);

    # Change column names to be more descriptive.
    if dist_catalog is not None:
      dist_catalog.rename(
        columns = dict(zip(self.all_cols,map(lambda x: "GWAS_" + x,self.all_cols))),
        inplace = True
      );

      # Re-order columns.
      lead_cols = ['ASSOC_MARKER','ASSOC_CHRPOS','ASSOC_TRAIT','GWAS_SNP','GWAS_CHRPOS','ASSOC_GWAS_DIST'];
      all_cols = dist_catalog.columns.tolist();
      other_cols = filter(lambda x: x not in lead_cols,all_cols);
      col_order = lead_cols + other_cols;
      dist_catalog = dist_catalog[col_order];

      # Remove unnecessary columns.
      dist_catalog = dist_catalog.drop(['GWAS_CHR','GWAS_POS'],axis=1);

    return dist_catalog; 
