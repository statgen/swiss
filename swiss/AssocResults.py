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

import pandas as pd
import numpy as np
from .utils import *
from bx.intervals.intersection import *
from Variant import *

LOAD_CHUNKSIZE = 500000

class ChromTree:
  def __init__(self):
    self.chrom = {}

  def add_position(self,chrom,pos):
    self.chrom.setdefault(chrom,IntervalTree()).add(pos,pos)

  def add_interval(self,chrom,start,end):
    self.chrom.setdefault(chrom,IntervalTree()).add(start,end)

  def find_overlap(self,chrom,start,end=None):
    node = self.chrom.get(chrom)
    if node is None:
      return None
    else:
      if end is None:
        end = start
      return node.find(start,end)

def sort_genome(dframe,chr_col,pos_col):
  dframe['int_chr'] = dframe[chr_col].map(chrom2chr)
  dframe = dframe.sort_values(['int_chr',pos_col])
  del dframe['int_chr']

  return dframe

def fix_chrom(x):
  try:
    x = str(int(float(x)))
  except:
    x = str(x)

  return x.replace("chr","")

def get_header(filepath,sep="\t"):
  with open(filepath) as f:
    header = f.readline().rstrip().split(sep)

  return header

def filter_imp_quality(dframe,rsq_col,threshold=0.3):
  try:
    threshold = float(threshold)
  except:
    error("Imputation quality threshold is not floatable, got: %s" % str(threshold))

  # RSQ needs to be float in order to threshold it correctly
  dframe[rsq_col] = dframe[rsq_col].astype("float")

  # Filter
  dframe = dframe[(dframe[rsq_col] >= threshold) | dframe[rsq_col].isnull()]

  return dframe

def nskip_double_pound(fpath):
  if fpath.endswith(".gz"):
    fp = gzip.open(fpath,"rt")
  else:
    fp = open(fpath)

  count = 0
  with fp:
    for line in fp:
      if line.startswith("##"):
        count += 1
        continue
      else:
        break

  return count

class AssocResults:
  def __init__(self,filepath=None,trait=None,df=None,pval_thresh=None,rsq_filter=None,query=None):
    """
    Construct an AssocResults object.
    :param filepath: Path to a file containing association results
    :param trait: Associated trait name
    :param df: Data frame of association results. Filepath will be ignored and the results from this data frame will be used.
    :param pval_thresh: If specified, will filter results at this p-value threshold
    :param rsq_filter: Filter on imputation accuracy
    :param filter: Arbitrary filter based on numexpr syntax. For example, "BETA > 0.5 & AL_FREQ < 0.001"
    :return: This object
    """

    self.filepath = filepath
    self.data = df
    self.trait = trait

    self.pval_col = "P-value"
    self.logp_col = "LOGPVALUE"
    self.vid_col = "MarkerName"
    self.chrom_col = "chr"
    self.pos_col = "pos"
    self.rsq_col = "RSQ"
    self.ref_col = "REF"
    self.alt_col = "ALT"
    self.trait_col = "TRAIT"

    # Internal column generated with EPACTS formatted IDs
    self.epacts_col = "SWISS_VARIANT"

    self.pval_thresh = pval_thresh
    self.rsq_filter = rsq_filter
    self.query = query

  def load(self,*args,**kwargs):
    data = self.data

    # If a data frame isn't already loaded, we need to read it from a file.
    if data is None:
      # Specify data types for loading.
      # Positions should be fine with a 32-bit int.
      dtypes = {
        self.pos_col: pd.np.uint32,
        self.pval_col: str
      }

      # Imputation accuracies only require a 32-bit float (maybe even a half float would be fine)
      header = get_header(self.filepath)
      if self.rsq_col in header:
        dtypes[self.rsq_col] = pd.np.float32

      # Pandas can't skip rows automatically if comments are denoted by > 1 character
      # So we need to skip them manually...
      skip_count = nskip_double_pound(self.filepath)

      df_iter = pd.read_table(self.filepath,na_values=["NA","None","."],iterator=True,chunksize=LOAD_CHUNKSIZE,skiprows=skip_count,dtype=dtypes,*args,**kwargs)

      chunks = []
      for chunk in df_iter:
        if self.logp_col not in chunk.columns:
          if self.pval_col not in chunk.columns:
            error("Cannot find p-value column (or log p-value column) in your data: please specify --pval-col or --logp-col")

          chunk[self.logp_col] = chunk[self.pval_col].map(convert_to_log10)

        if self.pval_thresh is not None:
          chunk = chunk[chunk[self.logp_col] < float(self.pval_thresh.log10())]
        else:
          chunk[self.pval_col] = np.random.uniform(size=chunk.shape[0])
          chunk[self.logp_col] = np.log10(chunk[self.pval_col])

        if self.rsq_filter is not None:
          chunk = filter_imp_quality(chunk,self.rsq_col,self.rsq_filter)

        if self.query is not None:
          chunk = chunk.query(self.query)

        chunks.append(chunk)

      data = pd.concat(chunks)

    else:
      # We still need to run our filters, even if a data frame was passed in.
      if self.pval_thresh is not None:
        # Here, they've already given us a p-value column. If they have high precision p-values, we can only hope they
        # passed a column of Decimal() objects. What's great is that -np.log10(Decimal(x)) returns Decimal, so the code
        # below works in either float or Decimal cases.
        if self.logp_col not in data:
          data[self.logp_col] = np.log10(data[self.pval_col])
        data = data[data[self.logp_col] < float(self.pval_thresh.log10())]
      else:
        # If no p-value threshold was given, it means give each variant a random p-value
        # to use when pruning (basically, a random SNP within each LD clump will end up getting
        # selected
        data[self.pval_col] = np.random.uniform(size=data.shape[0])
        data[self.logp_col] = np.log10(data[self.pval_col])

      if self.rsq_filter is not None:
        data = filter_imp_quality(data,self.rsq_col,self.rsq_filter)

      if self.query is not None:
        data = data.query(self.query)

    for col in ('vid_col',):
      if self.__dict__[col] not in data.columns:
        error("column '%s' does not exist in your association results data" % self.__dict__[col])

    has_epacts_cols = (self.chrom_col in data) &\
                      (self.pos_col in data) &\
                      (self.ref_col in data) &\
                      (self.alt_col in data)

    if has_epacts_cols:
      # If the data file has the needed columns, we can create an EPACTS formatted ID for each variant.
      data.loc[:, self.epacts_col] = data[self.chrom_col].astype("str").str.replace("chr", "") +\
                                     ":" + \
                                     data[self.pos_col].astype("str") +\
                                     "_" + \
                                     data[self.ref_col] +\
                                     "/" + \
                                     data[self.alt_col]
    else:
      # Well, we don't have enough to assemble an EPACTS ID from columns. So we need to construct it from the
      # marker column, if possible.
      if self.chrom_col not in data:
        data[self.chrom_col] = ""

      if self.pos_col not in data:
        data[self.pos_col] = -1

      data[self.epacts_col] = ""

      for index, row in data.iterrows():
        variant = data.at[index,self.vid_col]
        chrom, pos, ref, alt, _ = parse_epacts(variant)
        chrom = chrom.replace("chr","")

        # If they didn't give us chrom/pos as columns, we need to fill it in.
        data.at[index,self.chrom_col] = chrom
        data.at[index,self.pos_col] = pos

        # Insert normalized EPACTS ID (chop off the junk at the very end, if there was any)
        data.at[index,self.epacts_col] = "{}:{}_{}/{}".format(chrom,pos,ref,alt)

    # Drop variants that do not have chromosome
    data = data[data[self.chrom_col].notnull()]

    # Try to fix chromosome column.
    # Should just be 1, 2, 3, X, Y and not chr4 or chrX
    data[self.chrom_col] = data[self.chrom_col].map(fix_chrom)

    # Drop SNPs that do not have a position.
    data = data[data[self.pos_col].notnull()]

    # Position must be an integer or it doesn't make sense.
    data[self.pos_col] = data[self.pos_col].astype('int')

    # Try to insert the trait as a column.
    if self.trait_col in data.columns:
      print "Using trait column %s in association reults file..." % self.trait_col
    else:
      data[self.trait_col] = self.trait

    # Make sure reference is updated
    self.data = data

  def liftover(self,build):
    pass

  # Drop a list of variants from the data.
  # This will remove any variant with the same name OR the same chrom/pos.
  # def drop_variants(self,variant_list):
  #   for v in variant_list:
  #     # Drop any variant with the same name (marker ID.)
  #     self.data = self.data[self.data[self.marker_col] != v.vid]
  #
  #     # Drop any variant with the same chrom/pos.
  #     data_chrpos = self.data[self.chrom_col] + ":" + self.data[self.pos_col].map(lambda x: str(int(x)))
  #     variant_chrpos = v.as_chrpos()
  #     self.data = self.data[data_chrpos != variant_chrpos]

  # Keep a list of variants and discard the rest.
  # Variants are kept if their names directly match, or they are at the same chrom/pos.
  # def keep_variants(self,variant_list):
  #   names = [v.vid for v in variant_list]
  #   chrpos = [v.as_chrpos() for v in variant_list]
  #
  #   data_chrpos = self.data[self.chrom_col] + ":" + self.data[self.pos_col].map(lambda x: str(int(x)))
  #
  #   name_match = self.data[self.marker_col].isin(names)
  #   chrpos_match = data_chrpos.isin(chrpos)
  #
  #   self.data = self.data[name_match | chrpos_match]

  def dist_clump(self,p_thresh=5e-08,dist=50E3):
    # TODO: fix this mess (why return the data re-ordered, and why even subset at all - just use .at and labels instead of .iat)
    # ...

    self.data = self.data[self.data[self.logp_col] < float(p_thresh.log10())]
    self.data = self.data.sort_values(self.logp_col)

    self.data['drop_row'] = False
    begin_cols = [self.epacts_col, self.pval_col, self.logp_col, self.chrom_col, self.pos_col, 'drop_row']
    col_order = begin_cols + list(set(self.data.columns.tolist()).difference(begin_cols))
    self.data = self.data[col_order]

    if self.data.shape[0] <= 0:
      return

    tree = ChromTree()
    for i in xrange(self.data.shape[0]):
      chrom = self.data.iat[i,3]
      pos = self.data.iat[i,4]

      if tree.find_overlap(chrom,pos - dist,pos + dist):
        self.data.iat[i,5] = True # drop column
      else:
        tree.add_position(chrom,pos)

    final = self.data[self.data.drop_row == False]
    final = sort_genome(final,self.chrom_col,self.pos_col)
    del final['drop_row']

    self.data = final
    return final

  def get_variants(self):
    snps = []
    for row_index, row in self.data.iterrows():
      s = Variant()

      s.vid = row[self.vid_col]
      s.epacts = row[self.epacts_col]
      s.chrom = row[self.chrom_col]
      s.pos = row[self.pos_col]
      s.chrpos = "%s:%s" % (s.chrom,s.pos)

      s.traits.append(row[self.trait_col])

      snps.append(s)

    return snps

  def has_rows(self):
    return self.data.shape[0] <= 0
