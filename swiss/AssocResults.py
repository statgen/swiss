#!/usr/bin/env python
import pandas as pd
from utils import *
from bx.intervals.intersection import *
from Variant import *

LOAD_CHUNKSIZE = 500000

class ChromTree:
  def __init__(self):
    self.chrom = {};

  def add_position(self,chrom,pos):
    self.chrom.setdefault(chrom,IntervalTree()).add(pos,pos);

  def add_interval(self,chrom,start,end):
    self.chrom.setdefault(chrom,IntervalTree()).add(start,end);

  def find_overlap(self,chrom,start,end=None):
    node = self.chrom.get(chrom);
    if node is None:
      return None;
    else:
      if end is None:
        end = start;
      return node.find(start,end);

def sort_genome(dframe,chr_col,pos_col):
  dframe['int_chr'] = dframe[chr_col].map(chrom2chr);
  dframe = dframe.sort(['int_chr',pos_col]);
  del dframe['int_chr'];

  return dframe;

def fix_chrom(x):
  try:
    x = str(int(float(x)));
  except:
    x = str(x);

  return x.replace("chr","");

def get_header(filepath,sep="\t"):
  with open(filepath) as f:
    header = f.readline().rstrip().split(sep)

  return header

def filter_imp_quality(dframe,rsq_col,threshold=0.3):
  try:
    threshold = float(threshold);
  except:
    error("Imputation quality threshold is not floatable, got: %s" % str(threshold));

  # RSQ needs to be float in order to threshold it correctly 
  dframe[rsq_col] = dframe[rsq_col].astype("float");

  # Filter
  dframe = dframe[(dframe[rsq_col] >= threshold) | dframe[rsq_col].isnull()];

  return dframe

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

    self.filepath = filepath;
    self.data = df;
    self.trait = trait;

    self.pval_col = "P-value";
    self.marker_col = "MarkerName";
    self.chrom_col = "chr";
    self.pos_col = "pos";
    self.rsq_col = "RSQ";

    self.pval_thresh = pval_thresh
    self.rsq_filter = rsq_filter
    self.query = query

  def load(self,*args,**kwargs):
    # If a data frame isn't already loaded, we need to read it from a file.
    if self.data is None:
      # Specify data types for loading.
      # Positions should be fine with a 32-bit int.
      dtypes = {
        self.pos_col : pd.np.uint32
      }

      # Imputation accuracies only require a 32-bit float (maybe even a half float would be fine)
      header = get_header(self.filepath)
      if self.rsq_col in header:
        dtypes[self.rsq_col] = pd.np.float32

      compr = 'gzip' if is_gzip(self.filepath) else None;
      df_iter = pd.read_table(self.filepath,compression=compr,na_values=["NA","None","."],iterator=True,chunksize=LOAD_CHUNKSIZE,dtype=dtypes,*args,**kwargs);

      chunks = []
      for chunk in df_iter:
        if self.pval_thresh is not None:
          chunk = chunk[chunk[self.pval_col] < self.pval_thresh]

        if self.rsq_filter is not None:
          chunk = filter_imp_quality(chunk,self.rsq_col,self.rsq_filter)

        if self.query is not None:
          chunk = chunk.query(self.query)

        chunks.append(chunk)

      self.data = pd.concat(chunks)

    else:
      # We still need to run our filters, even if a data frame was passed in. 
      if self.pval_thresh is not None:
        self.data = self.data[self.data[self.pval_col] < self.pval_thresh]

      if self.rsq_filter is not None:
        self.data = filter_imp_quality(self.data,self.rsq_col,self.rsq_filter)

      if self.query is not None:
        self.data = self.data.query(self.query)

    for col in ('pval_col','marker_col','chrom_col','pos_col'):
      if self.__dict__[col] not in self.data.columns:
        error("column '%s' does not exist in your association results data - please set column names with --snp-col, --pval-col, etc." % self.__dict__[col]);

    # Drop SNPs that do not have chromosome
    self.data = self.data[self.data[self.chrom_col].notnull()];

    # Try to fix chromosome column. 
    # Should just be 1, 2, 3, X, Y and not chr4 or chrX
    self.data[self.chrom_col] = self.data[self.chrom_col].map(fix_chrom);

    # Drop SNPs that do not have a position.
    self.data = self.data[self.data[self.pos_col].notnull()];

    # Position must be an integer or it doesn't make sense.
    self.data[self.pos_col] = self.data[self.pos_col].astype('int');

    # Try to insert the trait as a column. 
    if 'TRAIT' in self.data.columns:
      warning("TRAIT column already exists in your association results file, can't put --trait into it!");
    else:
      self.data['TRAIT'] = self.trait;

  def liftover(self,build):
    pass

  # Drop a list of variants from the data. 
  # This will remove any variant with the same name OR the same chrom/pos. 
  def drop_variants(self,variant_list):
    for v in variant_list:
      # Drop any variant with the same name (marker ID.) 
      self.data = self.data[self.data[self.marker_col] != v.name];

      # Drop any variant with the same chrom/pos. 
      data_chrpos = self.data[self.chrom_col] + ":" + self.data[self.pos_col].map(lambda x: str(int(x)));
      variant_chrpos = v.as_chrpos();
      self.data = self.data[data_chrpos != variant_chrpos];

  # Keep a list of variants and discard the rest. 
  # Variants are kept if their names directly match, or they are at the same chrom/pos. 
  def keep_variants(self,variant_list):
    names = [v.name for v in variant_list];
    chrpos = [v.as_chrpos() for v in variant_list];

    data_chrpos = self.data[self.chrom_col] + ":" + self.data[self.pos_col].map(lambda x: str(int(x)));
    
    name_match = self.data[self.marker_col].isin(names);
    chrpos_match = data_chrpos.isin(chrpos);

    self.data = self.data[name_match | chrpos_match];

  def dist_clump(self,p_thresh=5e-08,dist=50E3):
    self.data = self.data[self.data[self.pval_col] < p_thresh];
    self.data = self.data.sort(self.pval_col);

    self.data['drop_row'] = False;
    begin_cols = [self.marker_col,self.pval_col,self.chrom_col,self.pos_col,'drop_row'];
    col_order = begin_cols + list(set(self.data.columns.tolist()).difference(begin_cols));
    self.data = self.data[col_order];

    if self.data.shape[0] <= 0:
      return;

    tree = ChromTree();
    for i in xrange(self.data.shape[0]):
      chrom = self.data.iat[i,2];
      pos = self.data.iat[i,3];

      if tree.find_overlap(chrom,pos - dist,pos + dist):
        self.data.iat[i,4] = True; # drop column
      else:
        tree.add_position(chrom,pos);

    final = self.data[self.data.drop_row == False];
    final = sort_genome(final,self.chrom_col,self.pos_col);
    del final['drop_row'];

    self.data = final;
    return final;

  def get_snps(self):
    snps = [];
    for row_index, row in self.data.iterrows():
      s = Variant();
      s.name = row[self.marker_col];
      s.chrom = row[self.chrom_col];
      s.pos = row[self.pos_col];
      s.chrpos = "%s:%s" % (s.chrom,s.pos);

      snps.append(s);

    return snps;

  def has_rows(self):
    if self.data.shape[0] <= 0:
      return False;
    else:
      return True;
