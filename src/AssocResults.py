#!/usr/bin/env python
import pandas as pd
from utils import *
from bx.intervals.intersection import *
from Variant import *

class ChromTree:
  def __init__(self):
    self.chrom = {};

  def add_position(self,chrom,pos):
    self.chrom.setdefault(chrom,IntervalTree()).add(pos,pos);

  def add_interval(self,chrom,start,end):
    self.chrom.setdefault(chrom,IntervalTree()).add(start,end);

  def find_overlap(self,chrom,start,end=None):
    node = self.chrom.get(chrom);
    if node == None:
      return None;
    else:
      if end == None:
        end = start;
      return node.find(start,end);

def sort_genome(dframe,chr_col,pos_col):
  dframe['int_chr'] = dframe[chr_col].map(chrom2chr);
  dframe = dframe.sort(['int_chr',pos_col]);
  del dframe['int_chr'];

  return dframe;

class AssocResults:
  def __init__(self,filepath=None,trait=None,df=None):
    self.filepath = filepath;
    self.data = df;
    self.trait = trait;

    self.pval_col = "P-value";
    self.marker_col = "MarkerName";
    self.chrom_col = "chr";
    self.pos_col = "pos";
    self.rsq_col = "RSQ";

  def load(self,*args,**kwargs):
    # Were we passed in a data frame that's already loaded, or a path to one?
    if self.data is None:
      # We don't have an already loaded data frame, so read it in.
      compr = 'gzip' if is_gzip(self.filepath) else None;
      self.data = pd.read_table(self.filepath,compression=compr,na_values=["NA","None","."],*args,**kwargs);

    for col in ('pval_col','marker_col','chrom_col','pos_col'):
      if self.__dict__[col] not in self.data.columns:
        error("column '%s' does not exist in your association results data - please set column names with --snp-col, --pval-col, etc." % self.__dict__[col]);

    # Try to fix chromosome column. 
    # Should just be 1, 2, 3, X, Y and not chr4 or chrX
    self.data[self.chrom_col] = self.data[self.chrom_col].map(lambda x: str(x).replace("chr",""));

    # Drop SNPs that do not have chromosome or position. 
    self.data = self.data[self.data[self.chrom_col].notnull()];
    
    self.data = self.data[self.data[self.pos_col].notnull()];
    self.data[self.pos_col] = self.data[self.pos_col].astype('int');

    # Try to insert the trait as a column. 
    if 'TRAIT' in self.data.columns:
      warning("TRAIT column already exists in your association results file, can't put --trait into it!");
    else:
      self.data['TRAIT'] = self.trait;

  def filter_imp_quality(self,threshold=0.3):
    try:
      threshold = float(threshold);
    except:
      error("Imputation quality threshold is not floatable, got: %s" % str(threshold));

    # RSQ needs to be float in order to threshold it correctly 
    self.data[self.rsq_col] = self.data[self.rsq_col].astype("float");

    start = self.data.shape[0];
    self.data = self.data[(self.data[self.rsq_col] >= threshold) | self.data[self.rsq_col].isnull()];
    end = self.data.shape[0];

    print "Imputation quality filter removed %i variants.." % (start - end);

  def do_filter(self,filter_expr):
    start = self.data.shape[0];
    self.data = self.data.query(filter_expr);
    end = self.data.shape[0];

    print "Filter '%s' removed %i variants.." % (filter_expr,start - end);

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
