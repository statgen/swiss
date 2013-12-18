#!/usr/bin/env python
from VCFastFinder import *
from AssocResults import *
from utils import *
from Variant import *
import pandas as pd
from termcolor import *

class LDClumper:
  def __init__(self,assoc_results,finder):
    self.assoc = assoc_results;
    self.finder = finder;

  def ld_clump(self,p_thresh=5e-08,ld_thresh=0.1,ld_dist=1E6):
    data = self.assoc.data;

    marker_col = self.assoc.marker_col;
    pval_col = self.assoc.pval_col;
    chr_col = self.assoc.chrom_col;
    pos_col = self.assoc.pos_col;

    ld_dist = int(ld_dist);
    
    # Start by selecting only those variants that are significant. 
    data = data[data[pval_col] < p_thresh];
    if data.shape[0] <= 0:
      return;

    # Sort by p-value. 
    data = data.sort(pval_col);

    print "\nStarting with significant %i SNPs:" % int(data.shape[0]);
    print sort_genome(data[[marker_col,pval_col,chr_col,pos_col]],chr_col,pos_col).to_string(index=False);
    print "";

    # Add a column to keep track of which variants were removed due to LD. 
    data['ld_with'] = "";

    # Add a column to keep track of which variants we could not clump (don't exist in LD files, for example.)
    data['failed_clump'] = "";

    # Loop over variants, removing those in LD and keeping the most significant from each batch.
    failed_ld_variants = [];
    for i in xrange(data.shape[0]):
      # Get the best variant at the top of the list. 
      if i < data.shape[0]:
        top_chr = data[chr_col].iat[i];
        top_pos = data[pos_col].iat[i];
      else:
        break;

      # Calculate LD for this SNP within the region. 
      ld_ok = self.finder.compute("%s:%s" % (top_chr,top_pos),top_chr,top_pos - ld_dist,top_pos + ld_dist,ld_thresh);

      # Remove variants in LD with this one. 
      if ld_ok:
        keep = [True] * data.shape[0];
        for j in xrange(i+1,data.shape[0]):
          row = data.iloc[j];
          row_chrpos = "%s:%s" % (row[chr_col],row[pos_col]);
          ld_value = self.finder.data.get(row_chrpos);

          if ld_value == None:
            continue;
          elif ld_value >= ld_thresh:
            keep[j] = False;

        # Which variants was the top SNP in LD with? 
        data['ld_with'].iloc[i] = ",".join(data[marker_col][[not x for x in keep]]);

        # This variant was successfully clumped
        data['failed_clump'].iloc[i] = "pass";

        # Remove variants in LD, move on to the next top SNP. 
        data = data[keep];
      else:
        failed_variant = Variant();
        failed_variant.name = data[marker_col].iat[i];
        failed_variant.chrom = data[chr_col].iat[i];
        failed_variant.pos = data[pos_col].iat[i];
        
        failed_ld_variants.append(failed_variant);

        data['failed_clump'].iloc[i] = "fail";
        print >> sys.stderr, colored("Warning: ",'yellow') + "could not calculate LD for variant %s - you should try a different source of LD information to properly clump these variants." % "%s:%s" % (top_chr,top_pos);

    data = sort_genome(data,chr_col,pos_col);
    self.assoc.data = data;
    return (self.assoc,failed_ld_variants);

