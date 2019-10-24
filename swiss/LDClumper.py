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

from AssocResults import *
from utils import *
from Variant import *
import pandas as pd
from termcolor import *

class LDClumper:
  def __init__(self,assoc_results,finder):
    self.assoc = assoc_results
    self.finder = finder

  def ld_clump(self,ld_thresh=0.1,ld_dist=1E6):
    data = self.assoc.data

    variant_col = self.assoc.epacts_col
    pval_col = self.assoc.pval_col
    logp_col = self.assoc.logp_col
    chr_col = self.assoc.chrom_col
    pos_col = self.assoc.pos_col

    ld_dist = int(ld_dist)

    # Sort by p-value.
    data = data.sort_values(logp_col)

    print "\nStarting with significant %i SNPs:" % int(data.shape[0])
    preview_cols = filter(lambda x: x in data.columns, [variant_col,pval_col,logp_col,chr_col,pos_col])
    print sort_genome(data[preview_cols],chr_col,pos_col).to_string(index=False)
    print ""

    # Add a column to keep track of which variants were removed due to LD.
    data['ld_with'] = ""

    # Add a column to keep track of the LD values with those removed variants.
    data['ld_with_values'] = ""

    # Add a column to keep track of which variants we could not clump (don't exist in LD files, for example.)
    data['failed_clump'] = ""

    # Loop over variants, removing those in LD and keeping the most significant from each batch.
    failed_ld_variants = []
    for i in xrange(data.shape[0]):
      # Get the best variant at the top of the list.
      if i < data.shape[0]:
        top_chr = data[chr_col].iat[i]
        top_pos = data[pos_col].iat[i]
        top_marker = data[variant_col].iat[i]
      else:
        break

      # Calculate LD for this SNP within the region.
      print "Working on LD clumping variant %s" % top_marker

      ld_ok = self.finder.compute(top_marker,top_chr,top_pos - ld_dist,top_pos + ld_dist,ld_thresh)

      # Remove variants in LD with this one.
      if ld_ok:
        keep = [True] * data.shape[0]
        ld_values = []

        for j in xrange(i+1,data.shape[0]):
          row = data.iloc[j]
          row_marker = row[variant_col]
          ld_tuple = self.finder.data.get(row_marker) # returns (dprime, rsq)

          if ld_tuple is None:
            continue

          r2 = ld_tuple[1]

          if r2 >= ld_thresh:
            keep[j] = False
            ld_values.append(r2)

        # Which variants was the top SNP in LD with?
        data['ld_with'].iloc[i] = ",".join(data[variant_col][[not x for x in keep]])

        # What were their LD values?
        data['ld_with_values'].iloc[i] = ",".join(map(lambda v: "%0.2f" % v,ld_values))

        # This variant was successfully clumped
        data['failed_clump'].iloc[i] = "pass"

        # Remove variants in LD, move on to the next top SNP.
        data = data[keep]
      else:
        failed_variant = Variant()
        failed_variant.vid = data[variant_col].iat[i]
        failed_variant.chrom = data[chr_col].iat[i]
        failed_variant.pos = data[pos_col].iat[i]

        failed_ld_variants.append(failed_variant)

        data['failed_clump'].iloc[i] = "fail"
        print >> sys.stderr, colored("Warning: ",'yellow') + "could not calculate LD for variant %s at %s" % (top_marker,"%s:%s" % (top_chr,top_pos))

    data = sort_genome(data,chr_col,pos_col)
    self.assoc.data = data
    return self.assoc, failed_ld_variants

