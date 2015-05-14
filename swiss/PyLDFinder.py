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

import gzip
import os
import subprocess
import sys
import hashlib
import re
import time
import traceback
import numpy as np
from LDRegionCache import *
from textwrap import fill
from utils import *
from itertools import imap, chain, izip
from collections import defaultdict

# def rsquared(x, y):
#   # Drop observations that are None or NaN.
#   good_is = [];
#   for i in xrange(len(x)):
#     if (x[i] is not None) and (y[i] is not None):
#       good_is.append(i);
#
#   x = [x[i] for i in good_is];
#   y = [y[i] for i in good_is];
#
#   n = len(x)
#
#   sum_x = float(sum(x))
#   sum_y = float(sum(y))
#
#   sum_x_sq = sum(map(lambda x: pow(x, 2), x))
#   sum_y_sq = sum(map(lambda x: pow(x, 2), y))
#
#   psum = sum(imap(lambda x, y: x * y, x, y))
#   num = psum - (sum_x * sum_y/n)
#   den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
#
#   if den == 0: return 0
#
#   return pow(num / den,2)

# Slightly faster implementation.
def rsquared(x,y):
  x = np.asarray(x,dtype=np.float32);
  y = np.asarray(y,dtype=np.float32);

  not_nan = ~ (np.isnan(x) | np.isnan(y));

  x = x[not_nan];
  y = y[not_nan];

  return np.corrcoef(x,y)[0,1] ** 2;

def geno_to_code(g):
  g = g[0] + g[2];
  if '.' in g:
    return None;
  elif g == '00':
    return 0;
  elif g == '01':
    return 1;
  elif g == '10':
    return 1;
  elif g == '11':
    return 2;
  else:
    raise Exception, "Error: genotypes from VCF have more than 4 possible states (more than 2 ordered alleles.)";

def haplo_to_code(g):
  if g[1] != "|":
    raise Exception, "Error: unable to calculate D or D' since genotypes are not phased";

  try:
    g0 = int(g[0]);
  except:
    g0 = None;

  try:
    g1 = int(g[2]);
  except:
    g1 = None;

  if g0 is not None and g0 > 1:
    raise Exception, "Error: more than 2 alleles when trying to calculate LD!"

  if g1 is not None and g1 > 1:
    raise Exception, "Error: more than 2 alleles when trying to calculate LD!"

  return [g0,g1];

def ld_rsquare_indexsnp_vcf(index_pos,vcf_file,region,min_r2=0,tabix_path="tabix"):
  try:
    index_pos = int(index_pos);
  except:
    raise ValueError, "Index position was not an integer";

  # First grab the index SNP's genotypes.
  chrom = region.split(":")[0];
  if not int(index_pos) >= 0:
    raise IOError, "Error computing LD: index SNP position %s is invalid.." % str(index_pos);

  index_region = "{0}:{1}-{1}".format(chrom,index_pos);
  index_chrpos = "{0}:{1}".format(chrom,index_pos);

  p = subprocess.Popen([tabix_path,vcf_file,index_region],stdout=subprocess.PIPE,stderr=subprocess.PIPE);
  (stdout,stderr) = p.communicate();

  if stdout == '':
    raise IOError, "Error: while calculating LD from VCF file: index SNP position %s does not exist in file.. " % str(index_pos);

  if stderr != '':
    raise IOError, "Error: while calculating LD from VCF file: tabix generated an error: \n%s" % stderr;

  index_rec = None;

  # Sometimes tabix returns multiple lines at the same position. For example, if an indel at a different position
  # has enough alleles to overlap the requested position.
  for line in stdout.rstrip().split(os.linesep):
    index_ref, index_alt = line.split("\t")[3:5];
    if len(index_ref) > 1 or len(index_alt) > 1:
      continue;

    index_rec = line.rstrip().split("\t");

  if index_rec is None:
    raise IOError, "Error: couldn't calculate LD for index variant at position %s, no SNPs present at this position" % index_chrpos;

  try:
    index_gts = map(geno_to_code,index_rec[9:]);
  except:
    # Something was wrong with the variant's genotypes - possibly not biallelic
    traceback.print_exc();
    raise ValueError, "Error: while calculating LD from VCF file: index variant is not biallelic or had invalid genotypes";

  # Now grab the other SNPs, and calculate r2 with each of them.
  p = subprocess.Popen([tabix_path,vcf_file,region],stdout=subprocess.PIPE,stderr=subprocess.PIPE);
  (stdout,stderr) = p.communicate();

  if stderr != '':
    raise IOError, "Error: while calculating LD from VCF file: tabix generated an error: \n%s" % stderr;

  seen = {};
  markers = 0;
  ld_data = {};

  records = stdout.split("\n");
  for rec in records:
    if rec == '' or rec is None:
      continue;

    rec = rec.split("\t");
    rec[-1] = rec[-1].rstrip();

    # We'll only calculate LD with other SNPs and not indels.
    (rec_ref,rec_alt) = rec[3:5];

    if len(rec_ref) != 1:
     continue;
    elif len(rec_alt) != 1:
     continue;

    # Did it pass filters?
    rec_pass = rec[6];
    if rec_pass != "PASS":
      continue;

    # SNP/chr/pos
    (chr,pos,snp) = rec[0:3];
    chrpos = "{0}:{1}".format(chr,pos);

    if pos == index_pos:
      continue;

    # Have we seen this variant already?
    if seen.get((chr,pos)) is not None:
      # Skip repeated variants - they're probably indels and we only handle SNPs here.
      continue;
    else:
      seen[(chr,pos)] = 1;

    # Genotypes, converted to 0/1/2 coding
    try:
      gts = map(geno_to_code,rec[9:]);
    except:
      continue;

    # Calculate r2.
    rsq = rsquared(index_gts,gts);

    # sys.stdout = sys.__stdout__;
    # sys.stderr = sys.__stderr__;
    #
    # import IPython
    # IPython.embed();

    if (rsq < 0) or (rsq > 1):
      raise ValueError, "Error in LD calculation, bad r2 value: %f" % rsq;

    # Store LD calculations
    if rsq > min_r2:
      ld_data[chrpos] = ("NA",rsq);

    markers += 1;

  if markers == 0:
    raise IOError, "Error: no valid markers in VCF file to compute LD from..";

  return ld_data;

# Wrapper function for calculating LD.
# method can be 'rsquare' or 'dprime'
# Called like:
# ld_from_vcf('rsquare',index_pos,vcf_file,region,tabix_path='tabix')
def ld_from_vcf(method,*args,**kargs):
  if method == 'rsquare':
    return ld_rsquare_indexsnp_vcf(*args,**kargs);
  # elif method == 'dprime':
  #   return ld_dprime_indexsnp_vcf(*args,**kargs);
  else:
    raise Exception, "Error: only 'rsquare' and 'dprime' methods available for calculating LD from VCF file.";

# TODO: fix to have error handling similar to ld_rsquare_indexsnp_vcf, raise instead of return None
# def ld_dprime_indexsnp_vcf(index_pos,vcf_file,region,tabix_path="tabix"):
#   # First grab the index SNP's genotypes.
#   chrom = region.split(":")[0];
#   if not int(index_pos) >= 0:
#     print >> sys.stderr, "Error computing LD: index SNP position %s is invalid.." % str(index_pos);
#     return;
#
#   index_region = "{0}:{1}-{1}".format(chrom,index_pos);
#
#   if 'chr' not in chrom:
#     chrom = 'chr' + chrom;
#
#   index_chrpos = "{0}:{1}".format(chrom,index_pos);
#
#   p = subprocess.Popen([tabix_path,vcf_file,index_region],stdout=subprocess.PIPE,stderr=subprocess.PIPE);
#   (stdout,stderr) = p.communicate();
#
#   if stdout == '':
#     print >> sys.stderr, "Error: while calculating LD from VCF file: index SNP position %s does not exist in file.. " % str(index_pos);
#     return;
#
#   if stderr != '':
#     print >> sys.stderr, "Error: while calculating LD from VCF file: tabix generated an error: \n%s" % stderr;
#     return;
#
#   # Insert index SNP phased genotypes into list in order
#   index_rec = stdout.rstrip().split("\t");
#   try:
#     index_gts = [g for g in chain.from_iterable(map(haplo_to_code,index_rec[9:]))];
#   except Exception as e:
#     # If we're here, either:
#     # 1) the index SNP was not biallelic, or
#     # 2) the index SNP had genotypes that were not phased
#     print >> sys.stderr, e.message;
#     return None;
#
#   (index_ref,index_alt) = index_rec[3:5];
#
#   if len(index_ref) != 1:
#     print >> sys.stderr, "Error: while calculating LD from VCF file: index SNP is not a SNP - ref allele was %s, alt allele was %s" % (index_ref,index_alt);
#     return;
#
#   if len(index_alt) != 1:
#     print >> sys.stderr, "Error: while calculating LD from VCF file: index SNP is not a SNP - ref allele was %s, alt allele was %s" % (index_ref,index_alt);
#     return;
#
#   # Now grab the other SNPs, and calculate r2 with each of them.
#   p = subprocess.Popen([tabix_path,vcf_file,region],stdout=subprocess.PIPE,stderr=subprocess.PIPE);
#   (stdout,stderr) = p.communicate();
#
#   if stderr != '':
#     print >> sys.stderr, "Error: while calculating LD from VCF file: tabix generated an error: \n%s" % stderr;
#     return;
#
#   seen = {};
#   markers = 0;
#   ld_data = {};
#
#   records = stdout.split("\n");
#   for rec in records:
#     if rec == '' or rec is None:
#       continue;
#
#     rec = rec.split("\t");
#     rec[-1] = rec[-1].rstrip();
#
# #     Commented out for LD with biallelic indels
# #      # Is this a SNP?
# #      (rec_ref,rec_alt) = rec[3:5];
#
# #      if len(rec_ref) != 1:
# #        continue;
# #      elif len(rec_alt) != 1:
# #        continue;
#
#     # Did it pass filters?
#     rec_pass = rec[6];
#     if rec_pass != "PASS":
#       continue;
#
#     # SNP/chr/pos
#     (chr,pos,snp) = rec[0:3];
#     if 'chr' not in chr:
#       chr = 'chr' + chr;
#     chrpos = "{0}:{1}".format(chr,pos);
#
#     if pos == index_pos:
#       continue;
#
#     # Have we seen this variant already?
#     if seen.get((chr,pos)) is not None:
#       print >> sys.stderr, "Warning: multiple variants at same position (%s) in VCF file, using the last variant" % chrpos;
#     else:
#       seen[(chr,pos)] = 1;
#
#     # Phased genotypes inserted into array in order
#     try:
#       gts = [g for g in chain.from_iterable(map(haplo_to_code,rec[9:]))];
#     except:
#       # If this SNP is either not biallelic or has unphased genotypes, skip it
#       continue;
#
#     # Number of non-missing alleles
#     nonmiss_gts = sum(map(lambda x: x is not None,gts));
#
#     # Calculate statistics from haplos
#     counts = defaultdict(int);
#     nonmiss_haplos = 0.0;
#     index_afs = [0,0];
#     afs = [0,0];
#     for haplo in izip(index_gts,gts):
#       if None in haplo:
#         continue;
#
#       counts[haplo] += 1;
#       index_afs[haplo[0]] += 1;
#       afs[haplo[1]] += 1;
#       nonmiss_haplos += 1;
#
#     # Normalize allele frequencies
#     index_afs = map(lambda x: x / nonmiss_haplos,index_afs);
#     afs = map(lambda x: x / nonmiss_haplos,afs);
#
#     # Normalize haplo counts
#     hfreq = defaultdict(int);
#     for k,v in counts.iteritems():
#       hfreq[k] = v / nonmiss_haplos;
#
#     # D statistic
#     d_stat = (hfreq[(0,0)]*hfreq[(1,1)] - hfreq[(1,0)]*hfreq[(0,1)]);
#
#     # D' statistic
#     if d_stat < 0:
#       d_max = max(-1 * index_afs[0] * afs[0], -1 * index_afs[1] * afs[1]);
#       d_prime = d_stat / d_max;
#     elif d_stat > 0:
#       d_max = min(index_afs[0] * afs[1], index_afs[1] * afs[0]);
#       d_prime = d_stat / d_max;
#     else:
#       d_prime = d_stat;
#
#     # Write out in the format expected by locuszoom.
#     ld_data[chrpos] = (d_prime,"NA");
#
#     markers += 1;
#
#   if markers == 0:
#     print >> sys.stderr, "Error: no valid markers in VCF file to compute LD from..";
#     return;
#
#   return ld_data;

class PyLDSettings:
  def __init__(self,vcf_path,tabix_path):
    # Check each path. If it doesn't exist, try to find it relative
    # to the m2zfast root directory. 
    for arg,value in locals().items():
      if arg == 'self':
        continue;
      
      path = find_systematic(value);
      if path is None or not os.path.exists(path):
        if arg == "tabix_path":
          error("cannot find tabix - please set the path in the configuration file, or make it available on your PATH.");
        else:
          error("path either does not exist or insufficient permissions to access it: %s" % str(value));
      else:
        exec "%s = \"%s\"" % (arg,path);
    
    if '.json' in vcf_path:
      import json
      with open(vcf_path) as jsin:
        self.vcf_path = json.load(jsin);
    else:
      self.vcf_path = vcf_path;

    self.tabix_path = tabix_path;

  def create_ld_cache_key(self):
    key_string = self.vcf_path;
    key = hashlib.sha512(key_string).hexdigest();
    return key;

class PyLDFinder():
  def __init__(self,pyld_settings,cache=None,cleanup=True,verbose=False):
    if not isinstance(pyld_settings,PyLDSettings):
      raise ValueError;
    
    self.data = {};
    self.snp = None;
    self.settings = pyld_settings;
    self.debug = False;
    self.start = None;
    self.stop = None;
    self.chr = None;
    self.cache = cache;
    self.cleanup = cleanup;
    self.verbose = verbose;

    self.calc_ok = True;

    self.min_r2 = None;

  def write(self,filename):
    try:
      f = open(filename,'w');
      print >> f, "snp1 snp2 dprime rsquare";
  
      if len(self.data) == 0:
        return False;
  
      for snp in self.data:
        print >> f, "%s %s %s %s" % (
          snp,
          self.snp,
          str(self.data.get(snp)[0]),
          str(self.data.get(snp)[1])
        );
  
      f.close();
    except:
      print >> sys.stderr, "Error: could not write computed LD to disk, permissions?";
      return False;
    
    return True;

  def _check_geno_paths(self):
    files = [];
    if type(self.settings.vcf_path) is dict:
      map(files.append,self.settings.vcf_path.itervalues());
    else:
      files.append(self.settings.vcf_path);
    
    for file in files:
      if not os.path.exists(file):
        msg = fill("Error: could not find required file to generate LD "
                  "%s. Check your conf file to make sure paths are "
                  "correct. " % file);
        die(msg);

  def _run_sequence(self):
    self._check_geno_paths();
    data = self._run_pyld();
  
    return data;

  def compute(self,snp,chr,start,stop,min_r2=0):
    self.snp = snp;
    self.start = max(0,start);
    self.stop = stop;
    self.chr = chr;
    self.data = None;
    self.min_r2 = min_r2;
    self.calc_ok = True;

    # If the cache has data for this SNP and region, use it.
    # Otherwise, compute it.
    if self.cache:
      if self.cache.hasRegion(snp,start,stop):
        self.data = self.cache.getAllLD(snp);
      else:
        self.data = self._run_sequence();
        self.cache.updateLD(snp,start,stop,self.data);
    else:
      self.data = self._run_sequence();

    # Complete successfully?
    return self.calc_ok;
  
  def _run_pyld(self):
    if type(self.settings.vcf_path) is dict:
      vcf = self.settings.vcf_path.get(self.chr);
      if vcf is None:
        self.calc_ok = False;
        print >> sys.stderr, "Error: no VCF file available for chromosome '%s' - maybe a X vs 23 issue?" % self.chr;
        return {};
    else:
      vcf = self.settings.vcf_path;

    region = "%s:%s-%s" % (self.chr,self.start,self.stop);
    index_pos = self.snp.split(":")[1];

    data = {};
    try:
      data = ld_from_vcf("rsquare",index_pos,vcf,region,min_r2=self.min_r2,tabix_path=self.settings.tabix_path);
    except Exception as e:
      if SWISS_DEBUG: raise

      print >> sys.stderr, e.message;
      self.calc_ok = False;

    return data;
