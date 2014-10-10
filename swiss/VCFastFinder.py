#!/usr/bin/env python
import gzip
import os
import subprocess
import sys
import hashlib
import re
from LDRegionCache import *
from textwrap import fill
from utils import *

class VCFastSettings:
  def __init__(self,vcf_path,vcfast_path):
    # Check each path. If it doesn't exist, try to find it relative
    # to the m2zfast root directory. 
    for arg,value in locals().items():
      if arg == 'self':
        continue;
      
      path = find_systematic(value);
      if path is None or not os.path.exists(path):
        if arg == "vcfast_path":
          error("cannot find vcfast - please set the path in the configuration file, or make it available on your PATH.");
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

    self.vcfast_path = vcfast_path;

  def createLDCacheKey(self):
    key_string = self.vcf_path;
    key = hashlib.sha512(key_string).hexdigest();
    return key;

class VCFastFinder():
  def __init__(self,vcfast_settings,cache=None,cleanup=True,verbose=False):
    if not isinstance(vcfast_settings,VCFastSettings):
      raise ValueError;
    
    self.data = {};
    self.snp = None;
    self.settings = vcfast_settings;
    self.debug = False;
    self.start = None;
    self.stop = None;
    self.chr = None;
    self.cache = cache;
    self.cleanup = cleanup;
    self.verbose = verbose;

    self.calc_ok = True;
    self.calc_fail = False;

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
        msg = fill("Error: could not find required file to generate LD using "
                  "VCFast: %s. Check your conf file to make sure paths are "
                  "correct. " % file);
        die(msg);

  def _run_sequence(self):
    self._check_geno_paths();
    data = self._run_vcfast();
  
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
    # Set by _run_vcfast
    return self.calc_ok;
  
  def _run_vcfast(self):
    if type(self.settings.vcf_path) is dict:
      vcf = self.settings.vcf_path.get(self.chr);
      if vcf is None:
        self.calc_ok = False;
        print >> sys.stderr, "Error: no VCF file available for chromosome '%s' - maybe a X vs 23 issue?" % self.chr;
        return {};
    else:
      vcf = self.settings.vcf_path;

    com = "{vcfast_path} index-LD --vcf {vcffile} --index {chrpos} --minMAF 0 --out - --win 99999999 --region {region} --min-r2 {min_r2}".format(
      vcfast_path = self.settings.vcfast_path,
      vcffile = vcf,
      chrpos = self.snp,
      region = "%s:%s-%s" % (self.chr,self.start,self.stop),
      min_r2 = self.min_r2
    );

    if self.verbose:
      print "Executing VCFast: %s" % com;

    proc = subprocess.Popen(com,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE);
    (stdout, stderr) = proc.communicate();

    if self.verbose:
      print stdout;
      print >> sys.stderr, stderr;

    if proc.returncode != 0:
      self.calc_ok = False;
      if not 'Cannot find a SNP at marker position' in stderr:
        print >> sys.stderr, "Error: VCFast did not complete successfully, errors were detected:";
        print >> sys.stderr, stderr;
  
    # Process LD info from stdout. 
    #CHROM  POS     ID      REF     ALT     AF      R2
    #10      94159   .       C       T       0.03619 0.00002

    start = False;
    data = {};
    for line in stdout.split("\n"):
      if line[0:6] == "#CHROM":
        start = True;
        continue;

      if start:
        e = line.split();
        if len(e) == 0:
          continue;

        chrpos = "%s:%s" % (e[0],e[1]);
        r2 = float(e[6]);

        if r2 > 1 or r2 < 0:
          self.calc_ok = False;
          raise ValueError, "Error calculating LD for variant %s, bad r2 value calculated (likely due to vcfast being used on unphased VCF), value was: %f" % (self.snp,r2);

        data[chrpos] = ("NA",r2);

    return data;

