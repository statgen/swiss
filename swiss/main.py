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

import copy
import os.path
import atexit
import __builtin__
import codecs
import shlex
import traceback
import difflib
import pandas as pd
import appdirs
from optparse import *
from termcolor import *
from glob import glob
from .utils import *
from PlinkLDFinder import PlinkLDFinder, PlinkLDSettings
from LDClumper import *
from AssocResults import *
from GWASCatalog import *
from itertools import *
from VerboseParser import *
from multiprocessing import Pool, cpu_count
from pprint import pprint
from textwrap import wrap
import swiss.conf
from swiss.conf.reader import read_conf
from swiss.conf.writer import write_conf

PROG_NAME = "Swiss"
PROG_VERSION = "1.0.0"
PROG_DATE = "11/26/2016"
PROG_AUTHOR = "Ryan Welch (welchr@umich.edu)"
PROG_URL = "https://github.com/welchr/Swiss"

__builtin__.SWISS_DEBUG = False

EPACTS_DTYPES = {
  "CALLRATE" : pd.np.float32,
  "MAF" : pd.np.float32,
  "BEG" : pd.np.float64,
  "END" : pd.np.float64
}

def find_likely_file(filepath):
  orig = os.path.split(filepath)[1]
  dirp = os.path.dirname(filepath)
  files = os.listdir(dirp)

  matches = difflib.get_close_matches(orig,files,1,0.95)
  if len(matches) > 0:
    return os.path.join(dirp,matches[0])
  else:
    return None

def show_diff(str1,str2):
  print "\n".join(difflib.ndiff([str1],[str2]))

class BasePair:
  def __init__(self,bp):
    self.bp = bp

  def as_kb(self):
    kb = float(self.bp) / 1000.0
    if int(kb) == kb:
      kb = int(kb)

    return "%ikb" % kb

def get_user_config_path():
  from os import path as p
  return p.join(appdirs.user_config_dir(),"swiss.yaml")

def get_conf(fpath=None):
  from os import path as p
  conf_file = fpath

  if fpath is None:
    # First try in user config dir. 
    test = get_user_config_path()
    if p.isfile(test):
      conf_file = test

    if conf_file is None:
      # Okay, wasn't there. Now try loading the default. 
      default = p.join(swiss.conf.__path__[0],"default.yaml")
      if p.isfile(default):
        conf_file = default

  # By now we should have found it. 
  if conf_file is not None:
    config = read_conf(conf_file)
    config["config_path"] = conf_file
  else:
    raise IOError("Cannot locate swiss configuration file")

  return config

def print_program_header():
  prog_string = "%s - %s (%s)" % (PROG_NAME,PROG_VERSION,PROG_DATE)
  max_length = max(map(len,[prog_string,PROG_AUTHOR,PROG_URL]))
  max_length += 4 # buffer

  # Top line
  def print_table_line(width=max_length):
    print "".join(chain("+",repeat('-',max_length),"+"))

  print_table_line()

  # Rows
  print "|%s|" % str.center(prog_string,max_length)
  print_table_line()
  print "|%s|" % str.center(PROG_AUTHOR,max_length)
  print_table_line()
  print "|%s|" % str.center(PROG_URL,max_length)
  print_table_line()

def find_gwas_catalog(conf,build,gwas_cat):
  import os.path as p

  entry = conf["gwas_catalogs"][build][gwas_cat]
  if p.isfile(entry):
    return entry
  else:
    # Relative to package install location
    import swiss
    default = p.join(swiss.__path__[0],entry)
    
    if p.isfile(default):
      return default
    else:
      raise IOError, "Can't find catalog: " + entry

def find_ld_source(conf,build,ld_source):
  import os.path as p

  entry = conf["ld_sources"][build][ld_source]
  if p.isfile(entry):
    return entry
  else:
    # Relative to package install location
    import swiss
    default = p.join(swiss.__path__[0],entry)
    
    if p.isfile(default):
      return default
    else:
      raise IOError, "Can't find LD source: " + entry

def get_settings(arg_string=None):
  usage = "swiss [options]"
  parser = VerboseParser(usage=usage)

  parser.add_option("--list-files",help="Show the locations of files in use by swiss.",default=False,action="store_true")

  # Association result input options. 
  parser.add_option("--assoc",help="[Required] Association results file.")
  parser.add_option("--multi-assoc",help="Designate that the results file is in EPACTS multi-assoc format.",action="store_true",default=False)
  parser.add_option("--trait",help="Description of phenotype for association results file. E.g. 'HDL' or 'T2D'")
  parser.add_option("--delim",help="Association results delimiter.",default="\t")
  parser.add_option("--build",help="Genome build your association results are anchored to.",default="hg19")
  parser.add_option("--variant-col",help="Variant column name in results file.",default="MARKER_ID")
  parser.add_option("--pval-col",help="P-value column name in results file.",default="PVALUE")
  parser.add_option("--chrom-col",help="Chromosome column name in results file.",default="CHR")
  parser.add_option("--pos-col",help="Position column name in results file.",default="POS")
  parser.add_option("--rsq-col",help="Imputation quality column name.",default="RSQ")
  parser.add_option("--trait-col",help="Trait column name. Can be omitted, in which case the value of --trait will be added as a column.",default=None)

  # Association result filtering results. 
  parser.add_option("--rsq-filter",help="Remove variants below this imputation quality.",default=None)
  parser.add_option("--filter",help="Give a general filter string to filter variants.",default=None)

  # Output options. 
  parser.add_option("--out",help="Prefix for output files.",default="swiss_output")

  # LD clumping options. 
  parser.add_option("--ld-clump",help="Clump association results by LD.",action="store_true",default=False)
  parser.add_option("--clump-p",help="P-value threshold for LD and distance based clumping.",default=5e-08)
  parser.add_option("--clump-ld-thresh",help="LD threshold for clumping.",default=0.2,type="float")
  parser.add_option("--clump-ld-dist",help="Distance from each significant result to calculate LD.",default=1E6,type="int")

  # Distance clumping options. 
  parser.add_option("--dist-clump",help="Clump association results by distance.",action="store_true",default=False)
  parser.add_option("--clump-dist",help="Distance threshold to use for clumping based on distance.",default=2.5e5,type="int")

  # LD source (GoT2D, 1000G, etc.) 
  parser.add_option("--ld-clump-source",help="Name of pre-configured LD source, or a VCF file from which to compute LD.",default="1000G_2012-03_EUR")
  parser.add_option("--list-ld-sources",help="Print a list of available LD sources for each genome build.",default=False,action="store_true")

  # GWAS catalog
  parser.add_option("--gwas-cat",help="GWAS catalog to use.",default="ebi")
  parser.add_option("--ld-gwas-source",help="Name of pre-configured LD source or VCF file to use when calculating LD with GWAS variants.",default="1000G_2012-03_EUR")
  parser.add_option("--list-gwas-cats",action="store_true",default=False,help="Give a listing of all valid GWAS catalogs and their descriptions.")
  parser.add_option("--list-gwas-traits",action="store_true",default=False,help="List all of the available traits in a selected GWAS catalog.")
  parser.add_option("--list-gwas-trait-groups",action="store_true",default=False,help="List all of the available groupings of traits in a selected GWAS catalog.")
  parser.add_option("--gwas-cat-p",help="P-value threshold for GWAS catalog variants.",default=5e-08,type="float")
  parser.add_option("--gwas-cat-ld",help="LD threshold for considering a GWAS catalog variant in LD.",default=0.1,type="float")
  parser.add_option("--gwas-cat-dist",help="Distance threshold for considering a GWAS catalog variant 'nearby'.",type="int",default=2.5e5)
  parser.add_option("--include-cols",help="List of columns to merge in from association results (grouped by variant.)",default=None)
  parser.add_option("--skip-overlap-check",help=SUPPRESS_HELP,default=False,action="store_true")
  parser.add_option("--do-overlap-check",help="Perform the check of whether the GWAS catalog has variants that are not in your --ld-gwas-source.",default=False,action="store_true")
  parser.add_option("--skip-gwas",help="Skip the step of looking for GWAS hits in LD with top variants after clumping.",default=False,action="store_true")

  # LD cache
  parser.add_option("--cache",help="Prefix for LD cache.",default="ld_cache")

  # Misc options
  parser.add_option("-T","--threads",default=1,type="int",help="Number of parallel jobs to run. Only works with --multi-assoc currently.")
  parser.add_option("--debug",default=False,action="store_true",help=SUPPRESS_HELP)

  if arg_string is None:
    (opts,args) = parser.parse_args()
  else:
    (opts,args) = parser.parse_args(shlex.split(arg_string))

  conf = get_conf()

  if opts.list_files:
    print "Configuration file was found at: " + conf["config_path"]
    print "The config file can be overridden by copying the default.yaml file to " + get_user_config_path()
    print "The default swiss data directory is: " + os.path.join(swiss.__path__[0],"data/")
    sys.exit(0)

  if opts.debug:
    __builtin__.SWISS_DEBUG = True
    pd.set_option('chained_assignment','warn')
  else:
    pd.set_option('chained_assignment',None)

  if SWISS_DEBUG:
    print >> sys.stderr, "Package paths in use: "
    for x in sys.path:
      print >> sys.stderr, x
    print >> sys.stderr, ""

  if opts.threads < 1:
    opts.threads = 1

  if opts.threads > cpu_count():
    print >> sys.stderr, "Warning: you set threads to %i, but the CPU count for this machine is %i..." % (opts.threads,cpu_count())

  if opts.list_gwas_cats:
    print_gwas()
    sys.exit(0)

  if not os.path.isfile(opts.gwas_cat):
    opts.gwas_cat_file = find_gwas_catalog(conf,opts.build,opts.gwas_cat)

    if opts.gwas_cat_file is None:
      error("Could not locate GWAS catalog file for conf entry '%s - %s'!" % (opts.build,opts.gwas_cat))
  else:
    opts.gwas_cat_file = opts.gwas_cat

  if opts.gwas_cat_file is not None:
    gcat = GWASCatalog(opts.gwas_cat_file)

    if opts.list_gwas_traits:
      print "Available traits for GWAS catalog '%s': \n" % opts.gwas_cat
      for group, traits in gcat.get_trait_group_pairs():
        print group
        print "".join(['-' for _ in xrange(len(group))])
        print ""

        for t in traits:
          print t

        print ""

      sys.exit(0)

    elif opts.list_gwas_trait_groups:
      print "Available trait groups for GWAS catalog '%s': " % opts.gwas_cat
      for group in sorted(gcat.get_trait_groups()):
        print group

      sys.exit(0)

  if opts.list_ld_sources:
    print_ld_sources()
    sys.exit(0)

  if opts.assoc is None:
    parser.print_help()
    print ""
    print >> sys.stderr, "Must specify association results file: --assoc"
    sys.exit(1)

  if not os.path.isfile(opts.assoc):
    parser.print_help()
    print ""
    print >> sys.stderr, "Association results file does not exist: %s" % opts.assoc
    sys.exit(1)

  if opts.ld_clump and opts.dist_clump:
    print >> sys.stderr, "Error: can't specify --ld-clump and --dist-clump together, choose one method to clump results"
    sys.exit(1)

  if opts.ld_clump:
    opts.clump_ld_thresh = float(opts.clump_ld_thresh)
    if opts.clump_ld_thresh < 0 or opts.clump_ld_thresh > 1:
      error("LD threshold for clumping must be >= 0 or <= 1.")

  if opts.clump_p in ("None","NA","NULL",""):
    opts.clump_p = None

  if opts.clump_p is not None:
    opts.clump_p = float(opts.clump_p)
    if opts.clump_p < 0 or opts.clump_p > 1:
      error("P-value threshold must be >= 0 or <= 1.")
  else:
    warning("No p-value threshold was specified! This will result in the best SNP being randomly selected.")

  opts.clump_ld_dist = int(float(opts.clump_ld_dist))

  # If one source is specified, but not the other, assume the user meant to use this for both. 
  # ld clump source XOR ld gwas source
  if (opts.ld_clump_source is None) != (opts.ld_gwas_source is None):
    if opts.ld_clump_source is None:
      opts.ld_clump_source = opts.ld_gwas_source
    elif opts.ld_gwas_source is None:
      opts.ld_gwas_source = opts.ld_clump_source
    
  # LD clumping source
  if not os.path.isfile(opts.ld_clump_source):
    if not conf["ld_sources"][opts.build].has_key(opts.ld_clump_source):
      # They specified something, but it's apparently not a file, and not a key. 
      # Maybe the file is a typo? 
      match = find_likely_file(opts.ld_clump_source)
      if match is not None:
        print "Couldn't find --ld-clump-source, but a file in the directory was a possible match: "
        print show_diff(opts.ld_clump_source,match)
        sys.exit(1)

    opts.ld_clump_source_file = find_ld_source(conf,opts.build,opts.ld_clump_source)

    if opts.ld_clump_source_file is None:
      error("Could not locate VCF file for conf entry '%s - %s'!" % (opts.build,opts.ld_clump_source))
  else:
    opts.ld_clump_source_file = opts.ld_clump_source
  
  # GWAS LD lookup source
  if not os.path.isfile(opts.ld_gwas_source):
    if not conf["ld_sources"][opts.build].has_key(opts.ld_gwas_source):
      # They specified something, but it's apparently not a file, and not a key. 
      # Maybe the file is a typo? 
      match = find_likely_file(opts.ld_gwas_source)
      if match is not None:
        print "Couldn't find --ld-gwas-source, but a file in the directory was a possible match: "
        print show_diff(opts.ld_gwas_source,match)
        sys.exit(1)
    
    opts.ld_gwas_source_file = find_ld_source(conf,opts.build,opts.ld_gwas_source)

    if opts.ld_gwas_source_file is None:
      error("Could not locate VCF file for conf entry '%s - %s'!" % (opts.build,opts.ld_gwas_source))
  else:
    opts.ld_gwas_source_file = opts.ld_gwas_source

  if not os.path.isfile(opts.assoc):
    error("Could not locate association results file (or insufficient permissions): %s" % opts.assoc)

  # File delimiter
  if opts.delim in ("tab","\t","\\t"):
    opts.delim = "\t"
  elif opts.delim in ("comma",","):
    opts.delim = ","
  elif opts.delim in ("space"," "):
    opts.delim = " "
  elif opts.delim == "whitespace":
    opts.delim = None

  try:
    opts.vcfast_path = find_systematic(conf["vcfast_path"])
  except:
    # TODO: if we ever use vcfast again, this will need to be patched
    pass

  opts.tabix_path = find_systematic(conf["tabix_path"])
  opts.plink_path = find_systematic(conf["plink_path"])

  out_exists = glob(os.path.join(opts.out,"*"))
  if len(out_exists) > 0:
    error("Output files already exist with this prefix: %s" % opts.out)

  # If multi-assoc is specified, the column names are already known.
  if opts.multi_assoc:
    opts.variant_col = "MARKER_ID"
    opts.pval_col = "PVALUE"
    opts.chrom_col = "#CHROM"
    opts.pos_col = "BEG"

  return (opts,args)

def print_gwas():
  conf = get_conf()

  print "%-10s %-40s" % ("Build","Catalog")
  print "%-10s %-40s" % ("-----","-------")
  for build, cats in conf["gwas_catalogs"].iteritems():
    print "%-10s %-40s" % (build,",".join(cats.keys()))

def print_ld_sources():
  conf = get_conf()

  print "%-10s %-40s" % ("Build","LD Sources")
  print "%-10s %-40s" % ("-----","----------")
  for build, ld_source in conf["ld_sources"].iteritems():
    ld_source_names = sorted(ld_source)
    print "%-10s %-40s" % (build,", ".join(ld_source_names))

class StreamTee(object):
  # Based on https://gist.github.com/327585 by Anand Kunal
  def __init__(self, stream1, stream2):
    self.stream1 = stream1
    self.stream2 = stream2
    self.__missing_method_name = None # Hack!

  def __getattribute__(self, name):
    return object.__getattribute__(self, name)

  def __getattr__(self, name):
    self.__missing_method_name = name # Could also be a property
    return getattr(self, '__methodmissing__')

  def __methodmissing__(self, *args, **kwargs):
    # Emit method call to the log copy
    callable2 = getattr(self.stream2, self.__missing_method_name)
    callable2(*args, **kwargs)

    # Emit method call to stdout (stream 1)
    callable1 = getattr(self.stream1, self.__missing_method_name)
    return callable1(*args, **kwargs)

def get_header(infile,sep="\t"):
  if infile.endswith(".gz"):
    f = gzip.open(infile)
  else:
    f = open(infile)

  with f:
    h = f.readline().split(sep)
    h[-1] = h[-1].rstrip()

    return h

# # Helper function to iterate over separate traits from an EPACTS multi assoc file.
# def multiassoc_epacts_iter(result_file):
#   header = get_header(result_file)
#
#   #CHROM    BEG    END       MARKER_ID    NS        AC  CALLRATE   GENOCNT      MAF  DHA.P   DHA.B  EstC.P  EstC.B  FAw3.P  FAw3.B  FAw3toFA.P  FAw3toFA.B  FAw6.P  FAw6.B  FAw6toFA.P
#
#   trait_ps = filter(lambda x: x.endswith(".P"),header)
#   trait_betas = filter(lambda x: x.endswith(".B"),header)
#   trait_cols = trait_ps + trait_betas
#   intro_cols = filter(lambda x: x not in trait_cols,header)
#
#   # Load first columns, since we'll always use them, along with the first trait column.
#   # base_trait = first trait to load
#   base_trait_p = trait_ps[0]
#   base_trait_b = trait_betas[0]
#   base_trait = base_trait_p.replace(".P","")
#   base_cols = intro_cols + [base_trait_p,base_trait_b]
#
#   print "\nLoading trait %s from association results file: %s" % (base_trait,result_file)
#
#   base_df = pd.read_table(result_file,compression = "gzip" if result_file.endswith(".gz") else None,na_values=["NA","None","."],usecols = base_cols)
#   base_df.rename(columns = {base_trait_p : "PVALUE",base_trait_b : "BETA"},inplace=True)
#
#   yield (base_trait,base_df)
#
#   for p in trait_ps[1:]:
#     trait = p.replace(".P","")
#     bcol = trait + ".B"
#
#     print "\nLoading trait %s from association results file: %s" % (trait,result_file)
#
#     df = pd.read_table(result_file,compression = "gzip" if result_file.endswith(".gz") else None,na_values=["NA","None","."],usecols = [p,bcol])
#     df.rename(columns = {p : "PVALUE",bcol : "BETA"},inplace=True)
#
#     base_df["PVALUE"] = df["PVALUE"]
#     base_df["BETA"] = df["BETA"]
#
#     yield (trait,base_df)

def multiassoc_epacts_iter(result_file,trait,pval_thresh=None,rsq_col=None,rsq_filter=None,query=None):
  header = get_header(result_file)

  traits = []
  for h in header:
    if h.endswith(".P"):
      traits.append(h.replace(".P",""))

  for trait in traits:
    dframe = multiassoc_epacts_load(result_file,trait,pval_thresh,rsq_col,rsq_filter,query)
    dframe.rename(columns = {
      trait + ".P" : "PVALUE",
      trait + ".B" : "BETA"
    },inplace=True)

    yield trait, dframe

def multiassoc_epacts_load(result_file,trait,pval_thresh=None,rsq_col=None,rsq_filter=None,query=None,verbose=True):
  if verbose:
    print "\nLoading trait %s from association results file: %s" % (trait,result_file)

  header = get_header(result_file)

  #CHROM    BEG    END       MARKER_ID    NS        AC  CALLRATE   GENOCNT      MAF  DHA.P   DHA.B  EstC.P  EstC.B  FAw3.P  FAw3.B  FAw3toFA.P  FAw3toFA.B  FAw6.P  FAw6.B  FAw6toFA.P

  # Check to make sure the trait requested is in the header.
  trait_ps = filter(lambda x: x.endswith(".P"),header)
  trait_betas = filter(lambda x: x.endswith(".B"),header)
  trait_cols = trait_ps + trait_betas
  intro_cols = filter(lambda x: x not in trait_cols,header)

  trait_names = map(lambda x: x.replace(".P",""),trait_ps)
  if trait not in trait_names:
    raise IOError, "Requested trait %s is not present in %s" % (trait,result_file)

  this_trait_cols = [trait + ".P",trait + ".B"]
  pval_col = trait + ".P"
  beta_col = trait + ".B"

  dtypes = EPACTS_DTYPES.copy()
  dtypes[beta_col] = pd.np.float32

  df_iter = pd.read_table(result_file,
    compression = "gzip" if result_file.endswith(".gz") else None,
    na_values=["NA","None","."],
    usecols = intro_cols + this_trait_cols,
    iterator = True,
    chunksize = 500000,
    dtype = dtypes
  )

  chunks = []
  for chunk in df_iter:
    if pval_thresh is not None:
      chunk = chunk[chunk[pval_col] < pval_thresh]

    if rsq_filter is not None:
      chunk = filter_imp_quality(chunk,rsq_col,rsq_filter)

    if query is not None:
      chunk = chunk.query(query)

    chunks.append(chunk)

  df = pd.concat(chunks)

  df.rename(columns = {
    trait + ".P" : "PVALUE",
    trait + ".B" : "BETA"
  },inplace=True)

  return df

def multiassoc_epacts_get_traits(result_file):
  header = get_header(result_file)

  # Check to make sure the trait requested is in the header.
  trait_ps = filter(lambda x: x.endswith(".P"),header)
  traits = map(lambda x: x.replace(".P",""),trait_ps)

  return traits

def merge_include_cols_gwas_hits(gwas_hits,results,include_cols,variant_col):
  include_cols = [i.strip() for i in include_cols.split(",")]
  include_cols = filter(lambda x: x in results.data.columns,include_cols)

  if len(include_cols) == 0:
    print "Warning: user specified --include-cols, but none of them existed in the association results!"
  else:
    assoc_incl_cols = results.data[[variant_col] + include_cols]
    assoc_incl_cols.rename(
      columns = dict(zip(include_cols,map(lambda x: "ASSOC_" + x,include_cols))),
      inplace = True
    )
    gwas_hits = pd.merge(gwas_hits,assoc_incl_cols,left_on="ASSOC_VARIANT",right_on=variant_col)
    del gwas_hits[variant_col]

  return gwas_hits

def fprint(str):
  print str

def run_process(assoc,trait,outprefix,opts):
  if isinstance(assoc,str):
    results = AssocResults(assoc,trait,pval_thresh=opts.clump_p,rsq_filter=opts.rsq_filter,query=opts.filter)
  else:
    results = AssocResults(trait=trait,df=assoc,pval_thresh=opts.clump_p,rsq_filter=opts.rsq_filter,query=opts.filter)

  results.vid_col = opts.variant_col
  results.chrom_col = opts.chrom_col
  results.pos_col = opts.pos_col
  results.pval_col = opts.pval_col
  results.rsq_col = opts.rsq_col
  if opts.trait_col is not None:
    results.trait_col = opts.trait_col
  results.load(sep=opts.delim)

  # LD finder for clumping
  vset = PlinkLDSettings(opts.ld_clump_source_file,opts.tabix_path,opts.plink_path)
  finder_clumping = PlinkLDFinder(vset,verbose=False,cache=None)

  # LD finder for GWAS catalog lookups
  vset_gwas = PlinkLDSettings(opts.ld_gwas_source_file,opts.tabix_path,opts.plink_path)
  finder_gwas = PlinkLDFinder(vset_gwas,verbose=False,cache=None)

  # GWAS catalog.
  gcat = GWASCatalog(opts.gwas_cat_file)

  if opts.do_overlap_check:
    print "\nIdentifying GWAS catalog variants that do not overlap with your --ld-gwas-source: %s" % opts.ld_gwas_source
    missing_vcf = gcat.variants_missing_vcf(opts.ld_gwas_source_file)
    missing_vcf.sort('PHENO',inplace=True)
    missing_vcf = sort_genome(missing_vcf,'CHR','POS')
    print colored('Warning: ','yellow') + "the following variants in the GWAS catalog are not present in your VCF file: "
    print missing_vcf["VARIANT EPACTS CHR POS PHENO GROUP".split()].to_string(index=False)
  else:
    print ""
    map(fprint,wrap("Skipping check of whether GWAS catalog variants are missing from your LD source. To enable, "
                    "use --do-overlap-check. Note this can take a fair amount of time depending on the "
                    "size of your VCF(s)."
    ))

  print "\nLoaded %i variants (after filters) from association results.." % results.data.shape[0]

  if opts.ld_clump:
    print "\nLD clumping results.."
    print "\nLD source: %s" % opts.ld_clump_source

    clumper = LDClumper(results,finder_clumping)
    clump_results = clumper.ld_clump(opts.clump_ld_thresh,opts.clump_ld_dist)

    if clump_results is not None:
      results_clumped, failed_ld_clump_variants = clump_results
    else:
      print "\nNo significant variants available for LD clumping, skipping.."
      return

    print "\nResults after clumping: "
    print_cols = [
      opts.variant_col,
      opts.chrom_col,
      opts.pos_col,
      opts.pval_col,
      'ld_with',
      'failed_clump'
    ]
    print results_clumped.data[print_cols].to_string(index=False)

    out_clump = outprefix + ".clump"

    print "\nWriting clumped results to: %s" % out_clump
    results_clumped.data.to_csv(out_clump,index=False,sep="\t",na_rep="NA")

    if not opts.skip_gwas:
      print "\nFinding clumped results in LD with GWAS catalog variants..."
      print "\nLD source: %s" % opts.ld_gwas_source
      gwas_ld, gwas_ld_failed_variants = gcat.variants_in_ld_multi(results_clumped,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist,opts.threads)

      if gwas_ld is not None:
        # If the user requested other columns be merged in with the gwas_ld, pull 'em out.
        if opts.include_cols:
          gwas_ld = merge_include_cols_gwas_hits(gwas_ld,results_clumped,opts.include_cols,opts.variant_col)

        out_ld_gwas = outprefix + ".ld-gwas.tab"
        print "\nWriting GWAS catalog variants in LD with clumped variants to: %s" % out_ld_gwas
        gwas_ld.to_csv(out_ld_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits were in LD with clumped variants."

      gwas_near = gcat.variants_nearby(results_clumped,opts.gwas_cat_dist)

      if gwas_near is not None:
        # If the user requested other association results columns be merged in with the gwas_hits, add them in.
        if opts.include_cols:
          gwas_near = merge_include_cols_gwas_hits(gwas_near,results_clumped,opts.include_cols,opts.variant_col)

        out_near_gwas = outprefix + ".near-gwas.tab"
        print "Writing GWAS catalog variants within %s of a clumped variant to: %s" % (BasePair(opts.gwas_cat_dist).as_kb(),out_near_gwas)
        gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits discovered within %s of any clumped results." % BasePair(opts.gwas_cat_dist).as_kb()

  elif opts.dist_clump:
    dist_ok = results.dist_clump(opts.clump_p,opts.clump_dist)

    if dist_ok is None:
      print "\nNo significant variants available for distance clumping, skipping.."
      return

    print "\nResults after clumping: "
    print_cols = [
      opts.variant_col,
      opts.chrom_col,
      opts.pos_col,
      opts.pval_col,
    ]

    print results.data[print_cols].to_string(index=False)
    out_clump = outprefix + ".clump"

    print "\nWriting clumped results to: %s" % out_clump
    results.data.to_csv(out_clump,index=False,sep="\t",na_rep="NA")

    if not opts.skip_gwas:
      print "\nFinding clumped results in LD with GWAS catalog variants..."
      print "\nLD source: %s" % opts.ld_gwas_source
      gwas_ld, gwas_ld_failed_variants = gcat.variants_in_ld_multi(results,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist,opts.threads)

      if gwas_ld is not None:
        # If the user requested other association results columns be merged in with the gwas_ld, add them in.
        if opts.include_cols:
          gwas_ld = merge_include_cols_gwas_hits(gwas_ld,results,opts.include_cols,opts.variant_col)

        out_ld_gwas = outprefix + ".ld-gwas.tab"
        print "\nWriting GWAS catalog variants in LD with clumped variants to: %s" % out_ld_gwas
        gwas_ld.to_csv(out_ld_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits were in LD with clumped variants."

      gwas_near = gcat.variants_nearby(results,opts.gwas_cat_dist)

      if gwas_near is not None:
        # If the user requested other association results columns be merged in with the gwas_hits, add them in.
        if opts.include_cols:
          gwas_near = merge_include_cols_gwas_hits(gwas_near,results,opts.include_cols,opts.variant_col)

        out_near_gwas = outprefix + ".near-gwas.tab"
        print "Writing GWAS catalog variants within %s of a clumped variant to: %s" % (BasePair(opts.gwas_cat_dist).as_kb(),out_near_gwas)
        gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits discovered within %s of any clumped results." % BasePair(opts.gwas_cat_dist).as_kb()

  else:
    if not opts.skip_gwas:
      print "\nFinding results in LD with GWAS catalog variants..."
      print "\nLD source: %s" % opts.ld_gwas_source
      gwas_ld, gwas_ld_failed_variants = gcat.variants_in_ld_multi(results,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist,opts.threads)

      if gwas_ld is not None:
        # If the user requested other columns be merged in with the gwas_ld, pull 'em out.
        if opts.include_cols:
          gwas_ld = merge_include_cols_gwas_hits(gwas_ld,results,opts.include_cols,opts.variant_col)

        print "Found %i GWAS catalog variants in LD with a clumped variant.." % gwas_ld.shape[0]

        out_ld_gwas = outprefix + ".ld-gwas.tab"
        print "\nWriting results to: %s" % out_ld_gwas
        gwas_ld.to_csv(out_ld_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits were in LD with clumped variants."

      gwas_near = gcat.variants_nearby(results,opts.gwas_cat_dist)

      if gwas_near is not None:
        # If the user requested other association results columns be merged in with the gwas_hits, add them in.
        if opts.include_cols:
          gwas_near = merge_include_cols_gwas_hits(gwas_near,results,opts.include_cols,opts.variant_col)

        out_near_gwas = outprefix + ".near-gwas.tab"
        print "Writing GWAS catalog variants within %s of a clumped variant to: %s" % (BasePair(opts.gwas_cat_dist).as_kb(),out_near_gwas)
        gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA")
      else:
        print "\nNo GWAS hits discovered within %s of any clumped results." % BasePair(opts.gwas_cat_dist).as_kb()
    else:
      print "Clumping was not performed, and --skip-gwas was enabled... so there's nothing to do?"

def proc_multi(trait,opts):
  log_obj = None
  try:
    # Setup log file.
    log_file = opts.out + ".%s" % trait + ".log"
    log_obj = codecs.open(log_file,'w','utf-8')

    sys.stdout = StreamTee(sys.stdout,log_obj)
    sys.stderr = StreamTee(sys.stderr,log_obj)

    # Run process for this trait.
    df = multiassoc_epacts_load(opts.assoc,trait,opts.clump_p,opts.rsq_col,opts.rsq_filter,opts.filter)
    out = opts.out + ".%s" % trait
    run_process(df,trait,out,opts)

  except:
    print >> sys.stderr, os.linesep + traceback.format_exc()

  finally:
    #sys.stdout.flush()
    #sys.stderr.flush()

    if log_obj is not None:
      log_obj.close()

def main(arg_string=None):
  print_program_header()
  print ""
  (opts,args) = get_settings(arg_string)

  # If we're running single threaded, everything will go to the same log file.
  # Otherwise, each thread will create its own log.
  if opts.threads == 1:
    # Setup log file.
    log_file = opts.out + ".log"
    log_obj = codecs.open(log_file,'w','utf-8')

    @atexit.register
    def close_log():
      try:
        log_obj.close()
      except:
        pass

    sys.stdout = StreamTee(sys.stdout,log_obj)
    sys.stderr = StreamTee(sys.stderr,log_obj)

  # Loop over traits if multi assoc file, otherwise just load it and do the single trait.
  if opts.multi_assoc:
    if opts.trait is not None:
      out = opts.out + ".%s" % opts.trait
      df = multiassoc_epacts_load(opts.assoc,opts.trait,opts.clump_p,opts.rsq_col,opts.rsq_filter,opts.filter)
      run_process(df,opts.trait,out,opts)

    else:
      if opts.threads == 1:
        print "Running in --multi-assoc mode, loading each trait from assoc file.."
        for trait, df in multiassoc_epacts_iter(opts.assoc,opts.trait,opts.clump_p,opts.rsq_col,opts.rsq_filter,opts.filter):
          out = opts.out + ".%s" % trait
          run_process(df,trait,out,opts)

      else:
        print "Starting threaded.. %i threads allowed simultaneously" % opts.threads

        pool = Pool(opts.threads)
        traits = multiassoc_epacts_get_traits(opts.assoc)

        # Set the thread count back to 1. We don't want the other potentially threaded components
        # to execute threaded (like GWAS catalog lookups.) For some reason you can't have multiprocessing
        # processes executing within multiprocesses.
        opts.threads = 1

        print "Found %i traits in multiassoc file.." % len(traits)

        for trait in traits:
          pool.apply_async(proc_multi,(trait,opts))

        pool.close()
        pool.join()

  else:
    if opts.trait is None and opts.trait_col is None:
      opts.trait = "NA"

    print "Loading results file: %s" % opts.assoc
    run_process(opts.assoc,opts.trait,opts.out,opts)

  return opts, args
      
if __name__ == "__main__":
  main()
