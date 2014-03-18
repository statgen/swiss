#!/usr/bin/env python
import copy
import os.path
import atexit
import __builtin__
import codecs
from optparse import *
from termcolor import *
from glob import glob
from utils import *
from VCFastFinder import *
from LDClumper import *
from AssocResults import *
from GWASCatalog import *
from itertools import *
from VerboseParser import *

SWISS_CONF = "conf/swiss.conf";
__builtin__.SWISS_DEBUG = True;

class BasePair:
  def __init__(self,bp):
    self.bp = bp;

  def as_kb(self):
    kb = float(self.bp) / 1000.0;
    if int(kb) == kb:
      kb = int(kb);

    return "%ikb" % kb;

class Conf(object):
  def __init__(self,conf_file):
    self._load(conf_file);

  def _load(self,file):
    conf_dict = {};
    execfile(file,conf_dict);

    for k,v in conf_dict.iteritems():
      exec "self.%s = v" % str(k);

def getConf(conf_file=SWISS_CONF):
  conf_file = find_relative(conf_file);
  conf = Conf(conf_file);
  return conf;

def get_settings():
  usage = "swiss [options]";
  parser = VerboseParser(usage=usage);

  # Association result input options. 
  parser.add_option("--assoc",help="[Required] Association results file.");
  parser.add_option("--trait",help="Description of phenotype for association results file. E.g. 'HDL' or 'T2D'");
  parser.add_option("--delim",help="Association results delimiter.",default="\t");
  parser.add_option("--build",help="Genome build your association results are anchored to.",default="hg19");
  parser.add_option("--snp-col",help="SNP column name in results file.",default="MarkerName");
  parser.add_option("--pval-col",help="P-value column name in results file.",default="P-value");
  parser.add_option("--chrom-col",help="Chromosome column name in results file.",default="CHR");
  parser.add_option("--pos-col",help="Position column name in results file.",default="POS");

  # Output options. 
  parser.add_option("--out",help="Prefix for output files.",default="swiss_output");

  # LD clumping options. 
  parser.add_option("--ld-clump",help="Clump association results by LD.",action="store_true",default=False);
  parser.add_option("--clump-p",help="P-value threshold for LD and distance based clumping.",default=5e-08);
  parser.add_option("--clump-ld-thresh",help="LD threshold for clumping.",default=0.2);
  parser.add_option("--clump-ld-dist",help="Distance from each significant result to calculate LD.",default=2E6);

  # Distance clumping options. 
  parser.add_option("--dist-clump",help="Clump association results by distance.",action="store_true",default=False);
  parser.add_option("--clump-dist",help="Distance threshold to use for clumping based on distance.",default=2.5e5);

  # LD source (GoT2D, 1000G, etc.) 
  parser.add_option("--ld-clump-source",help="Name of pre-configured LD source, or a VCF file from which to compute LD.",default="GOT2D_2011-11");
  parser.add_option("--list-ld-sources",help="Print a list of available LD sources for each genome build.",default=False,action="store_true");

  # GWAS catalog
  parser.add_option("--gwas-cat",help="GWAS catalog to use.",default="fusion");
  parser.add_option("--ld-gwas-source",help="Name of pre-configured LD source or VCF file to use when calculating LD with GWAS variants.",default="GOT2D_2011-11");
  parser.add_option("--list-gwas-cats",action="store_true",default=False,help="Give a listing of all valid GWAS catalogs and their descriptions.");
  parser.add_option("--gwas-cat-p",help="P-value threshold for GWAS catalog variants.",default=5e-08);
  parser.add_option("--gwas-cat-ld",help="LD threshold for considering a GWAS catalog variant in LD.",default=0.1);
  parser.add_option("--gwas-cat-dist",help="Distance threshold for considering a GWAS catalog variant 'nearby'.",default=2.5e5);
  parser.add_option("--include-cols",help="List of columns to merge in from association results (grouped by variant.)",default=None);

  # LD cache
  parser.add_option("--cache",help="Prefix for LD cache.",default="ld_cache");

  (opts,args) = parser.parse_args();

  conf = getConf();

  if opts.list_gwas_cats:
    print_gwas();
    sys.exit(0);

  if opts.list_ld_sources:
    print_ld_sources();
    sys.exit(0);

  if opts.assoc is None:
    parser.print_help();
    print "";
    print >> sys.stderr, "Must specify association results file: --assoc";
    sys.exit(1);

  if not os.path.isfile(opts.assoc):
    parser.print_help();
    print "";
    print >> sys.stderr, "Association results file does not exist: %s" % opts.assoc;
    sys.exit(1);

  if opts.ld_clump:
    opts.clump_ld_thresh = float(opts.clump_ld_thresh);
    if opts.clump_ld_thresh < 0 or opts.clump_ld_thresh > 1:
      error("LD threshold for clumping must be >= 0 or <= 1.");

  opts.clump_p = float(opts.clump_p);
  if opts.clump_p < 0 or opts.clump_p > 1:
    error("P-value threshold must be >= 0 or <= 1.");

  opts.clump_ld_dist = int(float(opts.clump_ld_dist));

  if not os.path.isfile(opts.gwas_cat):
    opts.gwas_cat_file = find_relative(conf.GWAS_CATALOGS[opts.build][opts.gwas_cat]);

    if opts.gwas_cat_file == None:
      error("Could not locate GWAS catalog file for conf entry '%s - %s'!" % (opts.build,opts.gwas_cat));
  else:
    opts.gwas_cat_file = opts.gwas_cat;

  # If one source is specified, but not the other, assume the user meant to use this for both. 
  # ld clump source XOR ld gwas source
  if (opts.ld_clump_source is None) != (opts.ld_gwas_source is None):
    if opts.ld_clump_source is None:
      opts.ld_clump_source = opts.ld_gwas_source;
    elif opts.ld_gwas_source is None:
      opts.ld_gwas_source = opts.ld_clump_source;
    
  # LD clumping source
  if not os.path.isfile(opts.ld_clump_source):
    opts.ld_clump_source_file = find_relative(conf.LD_SOURCES[opts.build][opts.ld_clump_source]);

    if opts.ld_clump_source_file == None:
      error("Could not locate VCF file for conf entry '%s - %s'!" % (opts.build,opts.ld_clump_source));
  else:
    opts.ld_clump_source_file = opts.ld_clump_source;
  
  # GWAS LD lookup source
  if not os.path.isfile(opts.ld_gwas_source):
    opts.ld_gwas_source_file = find_relative(conf.LD_SOURCES[opts.build][opts.ld_gwas_source]);

    if opts.ld_gwas_source_file == None:
      error("Could not locate VCF file for conf entry '%s - %s'!" % (opts.build,opts.ld_gwas_source));
  else:
    opts.ld_gwas_source_file = opts.ld_gwas_source;

  if not os.path.isfile(opts.assoc):
    error("Could not locate association results file (or insufficient permissions): %s" % opts.assoc);

  if opts.delim in ("tab","\t","\\t"):
    opts.delim = "\t";
  elif opts.delim in ("comma",","):
    opts.delim = ",";
  elif opts.delim in ("space"," "):
    opts.delim = " ";
  elif opts.delim == "whitespace":
    opts.delim = None;

  opts.vcfast_path = find_systematic(conf.VCFAST_PATH);
  opts.tabix_path = find_systematic(conf.TABIX_PATH);

  out_exists = glob(os.path.join(opts.out,"*"));
  if len(out_exists) > 0:
    error("Output files already exist with this prefix: %s" % opts.out);

  return (opts,args);

def print_gwas():
  conf = getConf();

  print "%-10s %-40s" % ("Build","Catalog");
  print "%-10s %-40s" % ("-----","-------");
  for build, cats in conf.GWAS_CATALOGS.iteritems():
    print "%-10s %-40s" % (build,",".join(cats.keys()));

def print_ld_sources():
  conf = getConf();

  print "%-10s %-40s" % ("Build","LD Sources");
  print "%-10s %-40s" % ("-----","----------");
  for build, ld_source in conf.LD_SOURCES.iteritems():
    ld_source_names = sorted(ld_source);
    print "%-10s %-40s" % (build,", ".join(ld_source_names));

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
    
def main():
  (opts,args) = get_settings();

  # Setup log file. 
  log_file = opts.out + ".log";
  log_obj = codecs.open(log_file,'w','utf-8');
  
  @atexit.register
  def close_log():
    try:
      log_obj.close();
    except:
      pass

  # TODO: this screws up ipython somehow
  # might be able to set it back when exception is thrown or something?
  # better long term solution might be to use logging...
  sys.stdout = StreamTee(sys.stdout,log_obj);
  sys.stderr = StreamTee(sys.stderr,log_obj);

  results = AssocResults(opts.assoc,opts.trait);
  results.marker_col = opts.snp_col;
  results.chrom_col = opts.chrom_col;
  results.pos_col = opts.pos_col;
  results.pval_col = opts.pval_col;
  results.load(sep=opts.delim);

  # LD finder for clumping 
  vset = VCFastSettings(opts.ld_clump_source_file,opts.vcfast_path);
#  cache_clump = LDRegionCache(vset.createLDCacheKey(),opts.cache + "_clump.db");
#  atexit.register(cache_clump.close);
  finder_clumping = VCFastFinder(vset,verbose=False,cache=None);
  
  # LD finder for GWAS catalog lookups
  vset_gwas = VCFastSettings(opts.ld_gwas_source_file,opts.vcfast_path);
#  cache_gwas = LDRegionCache(vset.createLDCacheKey(),opts.cache + "_gwas.db");
#  atexit.register(cache_gwas.close);
  finder_gwas = VCFastFinder(vset_gwas,verbose=False,cache=None);
  
  # GWAS catalog. 
  gcat = GWASCatalog(opts.gwas_cat_file);

  print "\nIdentifying GWAS catalog variants that do not overlap with your --ld-gwas-source: %s" % opts.ld_gwas_source;
  missing_vcf = gcat.variants_missing_vcf(opts.ld_gwas_source_file);
  missing_vcf.sort('PHENO',inplace=True);
  missing_vcf = sort_genome(missing_vcf,'CHR','POS');
  print colored('Warning: ','yellow') + "the following variants in the GWAS catalog are not present in your VCF file: ";
  print missing_vcf["SNP CHR POS PHENO Group".split()].to_string(index=False);
  
  print "\nLoaded %i variants from association results.." % results.data.shape[0];

  if opts.ld_clump:
    print "\nLD clumping results..";
    print "\nLD source: %s" % opts.ld_clump_source;
    clumper = LDClumper(results,finder_clumping);
    clump_results = clumper.ld_clump(opts.clump_p,opts.clump_ld_thresh,opts.clump_ld_dist);

    if clump_results is not None:
      results_clumped, failed_ld_clump_variants = clump_results;
    else:
      print "\nNo significant variants available for LD clumping, skipping..";
      return;

    print "\nResults after clumping: ";
    print_cols = [
      opts.snp_col,
      opts.chrom_col,
      opts.pos_col,
      opts.pval_col,
      'ld_with',
      'failed_clump'
    ];
    print results_clumped.data[print_cols].to_string(index=False);
    out_clump = opts.out + ".clump";
    print "\nWriting clumped results to: %s" % out_clump;
    results_clumped.data.to_csv(out_clump,index=False,sep="\t");

    print "\nFinding clumped results in LD with GWAS catalog variants...";
    print "\nLD source: %s" % opts.ld_gwas_source;
    gwas_hits, gwas_ld_failed_variants = gcat.variants_in_ld(results_clumped,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist);

    # If the user requested other columns be merged in with the gwas_hits, pull 'em out. 
    if opts.include_cols:
      include_cols = [i.strip() for i in opts.include_cols.split(",")];
      include_cols = filter(lambda x: x in results_clumped_nofail.data.columns,include_cols);

      if len(include_cols) == 0:
        print "Warning: user specified --include-cols, but none of them existed in the association results!";
      else:
        assoc_incl_cols = results_clumped_nofail.data[[opts.snp_col] + include_cols];
        assoc_incl_cols.rename(
          columns = dict(zip(include_cols,map(lambda x: "ASSOC_" + x,include_cols))),
          inplace = True
        );
        gwas_hits = pd.merge(gwas_hits,assoc_incl_cols,left_on="ASSOC_MARKER",right_on=opts.snp_col);
        del gwas_hits[opts.snp_col];

    out_ld_gwas = opts.out + ".ld-gwas.tab";
    print "\nWriting results to: %s" % out_ld_gwas;
    gwas_hits.to_csv(out_ld_gwas,index=False,sep="\t");

    results_clumped_ldfail = copy.deepcopy(results_clumped);
    results_clumped_ldfail.keep_variants(gwas_ld_failed_variants);

    gwas_near = gcat.variants_nearby(results_clumped_ldfail,opts.gwas_cat_dist);
    if gwas_near.shape[0] > 0:
      print "\nFor those variants for which LD buddies could not be computed, there were %i variants within %s of a GWAS hit." % (gwas_near.shape[0],BasePair(opts.gwas_cat_dist).as_kb());
      out_near_gwas = opts.out + ".near-gwas.tab";
      print "These variants were written to: %s" % out_near_gwas;
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t");
    else:
      print "\nFor those variants that did not exist in the VCF file (for computing LD buddies for GWAS hits), there were no GWAS hits within %s." % BasePair(opts.gwas_cat_dist).as_kb();

  elif opts.dist_clump:
    results.dist_clump(opts.clump_p,opts.clump_dist);

    print "\nResults after clumping: ";
    print_cols = [
      opts.snp_col,
      opts.chrom_col,
      opts.pos_col,
      opts.pval_col,
    ];
    print results.data[print_cols].to_string(index=False);
    out_clump = opts.out + ".clump.tab";
    print "\nWriting clumped results to: %s" % out_clump;
    results.data.to_csv(out_clump,index=False,sep="\t");

    gwas_near = gcat.variants_nearby(results,opts.gwas_cat_dist);
    if gwas_near.shape[0] > 0:
#      print "\nGWAS hits within %s of clumped results: " % BasePair(opts.gwas_cat_dist).as_kb();
#      print_cols = "ASSOC_MARKER ASSOC_CHRPOS ASSOC_TRAIT GWAS_SNP GWAS_CHRPOS ASSOC_GWAS_DIST BUILD PHENO GENE_LABEL POPULATION CITATION".split();
#      print gwas_near[print_cols].to_string(index=False);

      print "\nFor those variants for which LD buddies could not be computed, there were %i variants within %s of a GWAS hit." % (gwas_near.shape[0],BasePair(opts.gwas_cat_dist).as_kb());
      out_near_gwas = opts.out + ".near-gwas.tab";
      print "These variants were written to: %s" % out_near_gwas;
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t");
    else:
      print "\nNo GWAS hits discovered within %s of any clumped results." % BasePair(opts.gwas_cat_dist).as_kb();

  else:
    print "\nFinding results in LD with GWAS catalog variants...";
    print "\nLD source: %s" % opts.ld_gwas_source;
    gwas_hits, gwas_ld_failed_variants = gcat.variants_in_ld(results,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist);

    # If the user requested other columns be merged in with the gwas_hits, pull 'em out. 
    if opts.include_cols:
      include_cols = [i.strip() for i in opts.include_cols.split(",")];
      include_cols = filter(lambda x: x in results.data.columns,include_cols);

      if len(include_cols) == 0:
        print "Warning: user specified --include-cols, but none of them existed in the association results!";
      else:
        assoc_incl_cols = results.data[[opts.snp_col] + include_cols];
        assoc_incl_cols.rename(
          columns = dict(zip(include_cols,map(lambda x: "ASSOC_" + x,include_cols))),
          inplace = True
        );
        gwas_hits = pd.merge(gwas_hits,assoc_incl_cols,left_on="ASSOC_MARKER",right_on=opts.snp_col);
        del gwas_hits[opts.snp_col];

    out_ld_gwas = opts.out + ".ld-gwas.tab";
    print "\nWriting results to: %s" % out_ld_gwas;
    gwas_hits.to_csv(out_ld_gwas,index=False,sep="\t");

    results_clumped_ldfail = copy.deepcopy(results_clumped);
    results_clumped_ldfail.keep_variants(gwas_ld_failed_variants);

    gwas_near = gcat.variants_nearby(results_clumped_ldfail,opts.gwas_cat_dist);
    if gwas_near.shape[0] > 0:
      out_near_gwas = opts.out + ".near-gwas.tab";
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t");
    else:
      print "\nFor those variants that did not exist in the VCF file (for computing LD buddies for GWAS hits), there were no GWAS hits within %s." % BasePair(opts.gwas_cat_dist).as_kb();
      
if __name__ == "__main__":
  main();
