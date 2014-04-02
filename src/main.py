#!/usr/bin/env python
import copy
import os.path
import atexit
import __builtin__
import codecs
import traceback
from optparse import *
from termcolor import *
from glob import glob
from utils import *
from VCFastFinder import *
from PyLDFinder import *
from LDClumper import *
from AssocResults import *
from GWASCatalog import *
from itertools import *
from VerboseParser import *
from multiprocessing import Pool, cpu_count

SWISS_CONF = "conf/swiss.conf";
__builtin__.SWISS_DEBUG = False;

if SWISS_DEBUG:
  pd.set_option('chained_assignment','warn');
else:
  pd.set_option('chained_assignment',None);

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
  parser.add_option("--multi-assoc",help="Designate that the results file is in EPACTS multi-assoc format.",action="store_true",default=False);
  parser.add_option("--trait",help="Description of phenotype for association results file. E.g. 'HDL' or 'T2D'");
  parser.add_option("--delim",help="Association results delimiter.",default="\t");
  parser.add_option("--build",help="Genome build your association results are anchored to.",default="hg19");
  parser.add_option("--snp-col",help="SNP column name in results file.",default="MarkerName");
  parser.add_option("--pval-col",help="P-value column name in results file.",default="P-value");
  parser.add_option("--chrom-col",help="Chromosome column name in results file.",default="CHR");
  parser.add_option("--pos-col",help="Position column name in results file.",default="POS");
  parser.add_option("--rsq-col",help="Imputation quality column name.",default="RSQ");

  # Association result filtering results. 
  parser.add_option("--rsq-filter",help="Remove variants below this imputation quality.",default=None);
  parser.add_option("--filter",help="Give a general filter string to filter variants.",default=None);

  # Output options. 
  parser.add_option("--out",help="Prefix for output files.",default="swiss_output");

  # LD clumping options. 
  parser.add_option("--ld-clump",help="Clump association results by LD.",action="store_true",default=False);
  parser.add_option("--clump-p",help="P-value threshold for LD and distance based clumping.",default=5e-08);
  parser.add_option("--clump-ld-thresh",help="LD threshold for clumping.",default=0.2);
  parser.add_option("--clump-ld-dist",help="Distance from each significant result to calculate LD.",default=1E6);

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

  # Misc options
  parser.add_option("-T","--threads",default=1,type="int",help="Number of parallel jobs to run. Only works with --multi-assoc currently.");

  (opts,args) = parser.parse_args();

  conf = getConf();

  if opts.threads < 1:
    opts.threads = 1;

  if opts.threads > cpu_count():
    print >> sys.stderr, "Warning: you set threads to %i, but the CPU count for this machine is %i..." % (opts.threads,cpu_count());

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

  # File delimiter
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

  # If multi-assoc is specified, the column names are already known.
  if opts.multi_assoc:
    opts.snp_col = "MARKER_ID";
    opts.pval_col = "PVALUE";
    opts.chrom_col = "#CHROM";
    opts.pos_col = "BEG";

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

def get_header(infile,sep="\t"):
  if infile.endswith(".gz"):
    f = gzip.open(infile);
  else:
    f = open(infile);

  with f:
    h = f.readline().split(sep);
    h[-1] = h[-1].rstrip();

    return h;

# Helper function to iterate over separate traits from an EPACTS multi assoc file.
def multiassoc_epacts_iter(result_file):
  header = get_header(result_file);

  #CHROM    BEG    END       MARKER_ID    NS        AC  CALLRATE   GENOCNT      MAF  DHA.P   DHA.B  EstC.P  EstC.B  FAw3.P  FAw3.B  FAw3toFA.P  FAw3toFA.B  FAw6.P  FAw6.B  FAw6toFA.P

  trait_ps = filter(lambda x: x.endswith(".P"),header);
  trait_betas = filter(lambda x: x.endswith(".B"),header);
  trait_cols = trait_ps + trait_betas;
  intro_cols = filter(lambda x: x not in trait_cols,header);

  # Load first columns, since we'll always use them, along with the first trait column.
  # base_trait = first trait to load
  base_trait_p = trait_ps[0];
  base_trait_b = trait_betas[0];
  base_trait = base_trait_p.replace(".P","");
  base_cols = intro_cols + [base_trait_p,base_trait_b];

  print "\nLoading trait %s from association results file: %s" % (base_trait,result_file);

  base_df = pd.read_table(result_file,compression = "gzip" if result_file.endswith(".gz") else None,na_values=["NA","None","."],usecols = base_cols);
  base_df.rename(columns = {base_trait_p : "PVALUE",base_trait_b : "BETA"},inplace=True);

  yield (base_trait,base_df);

  for p in trait_ps[1:]:
    trait = p.replace(".P","");
    bcol = trait + ".B";

    print "\nLoading trait %s from association results file: %s" % (trait,result_file);

    df = pd.read_table(result_file,compression = "gzip" if result_file.endswith(".gz") else None,na_values=["NA","None","."],usecols = [p,bcol]);
    df.rename(columns = {p : "PVALUE",bcol : "BETA"},inplace=True);

    base_df["PVALUE"] = df["PVALUE"];
    base_df["BETA"] = df["BETA"];

    yield (trait,base_df);

def multiassoc_epacts_load(result_file,trait):
  print "\nLoading trait %s from association results file: %s" % (trait,result_file);

  header = get_header(result_file);

  #CHROM    BEG    END       MARKER_ID    NS        AC  CALLRATE   GENOCNT      MAF  DHA.P   DHA.B  EstC.P  EstC.B  FAw3.P  FAw3.B  FAw3toFA.P  FAw3toFA.B  FAw6.P  FAw6.B  FAw6toFA.P

  # Check to make sure the trait requested is in the header.
  trait_ps = filter(lambda x: x.endswith(".P"),header);
  trait_betas = filter(lambda x: x.endswith(".B"),header);
  trait_cols = trait_ps + trait_betas;
  intro_cols = filter(lambda x: x not in trait_cols,header);

  trait_names = map(lambda x: x.replace(".P",""),trait_ps);
  if trait not in trait_names:
    raise IOError, "Requested trait %s is not present in %s" % (trait,result_file);

  this_trait_cols = [trait + ".P",trait + ".B"];
  df = pd.read_table(result_file,compression = "gzip" if result_file.endswith(".gz") else None,na_values=["NA","None","."],usecols = intro_cols + this_trait_cols);
  df.rename(columns = {
    trait + ".P" : "PVALUE",
    trait + ".B" : "BETA"
  },inplace=True);

  return df;

def multiassoc_epacts_get_traits(result_file):
  header = get_header(result_file);

  # Check to make sure the trait requested is in the header.
  trait_ps = filter(lambda x: x.endswith(".P"),header);
  traits = map(lambda x: x.replace(".P",""),trait_ps);

  return traits;

def run_process(assoc,trait,outprefix,opts):
  if isinstance(assoc,str):
    results = AssocResults(assoc,trait);
  else:
    results = AssocResults(trait=trait,df=assoc);

  results.marker_col = opts.snp_col;
  results.chrom_col = opts.chrom_col;
  results.pos_col = opts.pos_col;
  results.pval_col = opts.pval_col;
  results.rsq_col = opts.rsq_col;
  results.load(sep=opts.delim);

  # Filter results on imputation quality, if requested.
  if opts.rsq_filter is not None:
    results.filter_imp_quality(opts.rsq_filter);

  # If the user specified an arbitrary filter, run it too.
  if opts.filter is not None:
    results.do_filter(opts.filter);

#   # LD finder for clumping
#   vset = VCFastSettings(opts.ld_clump_source_file,opts.vcfast_path);
# #  cache_clump = LDRegionCache(vset.createLDCacheKey(),opts.cache + "_clump.db");
# #  atexit.register(cache_clump.close);
#   finder_clumping = VCFastFinder(vset,verbose=False,cache=None);
#
#   # LD finder for GWAS catalog lookups
#   vset_gwas = VCFastSettings(opts.ld_gwas_source_file,opts.vcfast_path);
# #  cache_gwas = LDRegionCache(vset.createLDCacheKey(),opts.cache + "_gwas.db");
# #  atexit.register(cache_gwas.close);
#   finder_gwas = VCFastFinder(vset_gwas,verbose=False,cache=None);

  # LD finder for clumping
  vset = PyLDSettings(opts.ld_clump_source_file,opts.tabix_path);
  finder_clumping = PyLDFinder(vset,verbose=False,cache=None);

  # LD finder for GWAS catalog lookups
  vset_gwas = PyLDSettings(opts.ld_gwas_source_file,opts.tabix_path);
  finder_gwas = PyLDFinder(vset_gwas,verbose=False,cache=None);

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
    out_clump = outprefix + ".clump";

    print "\nWriting clumped results to: %s" % out_clump;
    results_clumped.data.to_csv(out_clump,index=False,sep="\t",na_rep="NA");

    print "\nFinding clumped results in LD with GWAS catalog variants...";
    print "\nLD source: %s" % opts.ld_gwas_source;
    gwas_hits, gwas_ld_failed_variants = gcat.variants_in_ld(results_clumped,finder_gwas,opts.gwas_cat_ld,opts.gwas_cat_dist);

    # If the user requested other columns be merged in with the gwas_hits, pull 'em out.
    if opts.include_cols:
      include_cols = [i.strip() for i in opts.include_cols.split(",")];
      include_cols = filter(lambda x: x in results_clumped.data.columns,include_cols);

      if len(include_cols) == 0:
        print "Warning: user specified --include-cols, but none of them existed in the association results!";
      else:
        assoc_incl_cols = results_clumped.data[[opts.snp_col] + include_cols];
        assoc_incl_cols.rename(
          columns = dict(zip(include_cols,map(lambda x: "ASSOC_" + x,include_cols))),
          inplace = True
        );
        gwas_hits = pd.merge(gwas_hits,assoc_incl_cols,left_on="ASSOC_MARKER",right_on=opts.snp_col);
        del gwas_hits[opts.snp_col];

    out_ld_gwas = outprefix + ".ld-gwas.tab";
    print "\nWriting GWAS catalog variants in LD with clumped variants to: %s" % out_ld_gwas;
    gwas_hits.to_csv(out_ld_gwas,index=False,sep="\t",na_rep="NA");

    gwas_near = gcat.variants_nearby(results_clumped,opts.gwas_cat_dist);

    if gwas_near is not None:
      out_near_gwas = outprefix + ".near-gwas.tab";
      print "Writing GWAS catalog variants within %s of a clumped variant to: %s" % (BasePair(opts.gwas_cat_dist).as_kb(),out_near_gwas);
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA");

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
    out_clump = outprefix + ".clump.tab";

    print "\nWriting clumped results to: %s" % out_clump;
    results.data.to_csv(out_clump,index=False,sep="\t",na_rep="NA");

    gwas_near = gcat.variants_nearby(results,opts.gwas_cat_dist);
    if gwas_near.shape[0] > 0:
      print "\nFor those variants for which LD buddies could not be computed, there were %i variants within %s of a GWAS hit." % (gwas_near.shape[0],BasePair(opts.gwas_cat_dist).as_kb());
      out_near_gwas = outprefix + ".near-gwas.tab";

      print "These variants were written to: %s" % out_near_gwas;
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA");
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

    print "Found %i GWAS catalog variants in LD with a clumped variant.." % gwas_hits.shape[0];

    out_ld_gwas = outprefix + ".ld-gwas.tab";
    print "\nWriting results to: %s" % out_ld_gwas;
    gwas_hits.to_csv(out_ld_gwas,index=False,sep="\t",na_rep="NA");

    gwas_near = gcat.variants_nearby(results,opts.gwas_cat_dist);

    if gwas_near is not None:
      print "\nThere were %i variants within %s of a GWAS hit." % (gwas_near.shape[0],BasePair(opts.gwas_cat_dist).as_kb());
      out_near_gwas = outprefix + ".near-gwas.tab";
      print "These variants were written to: %s" % out_near_gwas;
      gwas_near.to_csv(out_near_gwas,index=False,sep="\t",na_rep="NA");

def proc_multi(trait,opts):
  log_obj = None;
  try:
    # Setup log file.
    log_file = opts.out + ".%s" % trait + ".log";
    log_obj = codecs.open(log_file,'w','utf-8');

    sys.stdout = StreamTee(sys.stdout,log_obj);
    sys.stderr = StreamTee(sys.stderr,log_obj);

    # Run process for this trait.
    df = multiassoc_epacts_load(opts.assoc,trait);
    out = opts.out + ".%s" % trait;
    run_process(df,trait,out,opts);

  except:
    print >> sys.stderr, traceback.print_exc();

  finally:
    if log_obj is not None:
      log_obj.close();

def main():
  (opts,args) = get_settings();

  # If we're running single threaded, everything will go to the same log file.
  # Otherwise, each thread will create its own log.
  if opts.threads == 1:
    # Setup log file.
    log_file = opts.out + ".log";
    log_obj = codecs.open(log_file,'w','utf-8');

    @atexit.register
    def close_log():
      try:
        log_obj.close();
      except:
        pass

    sys.stdout = StreamTee(sys.stdout,log_obj);
    sys.stderr = StreamTee(sys.stderr,log_obj);

  # Loop over traits if multi assoc file, otherwise just load it and do the single trait.
  if opts.multi_assoc:
    if opts.threads == 1:
      print "Running in --multi-assoc mode, loading each trait from assoc file..";
      for trait, df in multiassoc_epacts_iter(opts.assoc):
        out = opts.out + ".%s" % trait;
        run_process(df,trait,out,opts);

    else:
      print "Starting threaded.. %i threads allowed simultaneously" % opts.threads;

      pool = Pool(opts.threads);
      traits = multiassoc_epacts_get_traits(opts.assoc);

      print "Found %i traits in multiassoc file.." % len(traits);

      for trait in traits:
        pool.apply_async(proc_multi,(trait,opts));

      pool.close();
      pool.join();

  else:
    print "Loading trait %s from results file: %s" % (opts.trait,opts.assoc);
    run_process(opts.assoc,opts.trait,opts.out,opts);
      
if __name__ == "__main__":
  main();
