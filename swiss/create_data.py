#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2016 Ryan Welch, The University of Michigan
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

from __future__ import print_function
import os, sys, re, gzip, tempfile, io
import math, time, traceback, sqlite3
import os.path as path
from tqdm import tqdm
from optparse import OptionParser
from collections import namedtuple
from six.moves.urllib.request import urlretrieve, urlcleanup
from six import itervalues, text_type
from itertools import chain
from toolz.itertoolz import partition_all

# Constants.
SQLITE_SNP_POS = "dbsnp"
SQLITE_TRANS = "trans"
SNP_HISTORY_URL = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/SNPHistory.bcp.gz"
RS_MERGE_ARCH_URL = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606/database/organism_data/RsMergeArch.bcp.gz"
NCBI_VCF_TEMPLATE_URL = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_{0}_{1}/VCF/00-All.vcf.gz"
GWAS_CAT_EBI_URL = "http://www.ebi.ac.uk/gwas/api/search/downloads/full"
GWAS_LOG_SIG = -math.log10(5e-08)

def tqdm_hook(t):
  """
  Args:
    t: tqdm object

  Wraps tqdm instance. Don't forget to close() or __exit__()
  the tqdm instance once you're done with it (easiest using `with` syntax).
  """

  last_b = [0]

  def inner(b=1,bsize=1,tsize=None):
    """
    Args:
      b : int, optional
          Number of blocks just transferred [default: 1].
      bsize : int, optional
          Size of each block (in tqdm units) [default: 1].
      tsize : int, optional
          Total size (in tqdm units). If [default: None] remains unchanged.
    """

    if tsize is not None:
      t.total = tsize

    t.update((b - last_b[0]) * bsize)
    last_b[0] = b

  return inner

def remove_file(filepath):
  try:
    os.remove(filepath)
  except:
    pass

def download_ebi_catalog(url,outpath):
  with tqdm(unit='B',unit_scale=True,miniters=1,desc=url.split('/')[-1]) as t:
    urlcleanup()
    urlretrieve(url,filename=outpath,reporthook=tqdm_hook(t),data=None)

  return outpath

def download_merge_arch(url=RS_MERGE_ARCH_URL,outpath="RsMergeArch.bcp.gz"):
  with tqdm(unit='B',unit_scale=True,miniters=1,desc=url.split('/')[-1]) as t:
    urlcleanup()
    urlretrieve(url,filename=outpath,reporthook=tqdm_hook(t),data=None)

  return outpath

def download_snp_history(url=SNP_HISTORY_URL,outpath="SNPHistory.bcp.gz"):
  with tqdm(unit='B',unit_scale=True,miniters=1,desc=url.split('/')[-1]) as t:
    urlcleanup()
    urlretrieve(url,filename=outpath,reporthook=tqdm_hook(t),data=None)

  return outpath

def download_dbsnp_vcf(dbsnp_build=None,genome_build=None,url=None,outpath=None):
  """
  Download the NCBI dbSNP VCF for a given human genome build and dbSNP build

  Args:
    dbsnp_build: b147
    genome_build: GRCh37p13
    url: Direct URL to file, e.g. ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b147_GRCh37p13/VCF/00-All.vcf.gz
    outpath: Constructed from genome_build and dbsnp_build. If not given, a random filename will be generated.

  Returns:
    Name of file into which we saved the data (will be constructed from builds, or random name)
  """

  if url is None:
    if not genome_build.startswith("GRC"):
      raise ValueError("Genome build should begin with GRC")

    if not dbsnp_build.startswith("b"):
      raise ValueError("dbSNP build should look like b147, b148, etc.")

    url = NCBI_VCF_TEMPLATE_URL.format(dbsnp_build,genome_build)

  if outpath is None:
    if genome_build is None or dbsnp_build is None:
      outpath = "dbsnp.vcf.gz"
    else:
      outpath = "human_9606_{}_{}_All.vcf.gz".format(dbsnp_build,genome_build)

  with tqdm(unit='B',unit_scale=True,miniters=1,desc=url.split('/')[-1]) as t:
    urlcleanup()
    urlretrieve(url,filename=outpath,reporthook=tqdm_hook(t),data=None)

  return outpath

def parse_year(x):
  try:
    tup = time.strptime(x,"%m/%d/%Y")
    fix = time.strftime("%Y",tup)
  except:
    return ""
  else:
    return fix

def parse_strongest_snp_risk(field):
  match = re.search("rs(\d+)-(\w)",field)
  if match is not None:
    digits, allele = match.groups()
    snp = "rs" + digits

    return snp, allele
  else:
    return None, None

def parse_or_beta(field):
  field = field.strip()
  if field == "NR":
    return float("nan")
  else:
    try:
      field = float(field)
    except:
      field = float("nan")

    return field

def parse_epacts(v,strict=True):
  """
  Try to parse an EPACTS ID into components.

  Args:
    strict (bool): if true, must match an EPACTS ID exactly (chr:pos_ref/alt)
      If false, then ref/alt can be missing, but at least chr:pos must be specified.
      In this case, ref/alt will be None in the returned tuple

  Returns:
    tuple: (chrom, pos, ref, alt, extra)
  """

  split = v.split("_")

  # Split chrom/pos
  # This is the minimum required information. If even this isn't present,
  # it's a bogus ID.
  try:
    chrom, pos = split[0].split(":")
  except:
    raise ValueError("EPACTS ID had no chrom/pos? " + v)

  # Position should be numeric
  try:
    int(pos)
  except:
    raise ValueError("Couldn't recognize position {} for variant {}".format(pos,v))

  # Try to split alleles if they were given
  try:
    ref, alt = split[1].split("/")
  except:
    if strict:
      raise ValueError("No ref/alt alleles found in EPACTS ID " + v)

    ref = None
    alt = None

  return [chrom,pos,ref,alt,"_".join(split[2:])]

Variant = namedtuple("Variant","rsid chrpos epacts chrom pos ref alt")

class SwissDB:
  def __init__(self,db_file):
    if not os.path.isfile(db_file):
      sys.exit("Error: could not locate SQLite database file: %s. Check conf file setting SQLITE_DB." % db_file)

    db = sqlite3.connect(db_file)

    db.execute("""
      CREATE TEMP VIEW trans_view AS SELECT rs_orig as rsid,chrpos,epacts,chrom,pos,ref,alt FROM {} p INNER JOIN {} t ON (t.rs_current = p.rsid)
    """.format("dbsnp","trans"))

    self.db = db

  def rsid_to_variant(self,variant):
    if not variant.startswith("rs"):
      raise ValueError("Expected variant ID to be an rsID, got {} instead".format(variant))

    cur1 = self.db.execute("SELECT * FROM trans_view WHERE rsid='{}'".format(variant))
    cur2 = self.db.execute("SELECT * FROM dbsnp WHERE rsid='{}'".format(variant))

    vobj = None
    res = 0
    for row in chain(cur1,cur2):
      rsid, chrpos, epacts, chrom, pos, ref, alt = row
      vobj = Variant(rsid,chrpos,epacts,chrom,pos,ref,alt)
      res += 1

    if res > 1:
      print("Warning: Variant {} has more than 1 position in database".format(variant), file=sys.stderr)

    return vobj

  def __del__(self):
    self.db.close()

def set_and_increment(d,k):
  v = d.get(k,0)
  v += 1
  d[k] = v

def append_cap(d,k,v,cap=5):
  n = d.setdefault(k,[])
  if len(n) < cap:
    n.append(v)

class FilterLog:
  def __init__(self,limit=5):
    self.data = dict()

  def add_warning(self,warning,example):
    root = self.data.setdefault("warning",{})
    entry = root.setdefault(warning,{})

    set_and_increment(entry,"count")
    append_cap(entry,"examples",example)

  def add_filter(self,ffilter,example):
    root = self.data.setdefault("filter",{})
    entry = root.setdefault(ffilter,{})

    set_and_increment(entry,"count")
    append_cap(entry,"examples",example)

  def summary(self):
    for dtype, root in self.data.iteritems():
      for istring, info in root.iteritems():
        print("{}: {}".format(dtype,istring))
        print("\t"*1 + ".. occurred {:,} times".format(info["count"]))
        print("\t"*1 + ".. examples: ")
        for example in info["examples"]:
          print("\t"*2 + example)

        print("")

def parse_gwas_catalog(filepath,dbpath,outpath):
  swiss_db = SwissDB(dbpath)
  log = FilterLog()

  with io.open(filepath,"r",encoding="utf-8") as f, io.open(outpath,'w',encoding="utf-8") as out:
    f.readline() # header

    new_header = u"\t".join("VARIANT EPACTS CHRPOS CHR POS REF ALT PHENO GROUP LOG_PVAL CITATION RISK_ALLELE RISK_AL_FREQ GENE_LABEL OR_BETA".split())
    print(new_header,file=out)

    seen_trait_snps = set()
    for line in tqdm(f,unit="line"):
      if line.strip() == "":
        continue

      e = line.split("\t")

      trait, rsids, log_pval, or_beta = (e[i] for i in (7,21,28,30))
      author, date, journal = (e[i] for i in (2,3,4))

      # Fix the OR/beta column to be either a float or NA.
      or_beta = parse_or_beta(or_beta)

      # Make citation for this entry
      year = parse_year(date)
      citation = "%s et al. %s %s" % (author,year,journal)

      # Do we have a risk allele for this variant?
      strongest_snp_and_risk = e[20]
      strongest_snp, risk_allele = parse_strongest_snp_risk(strongest_snp_and_risk)

      # Risk allele frequency (but only for the strongest SNP, remember)
      risk_al_freq = e[26].strip()
      if risk_al_freq == "NR" or risk_al_freq == "":
        risk_al_freq = float("nan")
      else:
        try:
          risk_al_freq = float(risk_al_freq)
        except:
          log.add_warning("Invalid risk allele frequency",text_type("trait {}, snps {}, freq was: {}").format(trait,rsids,risk_al_freq))
          risk_al_freq = float("nan")

      # Genes labeled for this region
      genes = ",".join(map(lambda x: x.strip(),e[13].split(",")))

      if "intergenic" in genes:
        genes = ",".join(map(lambda x: x.strip(),e[14].split(",")))

      # Is the GWAS p-value significant?
      try:
        log_pval = float(log_pval)
        if log_pval < GWAS_LOG_SIG:
          log.add_filter("Failed p-value threshold"," ".join(map(text_type,[rsids,trait,log_pval])))
          continue
      except ValueError:
        log.add_filter("Invalid log p-value"," ".join(map(text_type,[rsids,trait,log_pval])))
        continue

      # Is the trait not blank?
      if trait.strip() == "":
        log.add_filter("Missing trait"," ".join(map(text_type,[rsids,trait,log_pval])))
        continue

      # There can be multiple SNPs on the same line for the same trait.
      rsids = [i.strip() for i in rsids.split(";")]

      # Sometimes, SNPs are specified as a haplotype with "rs1:rs2"
      for i in range(len(rsids)):
        isnp = rsids[i]

        if isnp.startswith("rs") and ':' in isnp:
          # It's a haplotype.
          haplo_snps = [s.strip() for s in isnp.split(":")]
          rsids.extend(haplo_snps)
          rsids.pop(i)

      for rsid in rsids:
        if not rsid.startswith("rs"):
          log.add_filter("Variant is not rsID",rsid)
          continue

        # The GWAS catalog sometimes includes corrupted rsIDs with strange unicode characters
        # This should clear everything that isn't 'rs' or digits
        rsid = re.sub("[^rs0-9]+", "", rsid)

        # Find the position for this SNP.
        vrecord = swiss_db.rsid_to_variant(rsid)

        # If it didn't have a chrom/pos in the database, we can't use it.
        if vrecord is None:
          log.add_filter("Couldn't find position/alleles for rsID",rsid)
          continue

        chrom, pos = vrecord.chrom, vrecord.pos
        ref, alt = vrecord.ref, vrecord.alt
        epacts = "{}:{}_{}/{}".format(chrom,pos,ref,alt)
        chrpos = "{}:{}".format(chrom,pos)

        # If we've already seen this association, we don't need to print it.
        key = "%s_%s_%s_%s" % (rsid,chrom,pos,trait)
        if key in seen_trait_snps:
          log.add_filter("SNP already known for this trait",key)
          continue
        else:
          seen_trait_snps.add(key)

        # Is this the SNP for which we have a risk allele?
        if rsid == strongest_snp:
          risk_frq_out = "{:.2f}".format(risk_al_freq) if not math.isnan(risk_al_freq) else "NA"
          or_beta_out = "{:.2f}".format(or_beta) if not math.isnan(or_beta) else "NA"
        else:
          risk_frq_out = "NA"
          or_beta_out = "NA"

        pos_s = str(pos)
        log_pval_s = "{:.2f}".format(log_pval)

        # We need to decompose the record into multiple records, one for each allele if ref or alt is multi-allelic.
        ref_alleles = ref.split(",")
        alt_alleles = alt.split(",")
        if len(ref_alleles) > 1 and len(alt_alleles) > 1:
          log.add_filter("SNP had multiple ref and alt alleles, cannot parse", key)
          continue

        if len(ref_alleles) > 1:
          for al in ref_alleles:
            # The GWAS catalog should have listed both alleles that were tested, instead of just the risk allele. For
            # example, if there is 1 REF allele, multiple ALT alleles, and the risk allele is the REF... then which ALT allele
            # was tested against the risk allele for association?
            # So we have to blank out the risk allele if it doesn't match any of the alleles. Same code below for alt alleles.
            risk_al_out = risk_allele if risk_allele in (al, alt) and rsid == strongest_snp else "NA"
            final_row = u"\t".join([rsid,epacts,chrpos,chrom,pos_s,al,alt,trait,trait,log_pval_s,citation,risk_al_out,risk_frq_out,genes,or_beta_out])
            print(final_row,file=out)
        elif len(alt_alleles) > 1:
          for al in alt_alleles:
            risk_al_out = risk_allele if risk_allele in (al, alt) and rsid == strongest_snp else "NA"
            final_row = u"\t".join([rsid,epacts,chrpos,chrom,pos_s,ref,al,trait,trait,log_pval_s,citation,risk_al_out,risk_frq_out,genes,or_beta_out])
            print(final_row,file=out)

  print("")
  log.summary()
  return outpath

class MergeHistoryNode(object):
  __slots__ = ["rsid","parents","child","deleted"]

  def __init__(self):
    self.rsid = None
    self.parents = []
    self.child = None
    self.deleted = False

  def __str__(self):
    parent_rsids = ",".join([x.rsid for x in self.parents])
    if self.child is not None:
      child_rsid = self.child.rsid
    else:
      child_rsid = "No child"

    return "rsid: %s | parents: %s | child: %s | deleted? %s" % (self.rsid,parent_rsids,child_rsid,self.deleted)

  def __repr__(self):
    return "id: %i | " % id(self) + str(self)

class MergeHistory:
  def __init__(self,rs_merge_arch,snp_history,snp_build):
    """
    snp_build is the dbsnp/UCSC SNP build number, e.g. snp138, snp142, etc.
    fpath is the path to the RsMergeArch.bcp.gz file
    """

    self.rs_merge_arch = rs_merge_arch
    self.snp_history = snp_history
    self.snp_build = int(snp_build.replace("snp","").replace("b",""))

    self.rsids = {}
    self.deleted = set()
    self.refsnp_trans = {}

    self._parse()

  def find_current(self,rs):
    if isinstance(rs,MergeHistoryNode):
      source_node = rs
    else:
      source_node = self.rsids.get(rs)

    if source_node is None:
      return None

    target = source_node
    chain = [target]

    # Walk forward from source to sink node
    while 1:
      next_t = target.child
      if next_t is not None:
        target = next_t
        chain.append(target)
      else:
        break

    # Now walk backward to find first node that isn't deleted
    for i in range(len(chain)-1,0,-1):
      if not chain[i].deleted:
        return chain[i].rsid

  def iter_nodes(self):
    for node in itervalues(self.rsids):
      yield node

  def find_all_parents(self,rs):
    if isinstance(rs,MergeHistoryNode):
      start = self.rsids.get(rs.rsid)
    else:
      start = self.rsids.get(rs)

    if start is None:
      return None

    to_visit = []
    all_parents = []

    for p in start.parents:
      to_visit.append(p)

    while len(to_visit) > 0:
      node = to_visit.pop()

      all_parents.append(node.rsid)

      for p in node.parents:
        to_visit.append(p)

    return all_parents

  def _parse(self,):
    with gzip.open(self.snp_history,"rt") as infp:
      for line in infp:
        self.deleted.add("rs" + line.split()[0])

    if self.rs_merge_arch.endswith(".gz"):
      f = gzip.open(self.rs_merge_arch,"rt")
    else:
      f = open(self.rs_merge_arch)

    for line in f:
      e = line.rstrip().split("\t")
      build = int(e[2])

      if build > self.snp_build:
        continue

      rsid_high = "rs" + e[0]
      rsid_low = "rs" + e[1]

      node_high = self.rsids.get(rsid_high)
      if node_high is None:
        node_high = MergeHistoryNode()
        node_high.rsid = rsid_high
        node_high.deleted = rsid_high in self.deleted

      node_low = self.rsids.get(rsid_low)
      if node_low is None:
        node_low = MergeHistoryNode()
        node_low.rsid = rsid_low
        node_low.deleted = rsid_low in self.deleted

      node_high.child = node_low
      node_low.parents.append(node_high)

      self.rsids[rsid_high] = node_high
      self.rsids[rsid_low] = node_low

def create_db(vcf_path,db_path,rsmerge_path,snphistory_path,dbsnp_build,chunksize=500000):
  print("Parsing histories...")
  history = MergeHistory(rsmerge_path,snphistory_path,dbsnp_build)

  print("Creating database...")
  con = sqlite3.connect(db_path)

  con.execute("PRAGMA SYNCHRONOUS=OFF")
  con.execute("PRAGMA TEMP_STORE=MEMORY")
  con.execute("PRAGMA PAGE_SIZE=4096")
  con.execute("PRAGMA CACHE_SIZE=500000")

  con.execute("DROP TABLE IF EXISTS dbsnp")
  con.execute("CREATE TABLE dbsnp (rsid TEXT, chrpos TEXT, epacts TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT)")
  con.execute("DROP TABLE IF EXISTS trans")
  con.execute("CREATE TABLE trans (rs_orig TEXT, rs_current TEXT)")

  with gzip.open(vcf_path,"rt") as vcf, con:
    print("Parsing VCF: %s" % vcf_path)
    for chunk in partition_all(chunksize,vcf):
      processed = []
      for line in chunk:
        if line.startswith("#"):
          continue

        chrom, pos, vid, ref, alt = line.split("\t")[0:5]
        rsid = vid.split(":")[0]
        epacts = "{}:{}_{}/{}".format(chrom,pos,ref,alt)

        if not rsid.startswith("rs"):
          continue

        rsid_trans = history.find_current(rsid)
        chrpos = "{}:{}".format(chrom,pos)

        use_rsid = rsid if rsid_trans is None else rsid_trans

        data = (use_rsid,chrpos,epacts,chrom,pos,ref,alt)

        processed.append(data)

      con.executemany("INSERT INTO dbsnp (rsid,chrpos,epacts,chrom,pos,ref,alt) VALUES (?,?,?,?,?,?,?)",processed)

    # Create indexes for important columns
    print("Creating indexes...")
    indexes = [
      "CREATE INDEX idx_snp ON dbsnp (rsid)",
      "CREATE INDEX idx_chrpos ON dbsnp (chrpos)",
      "CREATE INDEX idx_chrom_and_pos ON dbsnp (chrom,pos)",
      "CREATE INDEX idx_epacts ON dbsnp (epacts)"
    ]

    for idx in indexes:
      print(idx)
      con.execute(idx)

    # Create SNP translation table
    print("Creating SNP translation table..")
    for chunk in partition_all(chunksize,history.iter_nodes()):
      rows = []
      for node in chunk:
        # Only want the "sink" nodes - the most recent rsids at the bottom of the tree
        if node.child is None:
          parents = history.find_all_parents(node)
          for p in parents:
            rows.append([p,node.rsid])

      con.executemany("INSERT INTO trans (rs_orig,rs_current) VALUES (?,?)",rows)

    print("Creating indexes...")
    indexes = [
      "CREATE INDEX idx_rs_orig ON trans (rs_orig)",
      "CREATE INDEX idx_rs_current ON trans (rs_current)"
    ]

    for idx in indexes:
      print(idx)
      con.execute(idx)

def get_settings():
  p = OptionParser()
  p.add_option("--gwas-cat",help="Build a gwas catalog file.",action="store_true",default=True)
  p.add_option("--db",help="Database name or path to existing database. If not given, name will default to genome_build.dbsnp_build.sqlite")
  p.add_option("--no-cleanup",help="Leave temporary files alone after creating database instead of deleting them.",default=False,action="store_true")
  p.add_option("--genome-build",help="Specify human genome build to use when downloading data, e.g. GRCh37p13",default="GRCh37p13")
  p.add_option("--dbsnp-build",help="Specify dbSNP build to use, e.g. b147",default="b147")
  p.add_option("--rs-merge-arch",help="(optional) Specify already downloaded RsMergeArch file")
  p.add_option("--snp-history",help="(optional) Specify already downloaded SNP history file")

  opts, args = p.parse_args()

  if not opts.dbsnp_build.startswith("b"):
    raise ValueError("Unrecognized dbSNP build, must begin with 'b': " + opts.dbsnp_build)

  if not opts.genome_build.startswith("GRCh"):
    raise ValueError("Unrecognized human genome build, must begin with 'GRCh': " + opts.genome_build)

  if opts.db is None:
    opts.db = "{}_{}.sqlite".format(opts.genome_build,opts.dbsnp_build)

  if opts.rs_merge_arch is not None and not os.path.isfile(opts.rs_merge_arch):
    raise IOError("Can't find file: " + opts.rs_merge_arch)

  if opts.snp_history is not None and not os.path.isfile(opts.snp_history):
    raise IOError("Can't find file: " + opts.snp_history)

  return opts, args

def mkpath(tmpdir,filename):
  return os.path.join(tmpdir,filename)

def main():
  opts, args = get_settings()
  temp_files = []

  print("Genome build: " + opts.genome_build)
  print("dbSNP build: " + opts.dbsnp_build)

  db_name = opts.db
  if not os.path.isfile(db_name):
    if opts.rs_merge_arch is None:
      print("Downloading rsID merge history...")
      rsmerge_path = download_merge_arch()
    else:
      rsmerge_path = opts.rs_merge_arch

    if opts.snp_history is None:
      print("Downloading rsID deletion history...")
      snphistory_path = download_snp_history()
    else:
      snphistory_path = opts.snp_history

    vcf_url = NCBI_VCF_TEMPLATE_URL.format(opts.dbsnp_build,opts.genome_build)
    print("Downloading NCBI dbSNP VCF for {} / {} @ {}".format(opts.dbsnp_build,opts.genome_build,vcf_url))
    vcf_path = download_dbsnp_vcf(opts.dbsnp_build,opts.genome_build)

    # Create sqlite database
    create_db(vcf_path,db_name,rsmerge_path,snphistory_path,opts.dbsnp_build)

    # Write a file so we know when the database was created.
    db_info = opts.db + ".info"
    with open(db_info,'w') as info_out:
      print(time.strftime("Database created at %H:%M:%S %Z on %B %d %Y"),file=info_out)

    print("\nDatabase successfully created: %s" % db_name)
  else:
    print("Skipping database creation, already exists: %s" % db_name,file=sys.stderr)

  # Should we also try to build a GWAS catalog?
  # Currently EBI only releases GRCh38 datasets
  tmp_gwascat = "from_ebi_gwascatalog.txt"
  if opts.gwas_cat:
    # Download the catalog
    print("Downloading GWAS catalog..")
    download_ebi_catalog(GWAS_CAT_EBI_URL,tmp_gwascat)

    # Do some filtering and reformatting, and looking up positions for this genome build
    print("\nParsing/reformatting GWAS catalog..")
    final_gwas_cat = "gwascat_ebi_%s.tab" % opts.genome_build
    parse_gwas_catalog(tmp_gwascat,db_name,final_gwas_cat)

    print("\nCreated GWAS catalog: %s" % final_gwas_cat)

    temp_files.append(tmp_gwascat)

  # Delete all of the temporary files/directories we created.
  if not opts.no_cleanup:
    for f in temp_files:
      remove_file(f)

if __name__ == "__main__":
  main()

