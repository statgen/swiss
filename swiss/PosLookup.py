#!/usr/bin/env python
import sqlite3
import os
import sys
from utils import *

def parse1000G(snp):
  if snp is None:
    return None;

  c = snp.split(":");
  if len(c) == 2:
    chr = chrom2chr(c[0][3:]);
    try:
      pos = int(c[1]);
    except:
      return None;
    
    return (chr,int(c[1]));
  else:
    return None;

class PosLookup:
  def __init__(self,db_file):  
    if not os.path.isfile(db_file):
      sys.exit("Error: could not locate SQLite database file: %s. Check conf file setting SQLITE_DB." % db_file);
      
    self.db = sqlite3.connect(db_file);
    self.execute = self.db.execute;
    
    self.execute("""
      CREATE TEMP VIEW snp_pos_trans AS SELECT rs_orig as snp,chr,pos FROM %s p INNER JOIN %s t ON (t.rs_current = p.snp);
    """ % (SQLITE_SNP_POS,SQLITE_TRANS));
    
    self.query = """
      SELECT snp,chr,pos FROM snp_pos_trans WHERE snp='%s';
    """;

  def __call__(self,snp):
    snp = str(snp);
    
    # If the SNP is a 1000G SNP, it already knows its chrom/pos by definition,
    # i.e. the SNP will be chr4:91941. 
    gcheck = parse1000G(snp);
    if gcheck:
      return gcheck;
    
    cur = self.execute(self.query % snp);
    chr = None;
    pos = None;

    res = 0;
    for row in cur:
      chr = row[1];
      pos = row[2];
      res += 1;

    region = "chr%s:%s" % (chr,pos);
    if res > 1:
      print >> sys.stderr, "Warning: SNP %s has more than 1 position in database, using: %s" % (str(snp),region);

    return (chr,pos);
