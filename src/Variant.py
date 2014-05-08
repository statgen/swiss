#!/usr/bin/env python
import re

# Filter a list of variants down to only those variants
# with a unique chr/pos. 
# Variants without chr/pos are dropped. 
# Be careful with this and indels! A SNP and an indel
# can have the same position.  
def unique_by_pos(variant_list):
  final = [];
  seen = set();
  for v in variant_list:
    if v.chrpos is None:
      continue;
    
    has = v.chrpos in seen;
    if has:
      continue;
    else:
      seen.add(v.chrpos);
      final.append(v);

  return final;

class Variant:
  def __init__(self,string=None,pos_callable=None):
    self.name = None;
    self.chrpos = None;
    self.chrom = None;
    self.pos = None;
    self.build = None;

    # Reference and alternative alleles. 
    # ref_al should be only 1 allele, e.g. G or GCAT
    # alt_al can be multiple alleles, comma-sep, e.g. T,TCA,TCATAGCTACGATAAT
    self.ref_al = None;
    self.alt_al = None;

    self.pos_callable = pos_callable;

    if string is not None:
      self.from_str(string);

  # chr4:939191
  def from_chrpos(self,chrpos):
    res = re.search("chr(.+?)\:(\d+)",chrpos);
    if res is not None:
      self.chrom = res.groups()[0];
      self.pos = int(res.groups()[1]);
      self.chrpos = "%s:%s" % (self.chrom,self.pos);
      self.name = self.chrpos;
    else:
      raise ValueError, "Invalid chrpos name: %s" % chrpos;

  # rs914141
  def from_rsid(self,rsid):
    self.name = rsid;

    if self.pos_callable is not None:
      (chrom,pos) = pos_callable(rsid);

      if chrom is not None and pos is not None:
        self.chrom = chrom;
        self.pos = pos;
        self.chrpos = "%s:%s" % (chrom,pos);

  # 11-1414141
  def from_dash(self,dash):
    res = re.search("(.+?)\-(\d+)",dash);
    if res is not None:
      self.chrom = res.groups()[0];
      self.pos = int(res.groups()[1]);
      self.chrpos = "%s:%s" % (self.chrom,self.pos);
      self.name = self.chrpos;
    else:
      raise ValueError, "Cannot parse 'chr-pos' style Variant name: %s" % chrpos;

  # Try to figure out format and call appropriate method (above)
  def from_str(self,string):
    if all([x in string for x in [':','chr']]):
      self.from_chrpos(string);
    elif re.search("(.+?)\-(\d+)",string) is not None:
      self.from_dash(string);
    elif re.search("rs(\d+)",string) is not None:
      self.from_rsid(string);
    else:
      # Doesn't seem to have a known format - 
      # maybe this is a weird Variant name like exm-19141?
      self.rsid = string;

  # Return this variant as an EPACTS marker ID. 
  def as_epacts(self):
    return "{chrom}:{pos}_{ref}/{alt}".format(
      chrom = self.chrom,
      pos = self.pos,
      ref = self.ref_al,
      alt = self.alt_al
    );

  def as_chrpos(self):
    return "{chrom}:{pos}".format(
      chrom = self.chrom,
      pos = self.pos
    );
  
  def __str__(self):
    return "%s | %s | %s | %s/%s" % (self.name,self.chrom,self.pos,self.ref_al,self.alt_al);

  def __repr__(self):
    return "%s | %s | %s | %s/%s @ id:%s" % (self.name,self.chrom,self.pos,self.ref_al,self.alt_al,id(self));


