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

import re

# Filter a list of variants down to only those variants
# with a unique chr/pos. 
# Variants without chr/pos are dropped. 
# Be careful with this and indels! A SNP and an indel
# can have the same position.  
def unique_by_pos(variant_list):
  final = []
  seen = set()
  for v in variant_list:
    if v.chrpos is None:
      continue
    
    has = v.chrpos in seen
    if has:
      continue
    else:
      seen.add(v.chrpos)
      final.append(v)

  return final

# Try to parse an EPACTS ID into components.
# Returns None if it couldn't be parsed.
# Otherwise a tuple of:
# (chrom, pos, ref, alt, extra)
def parse_epacts(v):
  split = v.split("_")

  # It should have at least 2 elements (chrom:pos, ref/alt)
  if len(split) < 2:
    return None

  # Split chrom/pos
  try:
    chrom, pos = split[0].split(":")
  except:
    return None

  # Position should be numeric
  try:
    long(pos)
  except:
    return None

  # Split the alleles
  try:
    ref, alt = split[1].split("/")
  except:
    return None

  return [chrom,pos,ref,alt,"_".join(split[2:])]

# Try to identify variant as SNP from EPACTS ID.
# Returns:
# -- True if definitely a SNP
# -- False if definitely not a SNP
# -- None if can't safely determine
def is_snp_epacts_heuristic(v):
  parsed = parse_epacts(v)
  if parsed is None:
    return None

  ref, alt = parsed[2:4]

  if len(ref) > 1:
    return False
  elif len(alt) > 1:
    return False

  return True

class Variant:
  def __init__(self,string=None,pos_callable=None):
    self.name = None
    self.chrpos = None
    self.chrom = None
    self.pos = None
    self.build = None
    self.traits = []

    # Reference and alternative alleles. 
    # ref_al should be only 1 allele, e.g. G or GCAT
    # alt_al can be multiple alleles, comma-sep, e.g. T,TCA,TCATAGCTACGATAAT
    self.ref_al = None
    self.alt_al = None

    self.pos_callable = pos_callable

    if string is not None:
      self.from_str(string)

  # chr4:939191
  def from_chrpos(self,chrpos):
    res = re.search("chr(.+?)\:(\d+)",chrpos)
    if res is not None:
      self.chrom = res.groups()[0]
      self.pos = int(res.groups()[1])
      self.chrpos = "%s:%s" % (self.chrom,self.pos)
      self.name = self.chrpos
    else:
      raise ValueError, "Invalid chrpos name: %s" % chrpos

  # rs914141
  def from_rsid(self,rsid):
    self.name = rsid

    if self.pos_callable is not None:
      (chrom,pos) = pos_callable(rsid)

      if chrom is not None and pos is not None:
        self.chrom = chrom
        self.pos = pos
        self.chrpos = "%s:%s" % (chrom,pos)

  # 11-1414141
  def from_dash(self,dash):
    res = re.search("(.+?)\-(\d+)",dash)
    if res is not None:
      self.chrom = res.groups()[0]
      self.pos = int(res.groups()[1])
      self.chrpos = "%s:%s" % (self.chrom,self.pos)
      self.name = self.chrpos
    else:
      raise ValueError, "Cannot parse 'chr-pos' style Variant name: %s" % chrpos

  # Try to figure out format and call appropriate method (above)
  def from_str(self,string):
    if all([x in string for x in [':','chr']]):
      self.from_chrpos(string)
    elif re.search("(.+?)\-(\d+)",string) is not None:
      self.from_dash(string)
    elif re.search("rs(\d+)",string) is not None:
      self.from_rsid(string)
    else:
      # Doesn't seem to have a known format - 
      # maybe this is a weird Variant name like exm-19141?
      self.rsid = string

  # Return this variant as an EPACTS marker ID. 
  def as_epacts(self):
    return "{chrom}:{pos}_{ref}/{alt}".format(
      chrom = self.chrom,
      pos = self.pos,
      ref = self.ref_al,
      alt = self.alt_al
    )

  def as_chrpos(self):
    return "{chrom}:{pos}".format(
      chrom = self.chrom,
      pos = self.pos
    )
  
  def __str__(self):
    return "%s | %s | %s | %s/%s" % (self.name,self.chrom,self.pos,self.ref_al,self.alt_al)

  def __repr__(self):
    return "%s | %s | %s | %s/%s @ id:%s" % (self.name,self.chrom,self.pos,self.ref_al,self.alt_al,id(self))


