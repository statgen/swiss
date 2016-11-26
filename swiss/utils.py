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

import __builtin__
import os, sys, decimal, fnmatch, gzip
from termcolor import colored

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

def call_ipdb():
  import sys
  sys.stdout = sys.__stdout__;
  sys.stderr = sys.__stderr__;

  import ipdb
  ipdb.set_trace();

def is_gzip(file,verbose=False):
  import gzip

  try:
    f = gzip.open(file);
  except:
    if verbose:
      print "(gzip) could not open file, not gzipped: %s" % file;

    return False;

  try:
    f.read(1024);
    b = True;
  except:
    b = False;
    if verbose:
      print "(gzip) could not read data from gzip file object: %s" % file;
  finally:
    f.close();

  return b;

def error(msg):
  if hasattr(__builtin__,'SWISS_DEBUG') and __builtin__.SWISS_DEBUG:
    raise Exception, msg;
  else:
    print >> sys.stderr, colored("Error: ",'red') + msg;
    sys.exit(1);

def warning(msg):
  print >> sys.stderr, colored("Warning: ",'yellow') + msg;

# Great way to implement python enums, courtesy of: 
# http://stackoverflow.com/questions/36932/whats-the-best-way-to-implement-an-enum-in-python
def enum(*sequential, **named):
  enums = dict(zip(sequential, range(len(sequential))), **named)
  reverse = dict((value, key) for key, value in enums.iteritems())
  enums['reverse'] = reverse
  return type('Enum', (), enums)

def singleton(cls):
  instances = {}
  def getinstance():
      if cls not in instances:
          instances[cls] = cls()
      return instances[cls]
  return getinstance

def which(f):
  for path in os.environ['PATH'].split(os.pathsep):
    if not os.path.exists(path):
      continue;

    if not os.path.isdir(path):
      continue;
    
    for file in os.listdir(path):
      if os.path.basename(f) == file:
        return os.path.join(path,file);

  return None;

def find_systematic(fpath):
  """
  Tries to find a file in the following order:
  1) the given file
  2) on the user's path
  """

  if fpath is None:
    return None

  if os.path.isfile(fpath):
    return os.path.abspath(fpath)
  
  whiched_file = which(fpath)
  if whiched_file is not None:
    return whiched_file
 
  return None

def is_number(x):
  try:
    decimal.Decimal(str(x));
  except:
    return False;

  return True;

# Terminates program with only a message.
def die(msg):
  print >> sys.stderr, msg;
  sys.exit(1);

# Tests to see if interval i1 is completely contained within interval i2.
# i1 should be a list or tuple [start,end], same with i2
def interval_contained(i1,i2):
  if i1[0] >= i2[0] and i1[1] <= i2[1]:
    return True;
  else:
    return False;

def chrom2chr(c):
  c = str(c);
  if c in ('X','chrX','chromX'):
    return 23;
  elif c in ('Y','chrY','chromY'):
    return 24;
  elif c in ('MT','mito','chrM','chrMT'):
    return 25;
  elif c in ('XY','chrXY'):
    return 26;
  else:
    try:
      c = c.replace("chrom","").replace("chr","")
      c = int(float(c));
    except:
      c = None;

    return c;

def chr2chrom(c):
  c = int(c);
  if c < 23 and c > 0:
    return str(c);
  elif c == 23:
    return 'X';
  elif c == 24:
    return 'Y';
  elif c == 25:
    return 'mito';
  elif c == 26:
    return 'XY';
  else:
    return None;
  
