#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2010 Ryan Welch, Randall Pruim
# 
# LocusZoom is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# LocusZoom is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import shelve
import sys
from FileLock import *

LD_CACHE_PROTOCOL = "1.3";

class LDRegionCache():
  def __init__(self,key,dbfile):
    self.dbfile = dbfile;
    self.key = key;
    self.root = None;

    # Load database.
    existed = False;
    if not os.path.isfile(dbfile):
      print >> sys.stderr, "Warning: could not locate LD cache database, tried: %s, creating new one.." % dbfile;
    else:
      print >> sys.stderr, "Accessing LD cache: %s" % dbfile;
      existed = True;
    
    try:
      self.db = shelve.open(dbfile);
      self.opened = True;
    except:
      self.opened = False;
      raise;

    # Check protocol. 
    proto = self.db.get('proto');
    if proto == None or proto != LD_CACHE_PROTOCOL:
      if existed:
        print >> sys.stderr, "Error: LD cache was created with an older version"\
        " of the program. Please delete '%s', or specify a new cache file, and re-run." % dbfile;
        sys.exit(1);
      else:
        self.db['proto'] = LD_CACHE_PROTOCOL;

    try:
      self.root = self.db.get(key);
    except:
      print >> sys.stderr, "Error: accessing LD cache failed. Please contact admin with following info: ";
      print >> sys.stderr, "Cache location: %s" % dbfile;
      print >> sys.stderr, "Available keys: "
      for key in self.db.keys():
        print str(key);
      print >> sys.stderr, "Protocol version: %s" % self.db.get('proto');
      sys.exit(1);
      
    if self.root == None:
      print >> sys.stderr, "Creating new LD cache branch for current LD source/population/build..";

      new_root = {};
      self.db[key] = new_root;
      self.root = new_root;

  def getLD(self,refsnp,othersnp):
    return self.root.get(refsnp).get('ld').get(othersnp);

  def getAllLD(self,refsnp):
    return self.root.get(refsnp).get('ld');

  # Add LD for a region into the cache.
  # Start - start position of interval calculated
  # Stop - stop position of interval calculated
  # Data - dictionary of LD values, where keys are 'other' SNPs, values are LD
  def updateLD(self,refsnp,start,stop,data):
    if data != None and len(data) > 0:
      refsnp_root = self.root.setdefault(refsnp,{});
      refsnp_root['start'] = start;
      refsnp_root['stop'] = stop;
      refsnp_root['ld'] = data;

  def getStart(self,refsnp):
    return self.root.get(refsnp).get('start');

  def getStop(self,refsnp):
    return self.root.get(refsnp).get('stop');

  def hasRegion(self,refsnp,start,stop):
    data = self.root.get(refsnp);
    if data:
      if data.get('start') <= start and data.get('stop') >= stop:
        return True;

    return False;

  def close(self):
    if self.opened:
      lock = FileLock(self.dbfile);
      lock.acquire();

      # Update database.
      if self.root != None:
        self.db[self.key] = self.root;
        self.db.close();

      # Release lock. 
      lock.release();
    
      # Set object to closed state. 
      self.opened = False;

  def __del__(self):
    self.close();
