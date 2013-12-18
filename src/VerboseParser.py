#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2010 Ryan Welch
# 
# Snipper is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Snipper is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

from optparse import OptionParser, SUPPRESS_HELP
from textwrap import *
from terminalsize import *
import sys

# Slightly modified OptionParser to make the help output look reasonable.
class VerboseParser(OptionParser):
  def __init__(self,*args,**kwargs):
    OptionParser.__init__(self,*args,**kwargs);

  def print_help(self):
    # Print usage.
    print "usage: %s" % self.usage;
    print "";

    console_width = 120;
    try:
      console_width = get_terminal_size(120)[0];
    except:
      pass
    
    # Print options.
    for option in self.option_list:
      if option.help == "SUPPRESSHELP":
        continue;

      if option.type != None:
        print fill(", ".join(option._short_opts + option._long_opts) + " <%s>" % option.type,\
          initial_indent="  ",subsequent_indent="  ",width=console_width);
      else:
        print fill(", ".join(option._short_opts + option._long_opts),\
          initial_indent="  ",subsequent_indent="  ",width=console_width);

      for line in option.help.split("\n"):
        print fill(line,initial_indent="    ",subsequent_indent="    ",width=console_width);
      print "";

