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

from optparse import OptionParser, SUPPRESS_HELP
from textwrap import *
from terminalsize import *
import sys

# Slightly modified OptionParser to make the help output look reasonable.
class VerboseParser(OptionParser):
  def __init__(self,*args,**kwargs):
    OptionParser.__init__(self,*args,**kwargs)

  def print_help(self):
    # Print usage.
    print "usage: %s" % self.usage
    print ""

    console_width = 120
    try:
      console_width = get_terminal_size(120)[0]
    except:
      pass
    
    # Print options.
    for option in self.option_list:
      if option.help == "SUPPRESSHELP":
        continue

      if option.type is not None:
        print fill(", ".join(option._short_opts + option._long_opts) + " <%s>" % option.type,\
          initial_indent="  ",subsequent_indent="  ",width=console_width)
      else:
        print fill(", ".join(option._short_opts + option._long_opts),\
          initial_indent="  ",subsequent_indent="  ",width=console_width)

      for line in option.help.split("\n"):
        print fill(line,initial_indent="    ",subsequent_indent="    ",width=console_width)
      
      if not option.default == ("NO","DEFAULT"):
        if option.type == "int":
          default_value = str(int(option.default))
        elif option.type in ('long','float'):
          default_value = "%0.4g" % option.default
        elif isinstance(option.default,str):
          if option.default == "\t":
            default_value = "tab"
          elif option.default == " ":
            default_value = "space"
          else:
            default_value = str(option.default)
        else:
          default_value = str(option.default)

        default_str = "Default value is: %s" % default_value
        print fill(default_str,initial_indent="    ",subsequent_indent="    ",width=console_width)

      print ""

