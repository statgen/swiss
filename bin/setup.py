#!/usr/bin/env python

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

import os, sys, stat
import os.path as path

# Identify where we're located.
this_dir = path.dirname(path.realpath(sys.argv[0]))
root_dir = path.join(this_dir,"..")

# Create the virtualenv
# The virtualenv script automatically installs the required packages
env_dir = path.join(this_dir,"../env")
env_script = path.join(this_dir,"virtualenv-swiss.py")
env_cmd = "%s %s" % (env_script,env_dir)
os.system(env_cmd)

# Set permissions so that anyone can read the libraries installed in the virtualenv
for root, dirs, files in os.walk(env_dir):
  os.chmod(root,os.stat(root).st_mode | stat.S_IROTH | stat.S_IXOTH)

  for f in files:
    fpath = path.join(root,f)

    if path.islink(fpath):
      continue

    os.chmod(fpath,os.stat(fpath).st_mode | stat.S_IROTH)

# Let's do a quick check to make sure we can import things from the virtualenv
activate_this = path.join(env_dir,'bin/activate_this.py')

if not os.path.isfile(activate_this):
  sys.exit("Error: no virtualenv found! Please contact the developer with any errors logged to your "
           "console at: https://github.com/welchr/Swiss/issues")

# Push the virtualenv onto the path
execfile(activate_this,dict(__file__=activate_this))

try:
  import pandas
except:
  print >> sys.stderr, "Error: tried to import pandas as a test, but failed. Please contact the developer with " \
                       "any errors logged to your console at: https://github.com/welchr/Swiss/issues"

  raise
