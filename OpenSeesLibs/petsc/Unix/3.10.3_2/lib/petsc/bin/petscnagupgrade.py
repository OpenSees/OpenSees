#!/usr/bin/env python
#!/bin/env python
#
#    Nags the user to update to the latest version
#
from __future__ import print_function
import os
import os.path, time,sys
import re
from distutils.version import LooseVersion as Version

def naggedtoday(file):
  if not os.path.exists(file): return 0
  if time.time() - os.path.getmtime(file) > 60*60*24: return 0
  return 1

def parse_version_h(pv):
  major    = int(re.compile(' PETSC_VERSION_MAJOR[ ]*([0-9]*)').search(pv).group(1))
  minor    = int(re.compile(' PETSC_VERSION_MINOR[ ]*([0-9]*)').search(pv).group(1))
  subminor = int(re.compile(' PETSC_VERSION_SUBMINOR[ ]*([0-9]*)').search(pv).group(1))
  if subminor != 0:           # Maintenance releases are numbered x.y.z
    return Version('%d.%d.%d' % (major, minor, subminor))
  else:                         # Feature releases are x.y
    return Version('%d.%d' % (major, minor))

def currentversion(petscdir):
  try:
    fd  = open(os.path.join(petscdir, 'include', 'petscversion.h'))
    pv = fd.read()
    fd.close()
    version = parse_version_h(pv)
  except:
    return
  try:
    import urllib2
    fd = urllib2.urlopen("https://bitbucket.org/petsc/petsc/raw/maint/include/petscversion.h",timeout = 2)
    #fd = urllib2.urlopen("http://www.mcs.anl.gov/petsc/petsc-current/include/petscversion.h",timeout = 2)
    pv = fd.read()
    fd.close()
    aversion = parse_version_h(pv)
  except:
    return
  if aversion > version:
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print("The version of PETSc you are using is out-of-date, we recommend updating to the new release")
    print(" Available Version: "+str(aversion)+"   Installed Version: "+str(version))
    print("http://www.mcs.anl.gov/petsc/download/index.html")
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  try:
    fd = open(os.path.join(petscdir,'.nagged'),"w")
    fd.close()
  except:
    return

  return 0
#
#
if __name__ ==  '__main__':
  if 'PETSC_DIR' in os.environ:
    petscdir = os.environ['PETSC_DIR']
  elif os.path.exists(os.path.join('.', 'include', 'petscversion.h')):
    petscdir  = '.'
  else:
    sys.exit(0)
  file     = os.path.join(petscdir,'.nagged')
  if not naggedtoday(file):
    currentversion(petscdir)


