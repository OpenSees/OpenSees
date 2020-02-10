#!/usr/bin/env python
#!/bin/env python

# Usage:
#  taucc -cc=g++ -pdt_parse=/../cxxparse -tau_lib_dir=/.../lib COMPILE_OPTIONS
#
#  Options:
#           -cc              : C/C++ compiler
#           -pdt_parse       : pdtoolkit parser for C++
#           -tau_lib_dir     : TAU library dir
#           -v,-verbose      : verbose mode - shows the exact commands invoked
#           -leave_tmp       : do not delete temporary files
#           -E               : run preprocessor only
#
from __future__ import print_function
import commands
import sys
import os
import string
import tempfile
def runcmd(cmd,verbose):
  if verbose:
    print(cmd)
  (status, output) = commands.getstatusoutput(cmd)
  if status:
    raise RuntimeError('Unable to run '+cmd+':\n'+output)
  elif output:
    print(output)

def getTauFlags(tau_lib_dir):
  fd,name=tempfile.mkstemp(prefix='taucc-')
  try:
    tau_makefile = str(os.path.join(tau_lib_dir, os.path.basename(os.environ['TAU_MAKEFILE'])))
  except:
    tau_makefile = str(os.path.join(tau_lib_dir, 'Makefile.tau-mpi-pdt'))
  buf  = 'include ' + tau_makefile
  buf += '''
tauflags:
	-@echo TAU_INSTRUMENT0R:${TAU_PREFIX_INSTALL_DIR}/${CONFIG_ARCH}/bin/tau_instrumentor
	-@echo TAU_DEFS:${TAU_DEFS}
	-@echo TAU_INCLUDE:${TAU_INCLUDE}
	-@echo TAU_LIBS:${TAU_LIBS}
	-@echo TAU_MPI_INC:${TAU_MPI_INC}
	-@echo TAU_MPI_LIBS:${TAU_MPI_LIBS}
	-@echo TAU_CXXLIBS:${TAU_CXXLIBS}
'''
  os.write(fd,buf)
  os.close(fd)
  cmd = 'make -f '+ name + ' tauflags'
  (status, output) = commands.getstatusoutput(cmd)
  if status:
    os.remove(name)
    raise RuntimeError('Unable to run '+cmd+':\n'+output)
  os.remove(name)
  # remove unnecessary stuff from output
  tau_defs =''
  tau_include=''
  tau_libs=''
  tau_mpi_libs=''
  for line in output.splitlines():
    if line.find('TAU_INSTRUMENT0R:') >= 0:  tau_instr = line.replace('TAU_INSTRUMENT0R:','')
    elif line.find('TAU_DEFS:') >= 0:  tau_defs = line.replace('TAU_DEFS:',' ')
    elif line.find('TAU_INCLUDE:') >= 0: tau_include = line.replace('TAU_INCLUDE:',' ')
    elif line.find('TAU_LIBS:') >= 0: tau_libs = line.replace('TAU_LIBS:',' ')
    elif line.find('TAU_MPI_INC:') >= 0: tau_mpi_inc = line.replace('TAU_MPI_INC:',' ')
    elif line.find('TAU_MPI_LIBS:') >= 0: tau_mpi_libs = line.replace('TAU_MPI_LIBS:',' ')
    elif line.find('TAU_CXXLIBS:') >= 0: tau_cxxlibs = line.replace('TAU_CXXLIBS:',' ')
  return tau_instr,tau_defs,tau_include,tau_libs,tau_mpi_inc,tau_mpi_libs,tau_cxxlibs

def main():

  sourcefiles=[]
  objfiles=[]
  libfiles=[]
  arglist=''
  pdt_parse='cxxparse'
  iscparse=0
  tau_instr='tau_instrumentor'
  cc='gcc'
  verbose=0
  compileonly=0
  leave_tmp = 0

  for arg in sys.argv[1:]:
    filename,ext = os.path.splitext(arg)
    argsplit =  arg.split('=')
    #look for sourcefiles, validate & add to a list
    if ext == '.c' or ext == '.cc' or ext == '.cpp' or ext == '.cxx' or ext == '.C':
      if os.path.isfile(arg):
        sourcefiles.append(arg)
    elif argsplit[0] == '-cc':
      cc = argsplit[1]
    elif argsplit[0] == '-pdt_parse':
      pdt_parse = argsplit[1]
      if pdt_parse.find('cparse') >=0: iscparse = 1
      elif pdt_parse.find('cxxparse') >=0: iscparse = 0
      else: sys.exit('Error: Unknown parser - use either cparse or cxxparse: '+pdt_parse)
    elif arg == '-c':
        compileonly = 1
    elif arg == '-E':
        compileonly = 1
        arglist += ' '+arg
    elif arg == '-leave_tmp':
      leave_tmp = 1
    elif argsplit[0] == '-tau_lib_dir':
      tau_lib_dir = argsplit[1]
    elif arg.startswith('-L') or arg.startswith('-l'):
      libfiles.append(arg)
    elif arg == '-v' or arg == '-verbose':
        verbose  = 1
        arglist += ' '+arg
    else:
      # Now make sure quotes are escaped properly
      # Group the rest of the arguments into a different list
      arg=arg.replace('"','\\"')
      arglist += ' '+arg

  # error if sourcefiles not specified with -c
  if sourcefiles == [] and compileonly:
    sys.exit('Error: no sourcefiles specified with -c')
  # obtain TAU info from TAU makefile
  tau_instr,tau_defs,tau_include,tau_libs,tau_mpi_inc,tau_mpi_libs,tau_cxxlibs = getTauFlags(tau_lib_dir)
  if sourcefiles != []:
    # Now Compile the sourcefiles
    for sourcefile in sourcefiles:
      root,ext = os.path.splitext(sourcefile)
      obj_file = root + '.o'
      objfiles.append(obj_file) # for use later at linktime
      if iscparse:
        if ext == '.c': pdt_file = root+ '.pdb'
        else: pdt_file = sourcefile+ '.pdb'
      else:
        if ext == '.cc' or ext == '.cpp' or ext == '.cxx' or ext == '.C' : pdt_file = root+ '.pdb'
        else: pdt_file = sourcefile+ '.pdb'
      tau_file = root +'.inst' + ext
      cmd1  = pdt_parse + ' ' + sourcefile + arglist + tau_defs + tau_include + tau_mpi_inc

      cmd2  = tau_instr + ' ' + pdt_file + ' ' + sourcefile +' -o '+ tau_file
      cmd2 += ' -c -rn PetscFunctionReturn -rv PetscFunctionReturnVoid\\(\\)'
      cmd3  = cc + ' -c ' + tau_file + ' -o ' + obj_file + arglist + tau_defs + tau_include +tau_mpi_inc


      runcmd(cmd1,verbose)
      runcmd(cmd2,verbose)
      if not leave_tmp: os.remove(pdt_file)
      runcmd(cmd3,verbose)
      if not leave_tmp: os.remove(tau_file)

  if not compileonly:
    objarg=''
    for objfile in objfiles:
      objarg += ' ' + objfile
    libarg=''
    for libfile in libfiles:
      libarg += ' '+libfile
    cmd1  = cc + ' ' + objarg +' '  + arglist + libarg + tau_mpi_libs + tau_libs + tau_cxxlibs
    runcmd(cmd1,verbose)
    if not leave_tmp: # delete the objfiles created
      for objfile in objfiles:
        os.remove(objfile)

if __name__ ==  '__main__':
  try:
    main()
  except Exception as e:
    sys.exit('ERROR: '+str(e))
