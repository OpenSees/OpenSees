#!/usr/bin/env python
#
#    Reads in the output trajectory data -ts_save_trajectory -ts_trajectory_type visualization
#    Also provides a way to plot the trajectory
#
#  Make sure $PETSC_DIR/bin is in your PYTHONPATH
#
from __future__ import print_function
import os
import PetscBinaryIO
import numpy as np
import matplotlib.pyplot as pyplot

def ReadTrajectory(directory):
  io = PetscBinaryIO.PetscBinaryIO()

  v = []
  t = []

  cnt = 0
  while 1:
    try:
      fh  = open(os.path.join(directory,'SA-%06d.bin'%cnt));
      cnt = cnt + 1
      objecttype = io.readObjectType(fh)
      v.append(io.readVec(fh))
      t.append(np.fromfile(fh, dtype=io._scalartype, count=1)[0])
    except:
      break

  names = []
  try:
    fh  = open(os.path.join(directory,'variablenames'))
    nstrings = np.fromfile(fh, dtype=io._inttype, count=1)[0]
    sizes = np.fromfile(fh, dtype=io._inttype, count=nstrings)
    for i in range(0,nstrings):
      s = np.fromfile(fh, dtype=np.byte, count=sizes[i])
      names.append("".join(map(chr, s))[0:-1])
  except:
    pass

  return (t,v,names)

def PlotTrajectories(t,v,names,subnames):
  print(names)
  sub = []
  for s in subnames:
    sub.append(names.index(s))
  w = []
  for i in v:
    w.append(i[sub])

  pyplot.plot(t,w)
  pyplot.legend(subnames)
  pyplot.show()

