#!/usr/bin/env python

from Pegasus.DAX3 import *
import sys
import os
import optparse


# --- global variables ----------------------------------------------------------------

prog_dir  = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0])))
prog_base = os.path.split(sys.argv[0])[1]   # Name of this program

opensees = None
bag_files = dict()
bag_jobs = dict()


# --- functions -----------------------------------------------------------------------


def add_opensees_job(dax):
   
    # input files to the DAX-level replica catalog
    input_files = list()
    names = ["Example3.1.tcl", "Example3.2.tcl"]

    for name in names:
        ifile = File(name)
        ifile.addPFN(PFN("file://" + os.getcwd() + "/inputs/" + name, "local"))
        dax.addFile(ifile)
        input_files.append(ifile)

    # job - materialProperties.tcl.XX model.tcl motion.tcl.YY dynamic WW
    job = Job(namespace="opensees", name="opensees", version="1.0")
    job.addArguments("Example3.2.tcl")

    for f in input_files:
        job.uses(f, link=Link.INPUT)

    out_name = "Node.out"
    f1 = File(out_name)
    job.uses(f1, link=Link.OUTPUT)
    bag_files[out_name] = f1

    dax.addJob(job)

    # add the job to the bag so we can look it up later
    bag_jobs["opensees"] = job
    

def add_plot_job(dax):
   
    # genPlot.m
    genPlot = File("genPlot.m")
    genPlot.addPFN(PFN("file://" + os.getcwd() + "/inputs/genPlot.m", "local"))
    dax.addFile(genPlot)

    # job
    job = Job(namespace="system", name="octave", version="1.0")

    job.uses(genPlot, link=Link.INPUT)

#    job.addArguments("--silent", "--eval", "genPlot;quit")
    job.addArguments("--silent","genPlot.m")

    node_file = bag_files["Node.out"]
    job.uses(node_file, link=Link.INPUT)

    f = File("Disp.jpg")
    job.uses(f, link=Link.OUTPUT)

    # add the file to the bag so we can look it up later
 #   bag_files[out_name] = f
            
    dax.addJob(job)

    dax.depends(parent=bag_jobs["opensees"], child=job)


# --- main ----------------------------------------------------------------------------

# Configure command line option parser
prog_usage = "usage: %s [options]" % (prog_base)
parser = optparse.OptionParser(usage=prog_usage)
parser.add_option("-p", "--num-mat-props", action = "store", dest = "num_mat_props",
                  type="int", default = "1",
                  help = "The number of material properties files")
parser.add_option("-m", "--num-motions", action = "store", dest = "num_motions",
                  type="int", default = "1",
                  help = "The number of motions")
(options, args) = parser.parse_args()

# Create a abstract workflow
dax = ADAG("dax")

# Add executables to the DAX-level replica catalog
opensees = Executable(namespace="opensees", name="opensees", version="1.0", os="linux", arch="x86_64", installed=True)

#opensees.addPFN(PFN("file:///apps/opensees/bin/OpenSees", "local"))
opensees.addPFN(PFN("file:///grid/app/nees/opensees/bin/OpenSees", "HCC_US_Fermigridosg1"))
dax.addExecutable(opensees)

octave = Executable(namespace="system", name="octave", version="1.0", os="linux", arch="x86_64", installed=True)
#octave.addPFN(PFN("file:///usr/bin/octave", "local"))
octave.addPFN(PFN("file:///grid/app/nees/bin/octave-3.2.4.sh", "HCC_US_Fermigridosg1"))
dax.addExecutable(octave)

add_opensees_job(dax)
add_plot_job(dax)

# Write the DAX to stdout
dax.writeXML(sys.stdout)



