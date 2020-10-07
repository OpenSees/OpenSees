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


def add_opensees_jobs(dax, num_materials, num_motions):
   
    # input files to the DAX-level replica catalog
    input_files = list()
    names = ["model.tcl", "SteelWSections.tcl", "analyzeModel.tcl"]

    for name in names:
        ifile = File(name)
        ifile.addPFN(PFN("file://" + os.getcwd() + "/inputs/" + name, "local"))
        dax.addFile(ifile)
        input_files.append(ifile)

    # material properties files
    mat_file = dict()
    for i in range(1, num_materials + 1):
        name = "matProperties.tcl.%d" % (i)
        mat_file[i] = File(name)
        mat_file[i].addPFN(PFN("file://" + os.getcwd() + "/inputs/" + name, "local"))
        dax.addFile(mat_file[i])

    # motion files
    mot_tcl_file = dict()
    mot_dat_file = dict()
    for i in range(1, num_motions + 1):

	    # tcl
        name = "motion.tcl.%d" % (i)
        mot_tcl_file[i] = File(name)
        mot_tcl_file[i].addPFN(PFN("file://" + os.getcwd() + "/inputs/" + name, "local"))
        dax.addFile(mot_tcl_file[i])
	
        # dat
        name = "motion.dat.%d" % (i)
        mot_dat_file[i] = File(name)
        mot_dat_file[i].addPFN(PFN("file://" + os.getcwd() + "/inputs/" + name, "local"))
        dax.addFile(mot_dat_file[i])

    # generate the jobs
    for n in range(1, num_materials + 1):
        for m in range(1, num_motions + 1):

            w = "%d.%d" %(n, m)

            # job - materialProperties.tcl.XX model.tcl motion.tcl.YY dynamic WW
            job = Job(namespace="opensees", name="opensees", version="1.0")
            job.addArguments("analyzeModel.tcl", "%d" % (n), "%d" % (m))

            for f in input_files:
                job.uses(f, link=Link.INPUT)

            job.uses(mat_file[n], link=Link.INPUT)
            job.uses(mot_tcl_file[m], link=Link.INPUT)
            job.uses(mot_dat_file[m], link=Link.INPUT)
            
            out_name = "NodeDisp.out.%s" % (w)
            f1 = File(out_name)
            job.uses(f1, link=Link.OUTPUT)
            # add the file to the bag so we can look it up later
            bag_files[out_name] = f1

            out_name = "NodeAccel.out.%s" % (w)
            f2 = File(out_name)
            job.uses(f2, link=Link.OUTPUT)
            # add the file to the bag so we can look it up later
            bag_files[out_name] = f2

            out_name = "NodeReaction.out.%s" % (w)
            f3 = File(out_name)
            job.uses(f3, link=Link.OUTPUT)
            # add the file to the bag so we can look it up later
            bag_files[out_name] = f3

            out_name = "NodeDrift.out.%s" % (w)
            f4 = File(out_name)
            job.uses(f4, link=Link.OUTPUT)
            # add the file to the bag so we can look it up later
            bag_files[out_name] = f4

            dax.addJob(job)
            # add the job to the bag so we can look it up later
            bag_jobs["opensees_%s" % (w)] = job

            # parents to this job
            #dax.depends(parent=bag_jobs["motion_job"], child=job)
    

def add_plot_job(dax, num_materials, num_motions):
   
    # genPlot.m
    genPlot = File("genPlot.m")
    genPlot.addPFN(PFN("file://" + os.getcwd() + "/inputs/genPlot.m", "local"))
    dax.addFile(genPlot)

    xTicks = File("setXTicks.m")
    xTicks.addPFN(PFN("file://" + os.getcwd() + "/inputs/setXTicks.m", "local"))
#    dax.addFile(setXTicks)
    dax.addFile(xTicks)
#    dax.addFile(genPlot)

   
    # job
    job = Job(namespace="system", name="octave", version="1.0")

    job.uses(genPlot, link=Link.INPUT)
    job.uses(xTicks, link=Link.INPUT)

#    job.addArguments("--silent", "--eval", "genPlot;quit")
    job.addArguments("--silent","genPlot.m", "%d" % (num_materials), "%d" % (num_motions))



    # input files
    for n in range(1, num_materials + 1):
        for m in range(1, num_motions + 1):
            w = "%d.%d" %(n, m)
            node_file = bag_files["NodeDisp.out.%s" % (w)]
            job.uses(node_file, link=Link.INPUT)
            node_file = bag_files["NodeAccel.out.%s" % (w)]
            job.uses(node_file, link=Link.INPUT)
            node_file = bag_files["NodeDrift.out.%s" % (w)]
            job.uses(node_file, link=Link.INPUT)
            node_file = bag_files["NodeReaction.out.%s" % (w)]
            job.uses(node_file, link=Link.INPUT)


    f = File("Accel.jpg")
    job.uses(f, link=Link.OUTPUT)
    f = File("Disp.jpg")
    job.uses(f, link=Link.OUTPUT)
    f = File("Reaction.jpg")
    job.uses(f, link=Link.OUTPUT)
    f = File("Drift.jpg")
    job.uses(f, link=Link.OUTPUT)

    # add the file to the bag so we can look it up later
 #   bag_files[out_name] = f
            
    dax.addJob(job)

    # parent jobs
    for n in range(1, num_materials + 1):
        for m in range(1, num_motions + 1):
            w = "%d.%d" %(n, m)
            dax.depends(parent=bag_jobs["opensees_%s" % (w)], child=job)


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

add_opensees_jobs(dax, options.num_mat_props, options.num_motions)
add_plot_job(dax, options.num_mat_props, options.num_motions)

# Write the DAX to stdout
dax.writeXML(sys.stdout)



