#!/usr/bin/env python

#
# Program to read out certain metal(s) of the
# .mod model file and store to new metal.mod file
#
# 25/11/14
#

import os.path
from sys import argv,exit
from numpy import in1d
import readline
import re # RegEx

datafileloc = argv[1]
datafile, extension = os.path.splitext(datafileloc)

metal = argv[2]
metals= metal.split(',')

outfile = metal+'.mod'

infile = open(datafileloc, 'r')
outfile= open(outfile, 'w+')

header_read = False
data_read   = False
model_read  = False
links_read  = False

# Searching for metal type in line in file
def findmetal(line):
    specids = re.findall("specid=([^\r\n]+)", str(line))
    return re.findall("([A-Z][a-z|0-9]?[A-Z|0-9]+)", str(specids))

def find_remember(line):
    if line.find("cutmod_set") > 0:
        keyvals = re.findall("([a-z]+)=([0-9|.|e|E|+|-]+)", str(line))
        return dict(keyvals)
    else:
        return {}

def find_substitute(line, keyvals):
    line_keyval = re.compile("([a-z]+)=([0-9|.|a-z|A-Z|_]+)")

    line_valstrip = re.compile("([a-z]+)")

    line_values = re.findall(line_valstrip, line)

    if any(in1d(keyvals.keys(), line_values)):
        print line
        #print re.sub("redshift=([0-9|.|a-z|A-Z|_]+)
        for key in keyvals:
            line = line.replace(" 1"+key, " "+keyvals[key]+key)
            line = line.replace("(1"+key, "("+keyvals[key]+key)
            line = line.replace("=1"+key, "="+keyvals[key]+key)

    return line

#def find_substitute_readdata(line,keyvals):
#    line_keyval = re.compile("([a-z]+)=([0-9|A-Z|a-z|\(|\)|\-|.|\[|\]|,|:]+)")
#    line_interior = re.compile("\([0-9|.|-]+([A-Z|a-z|\(|\)|\-|+|e|.|\[|\]|,|:]+)\)")
#
#    line_values = re.findall(line_interior, line)
#
#    if any(in1d(keyvals.keys(), line_values)):
#        for key in keyvals:
#            line = line.replace("1"+key, keyvals[key]+key)
#
#    return line
#
#
#    #voigt   ion=28Si_IV  12.3425240     redshift=1tvala      bturb=1bvala         1ta          specid=keckSiIV,vltSiIV,keckSiIVb,vltSiIVb

# Set of parameters defined in the input model file
# as # cmod_set param=value
custom_set = {}
# Printing lines if they impose restrictions on the
# parameters that are stored in custom_set
conditions = ['lim', 'emission', 'absorption']

for line in infile:

    lineargs = re.findall("([a-z]+)", line)

    if "data read" in line:
        header_read = True
        outfile.write(line)
    if "data end" in line:
        data_read   = True
        outfile.write(line)
        outfile.write("model read\n")
    if "model end" in line:
        model_read  = True
        outfile.write(line)
    if "links end" in line:
        links_read  = True
        outfile.write(line)

    if not header_read:
        outfile.write(line)

    elif (header_read & (not model_read)):
        # Read in cmod_set 
        # so values are stored
        custom_set.update(find_remember(line))
        if any(in1d(metals, findmetal(line))):
            # If we have any values on the form:
            #   1ta
            #   (1ta
            #   =1ta
            # which corresponds to a variable that has been set
            # and hopefully, it was defined under
            # model read
            # # cutmod_set ta=1.5e+4
            # this value will be substituted into the output file as
            #   1.5e+4ta
            #   (1.5e+4ta
            #   =1.5e+4ta
            line = find_substitute(line, custom_set)
            print find_substitute(line, custom_set)
            outfile.write(line)

    if ((not model_read) & any(in1d(conditions, lineargs))):
        outfile.write(line)



outfile.close()
infile.close()

print "Done, see newly created file: %s" % (str(outfile))
