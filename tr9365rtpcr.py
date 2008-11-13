#!/usr/local/bin/python

#
# Program: tr9365rtpcr.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RT-PCR input files into input files
#	for the assayload/gelload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr9365rtpcr.py
#
# Envvars:
#
# Inputs:
#
#       guo_results.txt, a tab-delimited file in the format:
#		field 1: MGI Mouse Symbol
#		field 2: MGI Marker ID
#		field 3: Probe Name
#		field 4: Metaphase II (lane 1)
#		field 5: 2-cell embryo (lane 2)
#		field 6: 4-cell embryo (lane 3)
#		field 7: 8-cell embryo (lane 4)
#		field 8: Morula (lane 5)
#		field 9: E3.5 blastocyst (lane 6)
#		field 10: E4.5 blastocyst (lane 7)
#		field 11: E4.5 inner cell mass (lane 8)
#
#	guo_lanes.txt, a tab-delimited file in the format:
#		field 1: Lane (1-8)
#		field 2: Lane Label
#		field 3: Genotype
#		field 4: Structure (printName:stage,printName:stage,...)
#		field 5: RNA
#		field 6: Age
#		field 7: Age Min/Max
#		field 8: Sex
#
#	guo_primers.txt, a tab-delimited file in the format:
#               field 1: Marker Symbol
#               field 2: MGI Marker Accession ID
#               field 3: Primer Name
#               field 4: Reference (J:#####)
#               field 5: Region Covered
#               field 6: Sequence 1
#               field 7: Sequence 2
#               field 8: Product Size
#               field 9: Notes
#               field 10: Nucleotide Sequence ID   (|-delimited)
#               field 11: Created By
#
# Outputs:
#
#       5 tab-delimited files:
#
#	RT_PCR_primer.txt
#	RT_PCR_probeprep.txt
#	RT_PCR_assay.txt
#	RT_PCR_gellane.txt
#	RT_PCR_gelband.txt
#
#       Error file
#
# Exit Codes:
#
# Assumes:
#
# Bugs:
#
# Implementation:
#

import sys
import os
import string

#globals

TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
NULL = ''

inLanesFile = ''	# file descriptor
inResultsFile = ''	# file descriptor
inPrimerFile = ''         # file descriptor
prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
gellaneFile = ''        # file descriptor
gelbandFile = ''        # file descriptor

datadir = os.environ['RTPCRDATADIR']

inLanesFileName = datadir + '/guo_lanes.txt'
inResultsFileName = datadir + '/guo_results.txt'
inPrimersFileName = datadir + '/output/newPrimer.txt'
prepFileName = datadir + '/RT_PCR_probeprep.txt'
assayFileName = datadir + '/RT_PCR_assay.txt'
gellaneFileName = datadir + '/RT_PCR_gellane.txt'
gelbandFileName = datadir + '/RT_PCR_gelband.txt'

# constants for primers
repeatUnit = NULL
moreProduct = '0'

# constants for probe prep
prepType = 'Not Specified'
hybridization = 'Not Applicable'
labelledWith = 'Not Applicable'
labelCoverage = 'Not Applicable'
visualizedWith = 'Not Specified'

# constants for assay
reference = 'J:140465'
assayType = 'RT-PCR'
assayNote = 'These results were obtained using TaqMan real time PCR. For each embryonic stage, 150 embryos were harvested and pooled together. The inner cell mass tissue was isolated using immunosurgery.'

createdBy = os.environ['CREATEDBY']

# constants for gel lanes
control = 'No'
sampleAmount = NULL
reporterGene = NULL
assayNote = NULL
ageNote = NULL
laneNote = NULL
rowNote = NULL
bandNote = NULL

# constants for gel rows
gelsize = NULL
gelunits = 'Not Specified'

mgiTypeKey = 8		# Assay
mgiPrefix = "MGI:"

# Purpose: prints error message and exits
# Returns: nothing
# Assumes: nothing
# Effects: exits with exit status
# Throws: nothing

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (string)
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    sys.exit(status)
 
# Purpose: initialize
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global inLanesFile, inResultsFile, inPrimersFile
    global primerFile, prepFile, assayFile, gellaneFile, gelbandFile
 
    try:
        inLanesFile = open(inLanesFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inLanesFileName)

    try:
        inResultsFile = open(inResultsFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResultsFileName)

    try:
        inPrimersFile = open(inPrimersFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inPrimersFileName)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        gellaneFile = open(gellaneFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % gellaneFileName)

    try:
        gelbandFile = open(gelbandFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % gelbandFileName)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    # store the primer name/id
    allPrimers = {}
    for p in inPrimersFile.readlines():
	tokens = string.split(p[:-1], TAB)
	primerName = tokens[2]
	primerID = tokens[11]
	allPrimers[primerName] = primerID

    # preps and assays

    assay = 1

    lineCount = 1
    for line in inResultsFile.readlines():

	if lineCount == 1:
	    lineCount = lineCount + 1
	    continue

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	accID = tokens[1]
	primerName = tokens[2]

	allBands = {}
	allBands[1] = tokens[3]
	allBands[2] = tokens[4]
	allBands[3] = tokens[5]
	allBands[4] = tokens[6]
	allBands[5] = tokens[7]
	allBands[6] = tokens[8]
	allBands[7] = tokens[9]
	allBands[8] = tokens[10]

	# write the probe prep information

	prepFile.write(str(assay) + TAB + \
	    allPrimers[primerName] + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assay) + TAB + \
	    accID + TAB + \
	    reference + TAB + \
	    assayType + TAB + \
	    reporterGene + TAB + \
	    assayNote + TAB + \
	    createdBy + CRT)

	lineCountB = 1

        inLanesFile = open(inLanesFileName, 'r')
        for lineB in inLanesFile.readlines():

	    if lineCountB == 1:
	        lineCountB = lineCountB + 1
	        continue

            # Split the line into tokens
            tokens = string.split(lineB[:-1], TAB)

	    laneLabel = tokens[1]
	    genotype = tokens[2]
	    structure = tokens[3]
	    rnaType = tokens[4]
	    age = tokens[5] + ' ' + tokens[6]
	    sex = tokens[7]

	    allStructures = string.split(structure, ',')

	    # write the gel lane and gel band information
    
	    lane = 1
	    row = 1

	    for s in allStructures:

	        structName, structTS = string.split(s, ':')

	        gellaneFile.write(str(assay) + TAB + \
		    str(lane) + TAB + \
		    laneLabel + TAB + \
		    genotype + TAB + \
		    rnaType + TAB + \
		    control + TAB + \
		    sampleAmount + TAB + \
		    sex + TAB + \
		    age + TAB + \
		    ageNote + TAB + \
		    laneNote + TAB + \
		    structName + TAB + \
		    structTS + CRT)

	    for s in allBands:
	        gelbandFile.write(str(assay) + TAB + \
	            str(lane) + TAB + \
	            str(row) + TAB + \
	            gelsize + TAB + \
	            gelunits + TAB + \
	            allBands[s] + TAB + \
	            rowNote + TAB + \
	            bandNote + CRT)

	    lane = lane + 1
	    row = row + 1
	    lineCountB = lineCountB + 1

	inLanesFile.close()
	assay = assay + 1
	lineCount = lineCount + 1

    #	end of "for line in inResultsFile.readlines():"

#
# Main
#

init()
process()
exit(0)

