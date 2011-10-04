#!/usr/local/bin/python

#
# Program: tr10840rtpcr1.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10840/SuraniTableS7.txt
#	into input files for the gelload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10840rtpcr1.py
#
# Envvars:
#
#	LOADFILE1
#	PROBEPREP_FILE
#	ASSAY_FILE
#	LANE_FILE
#	BAND_FILE
#	REFERENCE
#	CREATEDBY
#
# Inputs:
#
#       SuraniTableS7.txt (LOADFILE1), a tab-delimited file in the format:
#
#               field 0: Marker ID
#               field 1: Probe Name
#               field 2: Lane Label
#               field 3: Genotype
#	        field 4: RNAType
#	        field 5: Gel Control
#		field 6: Age
#		field 7: Lane Note
#		field 8: Structure Name
#		field 9: Theiler Stage
#		field 10 Band Strength
#		field 11: Band Note
#
# Outputs:
#
#       4 tab-delimited files:
#
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
import db

#globals

TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
NULL = ''

inAssayFile = ''	# file descriptor
inLaneFile = ''	# file descriptor
inBandFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
laneFile = ''       # file descriptor
bandFile = ''        # file descriptor

inAssay = os.environ['LOADFILE1']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
laneFileName = os.environ['LANE_FILE']
bandFileName = os.environ['BAND_FILE']

assayNote = ''
reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# constants for probe prep
assayType = 'RT-PCR'
prepType = 'Not Specified'
hybridization = 'Not Specified'
labelledWith = 'Not Applicable'
visualizedWith = 'Not Applicable'

# gel lane
sampleAmount = ''
sex = 'Not Specified'
ageNote = ''

# gel row/band
bandSize = ''
bandUnits = 'Not Specified'
rowNote = ''

# structure lookup
structureLookup = {}

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
    global inAssayFile
    global prepFile, assayFile, laneFile, bandFile
 
    try:
        inAssayFile = open(inAssay, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssay)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        laneFile = open(laneFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % laneFileName)

    try:
        bandFile = open(bandFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % bandFileName)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0
    laneKey = 0
    bandKey = 1
    prevMarker = ''
    prevProbe = ''

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	markerID = tokens[0]
	probeID = tokens[1]
	laneLabel = tokens[2]
	genotype = tokens[3]
	rnaType = tokens[4]
	gelControl = tokens[5]
	age = tokens[6]
	laneNote = tokens[7]
	structureName = tokens[8]
	structureTheilerStage = tokens[9]
	bandStrength = tokens[10]
	bandNote = tokens[11]

	if markerID != prevMarker:
	    prevMarker = markerID
	    prevProbe = ''

	if probeID != prevProbe:

	    assayKey = assayKey + 1
	    prevProbe = probeID
            laneKey = 0

	    prepFile.write(str(assayKey) + TAB + \
	        probeID + TAB + \
	        prepType + TAB + \
	        hybridization + TAB + \
	        labelledWith + TAB + \
	        visualizedWith + CRT)

	    # write the assay information

            assayFile.write(str(assayKey) + TAB + \
                markerID + TAB + \
                reference + TAB + \
                assayType + TAB + \
                TAB + \
                assayNote + TAB + \
                createdBy + CRT)

        laneKey = laneKey + 1
	laneFile.write(str(assayKey) + TAB + \
	    str(laneKey) + TAB + \
	    laneLabel + TAB + \
	    genotype + TAB + \
	    rnaType + TAB + \
	    gelControl + TAB + \
	    sampleAmount + TAB + \
	    sex + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    laneNote + TAB + \
	    structureName + TAB + \
	    structureTheilerStage + CRT)

	bandKey = 1
	bandFile.write(str(assayKey) + TAB + \
	    str(laneKey) + TAB + \
	    str(bandKey) + TAB + \
	    bandSize + TAB + \
	    bandUnits + TAB + \
	    bandStrength + TAB + \
	    rowNote + TAB +
	    bandNote + CRT)

    inAssayFile.close()
    assayFile.close()
    laneFile.close()
    bandFile.close()

    # end of "for line in inAssayFile.readlines():"

# Main
#

init()
process()
exit(0)

