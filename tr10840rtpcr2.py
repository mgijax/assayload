#!/usr/local/bin/python

#
# Program: tr10840rtpcr2.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10840/SuraniTableS6Bands.txt
#		10840/SuraniTableS6Lanes.txt
#	into input files for the gelload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10840rtpcr2.py
#
# Envvars:
#
#	LOADFILE2
#	LOADFILE3
#	PRIMERNAME
#	PROBEPREP_FILE
#	ASSAY_FILE
#	LANE_FILE
#	BAND_FILE
#	REFERENCE
#	CREATEDBY
#
# Inputs:
#
#       SuraniTableS6Bands.txt (LOADFILE2), a tab-delimited file in the format:
#
#               field 0: Marker ID
#               field 1: Probe Name
#               field 2: GelLane/BandStrength
#               field 3: GelLane/BandNote
#
#	row 0:
#		field 2: Gel Lane + ' BandStength'
#		field 3: Gel Lane + ' BandNote'
#
#	row every 2:
#		field 2: BandStrength
#		field 3: BandNote
#
#		field 4: BandStrength
#		field 5: BandNote
#
#		field 6: BandStrength
#		field 7: BandNote
#		...
#
#	example:
#		
#		each assay as 36 gel lanes
#		each gel lane has 1 band
#
#               Marker       Probe          Gel Lane
#		                            Band Strength  Band Note
#
#		MGI:1100518  Mm00484742_m1  2-cell-A1-1
#                                           Absent         The Ct value of this reaction was 40.
#
#		MGI:1100518  Mm00484742_m1  2-cell-A1-2
#		                            Absent         The Ct value of this reaction was 40.
#
#		MGI:1100518  Mm00484742_m1  2-cell-A2-1
#                                           Absent         The Ct value of this reaction was 40.
#
#		MGI:1100518  Mm00484742_m1  2-cell-A2-2
#                                           Absent         The Ct value of this reaction was 40.
#
#		MGI:105081   Mm03053485_s1  2-cell-A1-1
#                                           Present        The Ct value of this reaction was 31.38.
#
#		MGI:105081   Mm03053485_s1  2-cell-A1-2
#		                            Present        The Ct value of this reaction was 31.38.
#
#		MGI:105081   Mm03053485_s1  2-cell-A2-1
#                                           Present        The Ct value of this reaction was 31.38.
#
#		MGI:105081   Mm03053485_s1  2-cell-A2-2
#                                           Present        The Ct value of this reaction was 31.38.
#
#       SuraniTableS6Lanes.txt (LOADFILE3), a tab-delimited file in the format:
#
#               field 0: GelLane # (means nothing)
#               field 1: GelLane ==> SuraniTableS6Band.txt:row 0:field 2,4,6,...
#               field 2: Age
#               field 3: LaneNote
#               field 4: StructureName
#               field 5: TheilerStage
#
#       newPrimer.txt (PRIMERNAME)
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
inLaneFile = ''		# file descriptor
inPrimerFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
laneFile = ''           # file descriptor
bandFile = ''           # file descriptor

inAssay = os.environ['LOADFILE2']
inLane = os.environ['LOADFILE3']
inPrimer = os.environ['PRIMERNAME']

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
hybridization = 'Not Applicable'
labelledWith = 'Not Applicable'
visualizedWith = 'Not Specified'

# gel lane
genotype = 'MGI:2169515'
rnaType = 'total'
gelControl = 'No'
sampleAmount = ''
sex = 'Not Specified'
ageNote = ''

# gel row/band
bandSize = ''
bandUnits = 'Not Specified'
rowNote = ''

# primer lookup
primerLookup = {}

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
    global inAssayFile, inLaneFile, inPrimerFile
    global prepFile, assayFile, laneFile, bandFile
    global primerLookup
 
    try:
        inAssayFile = open(inAssay, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssay)

    try:
        inLaneFile = open(inLane, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inLane)

    try:
        inPrimerFile = open(inPrimer, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inPrimer)

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

    for line in inPrimerFile.readlines():
        tokens = string.split(line[:-1], TAB)
	primerName = tokens[2]
	primerID = tokens[11]
	primerLookup[primerName] = []
	primerLookup[primerName].append(primerID)
    inPrimerFile.close()

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0

    # bandFile

    row = 0
    for line in inAssayFile.readlines():

	row = row + 1

        # Split the line into tokens

	# 36 instances * 2 = 72 + 2 (marker, probe)
	# geltokens holds the gelLane "id"
	if row == 1:
            geltokens = string.split(line[:-1], TAB)
	    for s in range(2,74,2):
	         geltokens[s] = string.replace(geltokens[s], '  BandStrength', '')
	    #print geltokens
	else:
            tokens = string.split(line[:-1], TAB)

	    markerID = tokens[0]
	    primerName = tokens[1]

	    assayKey = assayKey + 1

	    primerID = primerLookup[primerName][0]

	    prepFile.write(str(assayKey) + TAB + \
	        primerID + TAB + \
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

	    #
	    # band stuff
	    #

	    laneKey = 0
	    for s in range(2,74,2):
	        laneKey = laneKey + 1
	        bandKey = 1
		bandStrength = tokens[s]
		bandNote = tokens[s+1]

		#print markerID, geltokens[s], bandStrength

	        bandFile.write(str(assayKey) + TAB + \
	            str(laneKey) + TAB + \
	            str(bandKey) + TAB + \
	            bandSize + TAB + \
	            bandUnits + TAB + \
	            bandStrength + TAB + \
	            rowNote + TAB +
	            bandNote + CRT)

    # end of bandFile

    #
    # laneFile
    #

    for assayKey in range(1,37):

	laneKey = 0

        for line in inLaneFile.readlines():

            tokens = string.split(line[:-1], TAB)

	    gelLane = tokens[1]
	    age = tokens[2]
	    laneNote = tokens[3]
	    structureName = tokens[4]
	    structureTheilerStage = tokens[5]

            laneKey = laneKey + 1

	    laneFile.write(str(assayKey) + TAB + \
	        str(laneKey) + TAB + \
	        gelLane + TAB + \
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

	# reset file pointer to the beginning
	inLaneFile.seek(0)

    inAssayFile.close()
    inLaneFile.close()

    assayFile.close()
    bandFile.close()
    laneFile.close()

    # end of "for line in inAssayFile.readlines():"

# Main
#

init()
process()
exit(0)

