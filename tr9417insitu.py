#!/usr/local/bin/python

#
# Program: tr9417insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate assayload.txt into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr9417insitu.py
#
# Envvars:
#
# Inputs:
#
#       assayload.txt, a tab-delimited file in the format:
#               field 1: Probe ID
#		field 2: MGI symbol
#		field 3: Marker ID
#		field 4: Specimen Label
#		field 5: Genotype ID
#		field 6: Tissue Strength
#		field 7: Figure Label	(see gxdimageload/tr9417fullsize.py/tr9417thumbnail.py)
#		field 8: Image File Name	(see gxdimageload/tr9417fullsize.py/tr9417thumbnail.py)
#		field 9: Image link	(see gxdimageload/tr9417fullsize.py/tr9417thumbnail.py)
#
#       TissueTranslation.txt, a tab-delimited file in the format:
#               field 1: Tissue Strength (from assayload, field 6)
#		field 2: Stage
#		field 3: MGI Structure
#		field 4: Strength
#		field 5: Pattern
#		field 6: Result Note
#
# Outputs:
#
#       4 tab-delimited files:
#
#	In_Situ_probeprep.txt
#	In_Situ_assay.txt
#	In_Situ_specimen.txt
#	In_Situ_results.txt
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

inInSituFile = ''	# file descriptor
inTissueFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['ASSAYLOADDATADIR']

inInSituFileName1 = datadir + '/assayload.txt'
inTissueFileName = datadir + '/TissueTranslation.txt'

prepFileName = datadir + '/In_Situ_probeprep.txt'
assayFileName = datadir + '/In_Situ_assay.txt'
specimenFileName = datadir + '/In_Situ_specimen.txt'
resultsFileName = datadir + '/In_Situ_results.txt'

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'digoxigenin'
labelCoverage = 'Not Specified'
visualizedWith = 'Alkaline phosphatase'

# constants for assay
reference = 'J:139177'
assayType = 'RNA in situ'
createdBy = os.environ['CREATEDBY']
assayNote = 'Specimen and probe information for this assay was downloaded from GenePaint.  Cryosections of fresh frozen material were fixed in 4% paraformaldehyde for 20 min. before further processing. A tyramidebiotin/streptavidin amplification step was included in the in situ hybridization procedure.'

# constants for specimen
specimenLabel1 = '1'
specimenLabel2 = '2'
genotype = 'MGI:2166653'	# NMRI
age = 'embryonic day 14.5'
ageNote = NULL
sex = 'Not Specified'
fixation = 'Fresh Frozen'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = NULL

tissueTrans = {}	# map tissue strength to rest of the line

# next available assay key
assayKey = 1

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
    global inInSituFile, inTissueFile
    global prepFile, assayFile, specimenFile, resultsFile
    global tissueTrans
 
    try:
        inInSituFile = open(inInSituFileName1, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName1)

    try:
        inTissueFile = open(inTissueFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inTissueFileName)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        specimenFile = open(specimenFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % specimenFileName)

    try:
        resultsFile = open(resultsFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % resultsFileName)

    lineNum = 0
    for line in inTissueFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)

	if lineNum == 1:
	    continue

	# else process an actual data line

	tissueStrenth = tokens[0]
	tissueTrans[tissueStrenth] = tokens

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    global assayKey

    # For each line in the input file

    specimenKey = 1
    lineNum = 0
    for line in inInSituFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)

	if lineNum == 1:
	    continue

	# else process an actual data line

	# field 1: Probe ID
	# field 2: MGI symbol
	# field 3: Marker ID
	# field 4: Specimen Label
	# field 5: Genotype ID
	# field 6: Tissue Strength
	# field 7: Figure Label	(see gxdimageload/tr9417fullsize.py/tr9417thumbnail.py)

	probeID = tokens[0]
	markerID = tokens[2]
	specimenLabel = tokens[3]
	genotypeID = tokens[4]
	tissueStrength = tokens[5]
	figureLabel = tokens[6]

	# create one assay per record

	# write the probe prep information

	prepFile.write(str(assayKey) + TAB + \
	    probeID + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assayKey) + TAB + \
	    markerID + TAB + \
	    reference + TAB + \
	    assayType + TAB + \
	    TAB + \
	    assayNote + TAB + \
	    createdBy + CRT)

        # write the specimen (one for each Assay)

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenLabel + TAB + \
	    genotypeID + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# write one result per tissue

	resultKey = 1
	allTissueStrength = string.split(tissueStrength, '|')
	for a in allTissueStrength:

		# field 1: Tissue Strength (from assayload, field 6)
		# field 2: Stage
		# field 3: MGI Structure
		# field 4: Strength
		# field 5: Pattern
		# field 6: Result Note
		# field 7: Images (jpg,jpg,...)

		t = tissueTrans[a]
		strength = t[3]
		pattern = t[4]
		tissue = t[2]
		theilerStage = t[1]
		resultNote = t[5]

		resultsFile.write(str(assayKey) + TAB + \
	            	str(specimenKey) + TAB + \
	            	str(resultKey) + TAB + \
		    	strength + TAB + \
		    	pattern + TAB + \
		    	tissue + TAB + \
		    	theilerStage + TAB + \
		    	resultNote + TAB + \
			figureLabel + CRT)
		resultKey = resultKey + 1

	assayKey = assayKey + 1

    # end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)

