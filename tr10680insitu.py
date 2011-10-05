#!/usr/local/bin/python

#
# Program: tr10680insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10680/MakeAssays.txt
#	into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10680insitu.py
#
# Envvars:
#
#	LOADFILE1
#	PROBEPREP_FILE
#	ASSAY_FILE
#	SPECIMEN_FILE
#	RESULTS_FILE
#	REFERENCE
#	CREATEDBY
#
# Inputs:
#
#       LoadFile1.txt, a tab-delimited file in the format:
#
#               field 0: Probe ID
#               field 1: Gene ID
#               field 2: Specimen Label
#               field 3: Age
#	        field 4: Strength        
#		field 5: Pattern 
#		field 6: Structure Name   
#		field 8: TS
#		field 9: ResultNote
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

inAssayFile = ''	# file descriptor
inSpecimenFile = ''	# file descriptor
inResultFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

inAssay = os.environ['LOADFILE1']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
specimenFileName = os.environ['SPECIMEN_FILE']
resultsFileName = os.environ['RESULTS_FILE']

# constants for probe prep
assayType = 'RNA in situ'
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'digoxigenin'
visualizedWith = 'Alkaline phosphatase'
assayNote = ''
reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# specimen
genotype = 'MGI:2166522'
ageNote = ''
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = ''

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
    global prepFile, assayFile, specimenFile, resultsFile
 
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
        specimenFile = open(specimenFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % specimenFileName)

    try:
        resultsFile = open(resultsFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % resultsFileName)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0
    assayLookup = []

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	probeID = tokens[0]
	markerID = tokens[1]
	specimenLabel = tokens[2]
	age = tokens[3]
	strength = tokens[4]
	pattern = tokens[5]
	structureName = tokens[6]
	structureTheilerStage = tokens[7]
	resultNote = tokens[8]

	# create one assay per record

	# write the probe prep information

	if probeID not in assayLookup:

	    assayKey = assayKey + 1

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

	    assayLookup.append(probeID)
            specimenKey = 1

	# one specimen per probe/marker

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenLabel + TAB + \
	    genotype + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# one result per specimen

	resultKey = 1

	resultsFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    str(resultKey) + TAB + \
	    strength + TAB + \
	    pattern + TAB + \
	    structureName + TAB + \
	    structureTheilerStage + TAB + \
	    resultNote + TAB + \
	    CRT)

        specimenKey = specimenKey + 1

    inAssayFile.close()
    assayFile.close()
    specimenFile.close()
    resultsFile.close()

    # end of "for line in inAssayFile.readlines():"

# Main
#

init()
process()
exit(0)

