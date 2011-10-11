#!/usr/local/bin/python

#
# Program: tr10872insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10872/MakeAssays.txt
#	into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10872insitu.py
#
# Envvars:
#
#	LOADFILE1
#	LOADFILE2
#       PROBEPREP_FILE
#       ASSAY_FILE
#       SPECIMEN_FILE
#       RESULTS_FILE
#
# Inputs:
#
#       Tamplindata.txt (LOADFILE1), a tab-delimited file in the format:
#
#               field 0: MGI Marker ID
#               field 1: Probe ID
#
#               field 2: E7.0
#               field 3: E7.0 image
#
#               field 4: E7.5
#               field 5: E7.5 image
#
#               field 6: E7.75
#               field 7: E7.75 image
#
#               field 8: E8.0
#               field 9: E8.0 image
#
#               field 10: E8.5
#               field 11: E8.5 image
#
#               field 12: E9.0
#               field 13: E9.0 image
#
#               field 14: E9.5
#               field 15: E9.5 image
#
#       allE.txt (LOADFILE2), a tab-delimited file in the format:
#
#               field 0: E ID (E7.0, E7.5, E7.75, E8.0, E8.5, E9.0, E9.5)
#		field 1: Strenght (their version)
#		field 2: Strenght (our version)
#		field 3: Pattern
#		field 4: Structure Name
#		field 5: Structure Theiler State
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

inAssayFile = ''	# file descriptor
inResultFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

inAssay = os.environ['LOADFILE1']
inResult = os.environ['LOADFILE2']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
specimenFileName = os.environ['SPECIMEN_FILE']
resultsFileName = os.environ['RESULTS_FILE']

assayType = 'RNA in situ'

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'digoxigenin'
visualizedWith = 'Alkaline phosphatase'
assayNote = ''
reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# specimen
genotype = 'MGI:2166377'
ageNote = ''
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Not Applicable'
specimenHybridization = 'whole mount'
specimenNote = 'These data are reported in Table S2.'

# structure lookup
structureLookup = {}

# results lookup
resultLookup = {}

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
    global inAssayFile, inResultFile
    global prepFile, assayFile, specimenFile, resultsFile
    global resultLookup
 
    try:
        inAssayFile = open(inAssay, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssay)

    try:
        inResultFile = open(inResult, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResult)

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

    for line in inResultFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	specimenLabel = tokens[0]
	value = tokens[1]
	strength = tokens[2]
	pattern = tokens[3]
	structureName = tokens[4]
	structureTheilerStage = tokens[5]
	resultNote = tokens[6]

	resultID = specimenLabel + ':' + value
	if not resultLookup.has_key(resultID):
	    resultLookup[resultID] = []
	resultLookup[resultID].append([strength, pattern, structureName, structureTheilerStage, resultNote])

    inResultFile.close()
    #for r in resultLookup:
#	print r[0], r[1], r[2], r[3], r[4]

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	markerID = tokens[0]
	probeID = tokens[1]

	# write the probe prep information

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

	#
	# specimens
	#
	# 2 = E7.0
	# 4 = E7.5
	# 6 = E7.75
	# 8 = 8.0
	# 10 = 8.5
	# 12 = 9.0
	# 14 = 9.5
	#

	specimenKey = 1

	for s in range(2,15,2):

	    # one specimen per probe/marker

	    # the value + image specifies the result it is attached to
	    value = tokens[s]
	    image = tokens[s+1]
	    image = string.replace(image, '|', ',')

	    if len(value) == 0:
		continue

	    age = 'embryonic day '

	    if s == 2:
		specimenLabel = 'E7.0'
		age = age + '7.0'
	    elif s == 4:
		specimenLabel = 'E7.5'
		age = age + '7.5'
	    elif s == 6:
		specimenLabel = 'E7.75'
		age = age + '7.75'
	    elif s == 8:
		specimenLabel = 'E8.0'
		age = age + '8.0'
	    elif s == 10:
		specimenLabel = 'E8.5'
		age = age + '8.5'
	    elif s == 12:
		specimenLabel = 'E9.0'
		age = age + '9.0'
	    elif s == 14:
		specimenLabel = 'E9.5'
		age = age + '9.5'

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

	    #
	    # find result data for specimenLabel/value combination
	    #
	    resultKey = 1
	    resultID = specimenLabel + ':' + value
	    #print specimenLabel, value, image
	    for r in resultLookup[resultID]:
	        #print r[0], r[1], r[2], r[3], r[4]
		strength = r[0]
		pattern = r[1]
		structureName = r[2]
		structureTheilerStage = r[3]
		resultNote = r[4]
	        resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
	            strength + TAB + \
	            pattern + TAB + \
	            structureName + TAB + \
	            structureTheilerStage + TAB + \
	            resultNote + TAB + \
	            image + CRT)

	        resultKey = resultKey + 1

            specimenKey = specimenKey + 1

	# one result per specimen

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

