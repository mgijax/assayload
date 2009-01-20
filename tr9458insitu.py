#!/usr/local/bin/python

#
# Program: tr9458insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate Tamplindata.txt into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr9458insitu.py
#
# Envvars:
#
# Inputs:
#
#       Tamplindata.txt, a tab-delimited file in the format:
#
#               field 1: MGI Gene Symbol
#		field 2: MGI Gene ID
#		field 3: MGI Probe ID
#
#		field 4: E7.25
#		field 5: E7.25 image
#		field 6: E7.25 specnote
#
#		field 7: E7.5
#		field 8: E7.5 image
#		field 9: E7.5 specnote
#
#		field 10: E7.75
#		field 11: E7.75 image
#		field 12: E7.75 specnote
#
#		field 13: E7.85
#		field 14: E7.85 image
#		field 15: E7.85 specnote
#
#		field 16: E7.95
#		field 17: E7.95 image
#		field 18: E7.95 specnote
#
#	E7.25key.txt
#	E7.5key.txt
#	E7.75key.txt
#	E8.5key.txt
#	E9.5key.txt
#	   tab-delimited file for each Ekey
#	
#       field 1: Ekey (equal to field 7, 10, 13, 16)
#       field 2: Strength
#       field 3: Pattern
#       field 4: Structure Print Name
#       field 5: Stage
#       field 6: Result Note
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

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

trdir = os.environ['TR_DIR']
datadir = os.environ['ASSAYLOADDATADIR']

inInSituFileName1 = os.environ['PROJECTDIR'] + '/Tamplindata.txt'

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
reference = 'J:141291'
assayType = 'RNA in situ'
createdBy = os.environ['CREATEDBY']
assayNote = ''

# constants for specimen
genotype = 'MGI:2166377'	# ICR
age = 'embryonic day'
ageNote = NULL
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Not Applicable'
specimenHybridization = 'whole mount'
specimenNote = NULL

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
    global inInSituFile
    global prepFile, assayFile, specimenFile, resultsFile
 
    try:
        inInSituFile = open(inInSituFileName1, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName1)

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

    global assayKey

    # For each line in the input file

    dictEkey = []
    specimenKey = 1
    lineNum = 0
    for line in inInSituFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# generate a directory of Ekeys
	#
        # field 7: E7.5
        # field 8: E7.5 image
        # field 9: E7.5 specnote
        #
        # field 10: E7.75
        # field 11: E7.75 image
        # field 12: E7.75 specnote
        #
        # field 13: E7.85
        # field 14: E7.85 image
        # field 15: E7.85 specnote
        #
        # field 16: E7.95
        # field 17: E7.95 image
        # field 18: E7.95 specnote
	#

	if lineNum == 1:
	    (e, age) = string.split(tokens[3], 'E')
	    x = (tokens[3], tokens[3] + 'key.txt', 3, 4, 5, age)
	    dictEkey.append(x)

	    (e, age) = string.split(tokens[6], 'E')
	    x = (tokens[6], tokens[6] + 'key.txt', 6, 7, 8, age)
	    dictEkey.append(x)

	    (e, age) = string.split(tokens[9], 'E')
	    x = (tokens[9], tokens[9] + 'key.txt', 9, 10, 11, age)
	    dictEkey.append(x)

	    (e, age) = string.split(tokens[12], 'E')
	    x = (tokens[12], tokens[12] + 'key.txt', 12, 13, 14, age)
	    dictEkey.append(x)

	    (e, age) = string.split(tokens[15], 'E')
	    x = (tokens[15], tokens[15] + 'key.txt', 15, 16, 17, age)
	    dictEkey.append(x)

	    continue


	# else process an actual data line

	#
	#  field 1: MGI Gene Symbol
	#  field 2: MGI Gene ID
	#  field 3: MGI Probe ID

	markerID = tokens[1]
	probeID = tokens[2]

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

	#
	# for each Ekey:
	#    generate the specimen row
	#

	for k in dictEkey:

	    specimenLabel = k[0]
	    figureLabel = specimenLabel
	    fpEkey = k[1]
	    specimenEkey = tokens[k[2]]
	    imageJPG = tokens[k[3]]
	    specimenNote = tokens[k[4]]
	    age = 'embryonic day ' + k[5]

	    if len(specimenEkey) == 0:
                continue

            # write the specimen (one for each Assay)

	    if specimenLabel == '':
	        continue

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

	    # write one result per tissue

	    resultKey = 1
            inResultFile = open(trdir + '/' + fpEkey, 'r')
	    lineNumResult = 0
	    for rline in inResultFile.readlines():

	        lineNumResult = lineNumResult + 1

		if lineNumResult == 1:
		    continue

                # Split the line into tokens
                rtokens = string.split(rline[:-1], TAB)

		# field 1: Ekey (first field of each E*.txt file)
		# field 2: Strength
		# field 3: Pattern
		# field 4: Structure Print Name
		# field 5: Stage
		# field 6: Result Note

		resultEkey = rtokens[0]
		strength = rtokens[1]
		pattern = rtokens[2]
		structureName = rtokens[3]
		theilerStage = rtokens[4]
		resultNote = rtokens[5]

		# insert a result for each match between the specimen and the result text file

		if resultEkey != specimenEkey:
		    continue

		resultsFile.write(str(assayKey) + TAB + \
	            	str(specimenKey) + TAB + \
	            	str(resultKey) + TAB + \
		    	strength + TAB + \
		    	pattern + TAB + \
		    	structureName + TAB + \
		    	theilerStage + TAB + \
		    	resultNote + TAB + \
			figureLabel + CRT)
		resultKey = resultKey + 1
            inResultFile.close()

	    specimenKey = specimenKey + 1

	assayKey = assayKey + 1

    # end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)

