#!/usr/local/bin/python

#
# Program: tr10976insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10976/Jackie_Load_Documents/GP2_loadfile.txt
#	into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10976insitu.py
#
# Envvars:
#
#	LOADFILE1
#       PROBEPREP_FILE
#       ASSAY_FILE
#       SPECIMEN_FILE
#       RESULTS_FILE
#
# Inputs:
#
#       GP2_loadfile.txt (LOADFILE1), a tab-delimited file in the format:
#	
#	field 1:  Assay #
#	field 2:  Specimen #
#	field 3:  Probe name
#	field 4:  Probe MGI ID
#	field 5:  Gene symbol
#	field 6:  Gene MGI ID
#	field 7:  Analysis ID
#	field 8:  Specimen Label
#	field 9:  Strain MGI ID
#	field 10: Strength
#	field 11: Pattern
#	field 12: Best Image
#	field 13: Image jpg
#	field 14: Figure Label
#	field 15: Note
#	
#	StructureLookupTable.txt (LOADFILE2):
#	field 1: Structure key (GP)
#	field 2: Structure name
#	field 3: Structture key (MGD)
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
inStructureFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

inAssay = os.environ['LOADFILE1']
inStructure = os.environ['LOADFILE2']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
specimenFileName = os.environ['SPECIMEN_FILE']
resultsFileName = os.environ['RESULTS_FILE']

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'digoxigenin'
visualizedWith = 'Alkaline phosphatase'

assayType = 'RNA in situ'
assayNote = 'Cryosections of fresh frozen material were fixed in 4% paraformaldehyde for 20 min. before further processing. A tyramide-biotin/streptavidin amplification step was included in the in situ hybridization procedure.  These data were downloaded from GenePaint.org on 1/30/12.'

reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# specimen
age = 'embryonic day 14.5'
ageNote = ''
sex = 'Not Specified'
fixation = 'Fresh Frozen'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = ''

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
    global inAssayFile, inStructureFile
    global prepFile, assayFile, specimenFile, resultsFile
    global structureLookup
 
    try:
        inAssayFile = open(inAssay, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssay)

    try:
        inStructureFile = open(inStructure, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inStructure)

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

    #
    # stringing with field 10, every 6th field is the gbKey
    #
    key = 10
    for line in inStructureFile.readlines():
        tokens = string.split(line[:-1], TAB)
	mgdKey = tokens[2]
	if not structureLookup.has_key(key):
		structureLookup[key] = []
	structureLookup[key].append(mgdKey)
	key = key + 6
    inStructureFile.close()
    #print structureLookup

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    # For each line in the input file

    assayLookup = []

    lineNum = 0
    for line in inAssayFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	if lineNum <= 2:
	    continue

	assayKey = tokens[0]
	specimenKey = tokens[1]
	probeID = tokens[3]
	markerID = tokens[5]

	specimenLabel = tokens[7]
	genotype = tokens[8]

	if assayKey not in assayLookup:

	     # write the probe prep information

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
	# iterate thru 10-15 X 102 times
	#
	resultsKey = 1
	for s in range(9,621,6):

#	field 10: Strength
#	field 11: Pattern
#	field 12: Best Image
#	field 13: Image jpg
#	field 14: Figure Label
#	field 15: Note

		strength = tokens[s]
		pattern = tokens[s+1]
		bestimage = tokens[s+2]
		imagejpg = tokens[s+3]
		figureLabel = tokens[s+4]
		resultNote = tokens[s+5]
		structureName = structureLookup[s+1][0]

		if strength == 'weak':
			strength = 'Weak'
		elif strength == 'medium':
			strength = 'Moderate'
		elif strength == 'strong':
			strength = 'Strong'
		elif strength == 'none':
			strength = 'Absent'

		if len(strength) == 0:
			continue

		if pattern == 'ubiquitous':
			pattern = 'Ubiquitous'
		elif pattern == 'regional':
			pattern = 'Regionally restricted'
		elif pattern == 'scattered':
			pattern = 'Scattered'
		elif pattern == 'none':
			pattern = 'Not Specified'
		else:
			pattern = 'Not Applicable'

        	resultsFile.write(str(assayKey) + TAB + \
                     	str(specimenKey) + TAB + \
		     	str(resultsKey) + TAB + \
			strength + TAB + \
			pattern + TAB + \
			structureName + TAB + \
			'22' + TAB + \
			resultNote + TAB + \
			imagejpg + CRT)

		resultsKey = resultsKey + 1

	assayLookup.append(assayKey)

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

