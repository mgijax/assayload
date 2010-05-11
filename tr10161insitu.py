#!/usr/local/bin/python

#
# Program: tr10161insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 10161/Assay/Assays.txt into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10161insitu.py
#
# Envvars:
#
# Inputs:
#
#       Assays.txt, a tab-delimited file in the format:
#
#               field 1: Probe ID
#               field 2: MGI Marker Accession ID
#               field 3: Reference (J:#####)
#               field 4: Assay Type
#               field 5: Reporter Gene
#               field 6: Assay Note
#               field 7: Created By
#
#       Specimen file, a tab-delimited file in the format:
#               field 1: MGI Marker Accession ID
#               field 2: Specimen Label
#               field 3: Genotype ID
#               field 4: Age
#               field 5: Age Note
#               field 6: Sex
#               field 7: Fixation
#               field 8: Embedding Method
#               field 9: Hybridization
#               field 10: Specimen Note
#
#       Specimen Results file, a tab-delimited file in the format:
#               field 1: Assay #
#               field 2: Specimen #
#               field 3: Result #
#               field 4: Strength
#               field 5: Pattern
#               field 6: MGI Structure Name
#               field 7: MGI Structure Theiler Stage
#               field 8: Result Note
#               field 9: Comma-Separated Image Names (if any)
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

trdir = os.environ['TR_DIR']
datadir = os.environ['ASSAYLOADDATADIR']

inAssayFileName = os.environ['PROJECTDIR'] + '/Assays/Assays.txt'
inSpecimenFileName = os.environ['PROJECTDIR'] + '/Assays/Specimens.txt'
inResultFileName = os.environ['PROJECTDIR'] + '/Assays/Results.txt'

prepFileName = datadir + '/In_Situ_probeprep.txt'
assayFileName = datadir + '/In_Situ_assay.txt'
specimenFileName = datadir + '/In_Situ_specimen.txt'
resultsFileName = datadir + '/In_Situ_results.txt'

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'digoxigenin'
visualizedWith = 'Alkaline phosphatase'

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
    global inAssayFile, inSpecimenFile, inResultFile
    global prepFile, assayFile, specimenFile, resultsFile
 
    try:
        inAssayFile = open(inAssayFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssayFileName)

    try:
        inSpecimenFile = open(inSpecimenFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inSpecimenFileName)

    try:
        inResultFile = open(inResultFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResultFileName)

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

    assayKey = 1
    assayLookup = {}

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	probeID = tokens[0]
	markerID = tokens[1]
	reference = tokens[2]
	assayType = tokens[3]
	reporter = tokens[4]
	assayNote = tokens[5]
	createdBy = tokens[6]

	# create one assay per record

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

	if not assayLookup.has_key(markerID):
	    assayLookup[markerID] = []
        assayLookup[markerID] = assayKey

	assayKey = assayKey + 1

    inAssayFile.close()

    # end of "for line in inAssayFile.readlines():"

    # write one specimen per assay

    specimenKey = 1
    prevMarkerID = ''
    specimenLookup = {}

    for sline in inSpecimenFile.readlines():

        # Split the line into tokens
        tokens = string.split(sline[:-1], TAB)

	markerID = tokens[0]
	specimenID = tokens[1]
	genotype = tokens[2]
	age = tokens[3]
	ageNote = tokens[4]
	sex = tokens[5]
	fixation = tokens[6]
	embedding = tokens[7]
	specimenHybridization = tokens[8]

	if markerID != prevMarkerID:
            specimenKey = 1
	    assayKey = assayLookup[markerID]
	    prevMarkerID = markerID

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenID + TAB + \
	    genotype + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    CRT)

	if not specimenLookup.has_key(specimenID):
	    specimenLookup[specimenID] = []
        specimenLookup[specimenID] = specimenKey

        specimenKey = specimenKey + 1

    inSpecimenFile.close()

    # write one results per specimen

    resultKey = 1
    prevAssayID = ''
    prevSpecimenID = ''

    for rline in inResultFile.readlines():

        # Split the line into tokens
        tokens = string.split(rline[:-1], TAB)

	markerID = tokens[0]
	specimenID = tokens[1]
	strength = tokens[2]
	pattern = tokens[3]
	structureName = tokens[4]
	structureTheilerStage = tokens[5]
	resultNote = tokens[6]
	imageNames = tokens[7]

	if markerID != prevMarkerID:
            specimenKey = 1
	    assayKey = assayLookup[markerID]
	    prevMarkerID = markerID

	if specimenID != prevSpecimenID:
            resultKey = 1
	    specimenKey = specimenLookup[specimenID]
	    prevSpecimenID = specimenID

	resultsFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    str(resultKey) + TAB + \
	    strength + TAB + \
	    pattern + TAB + \
	    structureName + TAB + \
	    structureTheilerStage + TAB + \
	    resultNote + TAB + \
	    imageNames + CRT)

        resultKey = resultKey + 1

    inResultFile.close()

    # end of "for line in inResultFile.readlines():"

#
# Main
#

init()
process()
exit(0)

