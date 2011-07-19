#!/usr/local/bin/python

#
# Program: tr10629insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10629/LoaderDocuments/LoadFile1.txt 
#		10629/LoaderDocuments/LoadFile2.txt 
#		10629/LoaderDocuments/LoadFile3.txt 
#	into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10629insitu.py
#
# Envvars:
#
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
#               field 1: Labled with
#               field 2: Visualized with
#               field 3: Gene Symbol
#               field 4: MGI Marker Accession ID
#
#       LoadFile2.txt, a tab-delimited file in the format:
#
#               field 0: Probe ID
#               field 1: Specimen Label
#               field 2: Genotype ID
#               field 3: Age
#               field 4: Age note
#               field 5: Sex
#               field 6: Fixation
#               field 7: Embedding
#               field 8: Hybridization
#               field 9: Specimen Note
#
#       LoadFile3.txt, a tab-delimited file in the format:
#
#               field 0: Specimen Label
#		field 1: EMAP_ID Structure_Name
#		field 2: Strength
#		field 3: Pattern
#		field 4: Result Note
#		field 5: Image Names
#
#	StructureLookup.txt, a tab-delimited file in the format:
#
#		field 0: EMAP id
#		field 1: GUDMAP name
#		field 2: MGI ID of the structure
#		field 3: Theiler Stage
#		field 4: GXD Structure Print Name
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
inStructureFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor
structureFile = ''      # file descriptor

inAssayName = os.environ['LOADFILE1']
inSpecimenFileName = os.environ['LOADFILE2']
inResultFileName = os.environ['LOADFILE3']
inStructureFileName = os.environ['STRUCTURE_LOOKUP']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
specimenFileName = os.environ['SPECIMEN_FILE']
resultsFileName = os.environ['RESULTS_FILE']

# constants for probe prep
assayType = 'RNA in situ'
prepType = 'RNA'
hybridization = 'Antisense'
assayNote = ''
reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# structure lookup
imageLookup = {}

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
    global inAssayFile, inSpecimenFile, inResultFile, inImageFile
    global prepFile, assayFile, specimenFile, resultsFile
    global imageLookup
 
    try:
        inAssayFile = open(inAssayName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssayName)

    try:
        inSpecimenFile = open(inSpecimenFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inSpecimenFileName)

    try:
        inResultFile = open(inResultFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResultFileName)

    try:
        inStructureFile = open(inStructureFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inStructureFileName)

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

    for line in inStructureFile.readlines():
	tokens = string.split(line[:-1], TAB)
	structureID = int(tokens[0])
	imageLookup[structureID] = []
	imageLookup[structureID].append(tokens)
    inStructureFile.close()

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
	labelledWith = tokens[1]
	visualizedWith = tokens[2]
	markerID = tokens[4]

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

	if not assayLookup.has_key(probeID):
	    assayLookup[probeID] = []
        assayLookup[probeID] = assayKey

	assayKey = assayKey + 1

    inAssayFile.close()

    # end of "for line in inAssayFile.readlines():"

    # write one specimen per assay

    specimenKey = 1
    assay2Lookup = {}
    assay3Lookup = {}

    for sline in inSpecimenFile.readlines():

        # Split the line into tokens
        tokens = string.split(sline[:-1], TAB)

	probeID = tokens[0]
	specimenID = tokens[1]
	genotype = tokens[2]
	age = tokens[3]
	ageNote = tokens[4]
	sex = tokens[5]
	fixation = tokens[6]
	embedding = tokens[7]
	specimenHybridization = tokens[8]
	specimenNote = tokens[9]

	assayKey = assayLookup[probeID]

	# assay by specimen id
	if not assay3Lookup.has_key(specimenID):
	    assay3Lookup[specimenID] = []
	    specimenKey = 1
	else:
	    specimenKey = assay3Lookup[specimenID] + 1

        assay3Lookup[specimenID] = specimenKey

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
	    specimenNote + CRT)

	# assay by specimen id
	if not assay2Lookup.has_key(specimenID):
	    assay2Lookup[specimenID] = []
        assay2Lookup[specimenID] = assayKey

    inSpecimenFile.close()

    # write one results per specimen

    resultLookup = {}

    for rline in inResultFile.readlines():

        # Split the line into tokens
        tokens = string.split(rline[:-1], TAB)

	specimenID = tokens[0]
	edinID = int(tokens[1])
	strength = tokens[2]
	pattern = tokens[3]
	resultNote = tokens[4]
	imageName = tokens[5]

	r = imageLookup[edinID]
	structureName = r[0][4]
	structureTheilerStage = r[0][3]

	assayKey = assay2Lookup[specimenID]
	specimentKey = assay3Lookup[specimenID]

	if not resultLookup.has_key(specimenID):
	    resultLookup[specimenID] = []
	    resultKey = 1
	else:
	    resultKey = resultLookup[specimenID] + 1

        resultLookup[specimenID] = resultKey

	resultsFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    str(resultKey) + TAB + \
	    strength + TAB + \
	    pattern + TAB + \
	    structureName + TAB + \
	    str(structureTheilerStage) + TAB + \
	    resultNote + TAB + \
	    imageName + CRT)

    inResultFile.close()

    # end of "for line in inResultFile.readlines():"
#
# Main
#

init()
process()
exit(0)

