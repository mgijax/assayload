#!/usr/local/bin/python

#
# Program: tr11471immuno.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 11471/dataload/Immunohistochemistry 
#	files into input files for the immunoload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr11471immuno.py
#
# Envvars:
#
# Inputs:
#
#	StructureLookUp.txt, a tab-delimited file in the format:
#		field 1: EMAP ID
#		field 2: Structure key
#
#       AntibodyPrep.txt, a tab-delimited file in the format:
#               field 1: Anitbody ID
#               field 2: Secondary Antibody
#               field 3: Labelled With
#
#       ImmunohistochemistryAssay.txt, a tab-delimited file in the format:
#
#               field 1: Anitbody ID
#               field 2: MGI Marker Accession ID
#               field 3: Reference (J:#####)
#               field 4: Assay Type
#               field 5: Reporter Gene
#               field 6: Assay Note
#               field 7: Created By
#
#       ImmunohistochemistrySpecimen.txt, a tab-delimited file in the format:
#               field 1: Anitbody ID
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
#       ImmunohistochemistryResults.txt, a tab-delimited file in the format:
#               field 1: Specimen Label
#               field 2: Strength
#               field 3: Pattern
#		field 4: EMAP id
#               field 5: Result Note
#               field 6: Image
#
# Outputs:
#
#       4 tab-delimited files:
#
#	Immuno_antibodyprep.txt
#	Immuno_assay.txt
#	Immuno_specimen.txt
#	Immuno_results.txt
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

inPrepFile = ''		# file descriptor
inAssayFile = ''	# file descriptor
inSpecimenFile = ''	# file descriptor
inResultFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

trdir = os.environ['TR_DIR']
datadir = os.environ['ASSAYLOADDATADIR']

inStructureFileName = os.environ['PROJECTDIR'] + '/dataload/StructureLookUp.txt'
inPrepFileName = os.environ['PROJECTDIR'] + '/dataload/Immunohistochemistry/AntibodyPrep.txt'
inAssayFileName = os.environ['PROJECTDIR'] + '/dataload/Immunohistochemistry/ImmunoAssay.txt'
inSpecimenFileName = os.environ['PROJECTDIR'] + '/dataload/Immunohistochemistry/ImmunoSpecimen.txt'
inResultFileName = os.environ['PROJECTDIR'] + '/dataload/Immunohistochemistry/ImmunoResults.txt'

prepFileName = datadir + '/Immuno_prep.txt'
assayFileName = datadir + '/Immuno_assay.txt'
specimenFileName = datadir + '/Immuno_specimen.txt'
resultsFileName = datadir + '/Immuno_results.txt'

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
    global inStructureFile, inPrepFile, inAssayFile, inSpecimenFile, inResultFile
    global prepFile, assayFile, specimenFile, resultsFile
 
    try:
        inStructureFile = open(inStructureFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inStructureFileName)

    try:
        inPrepFile = open(inPrepFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inPrepFileName)

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
    prepLookup = {}
    assayLookup = {}

    # emap ID -> structure name, structure theiler stage
    emapLookup = {}
    for line in inStructureFile.readlines():
        tokens = string.split(line[:-1], TAB)
	results = db.sql('select s.printName as structure, t.stage from GXD_Structure s, GXD_TheilerStage t where s._Structure_key = %s and s._Stage_key = t._Stage_key' % tokens[1], 'auto')
	key = tokens[0]
	value = results[0]
	emapLookup[key] = []
	emapLookup[key].append(value)
    inStructureFile.close()
    #print emapLookup

    # create lookup of images (_image_key, figureLabel)

    imageLookup = {}
    results = db.sql('''
	select ii._ImagePane_key, i.figureLabel 
	from IMG_Image i, IMG_ImagePane ii
	where i.figureLabel like "GUDMAP:%"
	and i._Refs_key = 172505
	and i._Image_key = ii._Image_key
	''', 'auto')
    for r in results:
	key = r['figureLabel']
	imageLookup[key] = []
	imageLookup[key].append(r['_ImagePane_key'])
    #print imageLookup

    # For each line in the input file

    for line in inPrepFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	prepID = tokens[0]
	secondary = tokens[1]
	labelledWith = tokens[2]

	prepLookup[prepID] = tokens
    #print prepLookup

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	prepID = tokens[0]
	markerID = tokens[1]
	reference = tokens[2]
	assayType = tokens[3]
	reporter = tokens[4]
	assayNote = tokens[5]
	createdBy = tokens[6]

	# create one assay per record

	# write the prep information

	prepFile.write(str(assayKey) + TAB + \
		prepID + TAB + \
		prepLookup[prepID][1] + TAB + \
		prepLookup[prepID][2] + CRT)

	# write the assay information

        assayFile.write(str(assayKey) + TAB + \
            markerID + TAB + \
            reference + TAB + \
            assayType + TAB + \
            TAB + \
            assayNote + TAB + \
            createdBy + CRT)

	if not assayLookup.has_key(prepID):
	    assayLookup[prepID] = []
        assayLookup[prepID] = assayKey

	assayKey = assayKey + 1

    inAssayFile.close()

    # end of "for line in inAssayFile.readlines():"

    # write one specimen per assay

    specimenKey = 1
    prevProbeID = ''
    specimenLookup = {}
    specimenProbeLookup = {}

    for sline in inSpecimenFile.readlines():

        # Split the line into tokens
        tokens = string.split(sline[:-1], TAB)

	prepID = tokens[0]
	specimenID = tokens[1]
	genotype = tokens[2]
	age = tokens[3]
	ageNote = tokens[4]
	sex = tokens[5]
	fixation = tokens[6]
	embedding = tokens[7]
	specimenHybridization = tokens[8]
	specimenNote = tokens[9]

	if prepID != prevProbeID:
            specimenKey = 1
	    assayKey = assayLookup[prepID]
	    prevProbeID = prepID

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

	if not specimenLookup.has_key(specimenID):
	    specimenLookup[specimenID] = []
        specimenLookup[specimenID] = specimenKey
	specimenProbeLookup[specimenID] = prepID

        specimenKey = specimenKey + 1

    inSpecimenFile.close()

    # write one results per specimen

    resultKey = 1
    prevSpecimenID = ''

    for rline in inResultFile.readlines():

        # Split the line into tokens
        tokens = string.split(rline[:-1], TAB)

	specimenID = tokens[0]
	strength = tokens[1]
	pattern = tokens[2]
	emapID = tokens[3]
	resultNote = tokens[4]
	imageName = tokens[5]

	emap = emapLookup[emapID][0]
	structureName = emap['structure']
	structureTheilerStage = emap['stage']

	if imageLookup.has_key(imageName):
	    imageName = imageLookup[specimenID][0]
        else:
	    imageName = ''

	if specimenID != prevSpecimenID:
            resultKey = 1
	    specimenKey = specimenLookup[specimenID]
	    assayKey = assayLookup[specimenProbeLookup[specimenID]]
	    prevSpecimenID = specimenID

	resultsFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    str(resultKey) + TAB + \
	    strength + TAB + \
	    pattern + TAB + \
	    structureName + TAB + \
	    str(structureTheilerStage) + TAB + \
	    resultNote + TAB + \
	    str(imageName) + CRT)

        resultKey = resultKey + 1

    inResultFile.close()

    # end of "for line in inResultFile.readlines():"

#
# Main
#

init()
process()
exit(0)

