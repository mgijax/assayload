#!/usr/local/bin/python

#
# Program: immunoload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Assays into Immuno Structures
#
#	GXD_AntibodyPrep
#	GXD_Assay
#	GXD_AssayNote***new
#	GXD_Specimen
#	GXD_InSituResults
#	GXD_ISResultStructure
#	ACC_Accession
#
# Requirements Satisfied by This Program:
#
# Usage:
#	immunoload.py
#
# Envvars:
#
# Inputs:
#
#	Antibody Prep file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Antibody MGI ID
#		field 3: Secondary Antibody
#		field 4: Labelled With
#
#	Assay file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: MGI Marker Accession ID
#		field 3: Reference (J:#####)
#		field 4: Assay Type
#		field 5: Reporter Gene
#		field 6: Assay Note
#		field 7: Created By
#
#	Specimen file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Specimen #
#		field 3: Specimen Label
#		field 4: Genotype ID
#		field 5: Age
#		field 6: Age Note
#		field 7: Sex
#		field 8: Fixation
#		field 9: Embedding Method
#		field 10: Hybridization
#		field 11: Specimen Note
#
#	Specimen Results file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Specimen #
#		field 3: Result #
#		field 4: Strength
#		field 5: Pattern
#		field 6: MGI Structure Name
#		field 7: MGI Structure Theiler Stage
#		field 8: Result Note
#		field 9: Comma-Separated Image Names (if any)
#
# Outputs:
#
#       BCP files:
#
#	GXD_AntibodyPrep.bcp		Antibody Prep records
#	GXD_Assay.bcp			Assay records
#	GXD_AssayNote.bcp		Assay Note records
#	GXD_Specimen.bcp		Specimens
#	GXD_InSituResult.bcp		InSitu Results
#	GXD_ISResultStructure.bcp	InSitu Result Structures
#       ACC_Accession.bcp               Accession records
#	Result_Image.txt		Results/Image text file (see gxdimageload/assocResultImage)
#		field 1: Result key
#		field 2: Figure Label
#	
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# Exit Codes:
#
# Assumes:
#
#	That no one else is adding records to the database.
#
# Bugs:
#
# Implementation:
#
# History
#
# 03/04/2014 lec
#	- TR11471/new/copied from insituload.py
#

import sys
import os
import string
import db
import mgi_utils
import agelib
import loadlib
import gxdloadlib
import gxdexpression

#globals

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
mode = os.environ['ASSAYLOADMODE']
datadir = os.environ['ASSAYLOADDATADIR']

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
bcpdelim = TAB		# bcp file delimiter

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor

# input files

inPrepFile = ''           # file descriptor
inAssayFile = ''          # file descriptor
inAssayNoteFile = ''      # file descriptor
inSpecimenFile = ''       # file descriptor
inResultsFile = ''        # file descriptor

inPrepFileName = datadir + '/Immuno_prep.txt'
inAssayFileName = datadir + '/Immuno_assay.txt'
inSpecimenFileName = datadir + '/Immuno_specimen.txt'
inResultsFileName = datadir + '/Immuno_results.txt'

# output files

outPrepFile = ''	# file descriptor
outAssayFile = ''	# file descriptor
outAssayNoteFile = ''	# file descriptor
outSpecimenFile = ''	# file descriptor
outResultStFile = ''	# file descriptor
outResultFile = ''	# file descriptor
outAccFile = ''         # file descriptor
outResultImageFile = '' # file descriptor

prepTable = 'GXD_AntibodyPrep'
assayTable = 'GXD_Assay'
assaynoteTable = 'GXD_AssayNote'
specimenTable = 'GXD_Specimen'
resultTable = 'GXD_InSituResult'
resultStTable = 'GXD_ISResultStructure'
accTable = 'ACC_Accession'
resultImageTable = 'GXD_InSituResultImage'

outPrepFileName = datadir + '/' + prepTable + '.bcp'
outAssayFileName = datadir + '/' + assayTable + '.bcp'
outAssayNoteFileName = datadir + '/' + assaynoteTable + '.bcp'
outSpecimenFileName = datadir + '/' + specimenTable + '.bcp'
outResultFileName = datadir + '/' + resultTable + '.bcp'
outResultStFileName = datadir + '/' + resultStTable + '.bcp'
outAccFileName = datadir + '/' + accTable + '.bcp'
outResultImageFileName = datadir + '/' + resultImageTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

# primary keys

antibodyPrepKey = 0	# GXD_AntibodyPrep._AntibodyPrep_key
assayKey = 0		# GXD_Assay._Assay_key
specimenKey = 0		# GXD_GelLane._GelLane_key
resultKey = 0		# GXD_GelRow._GelRow_key
accKey = 0              # ACC_Accession._Accession_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart

# accession constants

assayMgiTypeKey = '8'   # Assay
mgiPrefix = "MGI:"      # Prefix for MGI accession ID
accLogicalDBKey = '1'   # Logical DB Key for MGI accession ID
accPrivate = '0'        # Private status for MGI accession ID (false)
accPreferred = '1'      # Preferred status MGI accession ID (true)

assayPrep = {}	# Assay ID/Probe Prep keys
assayAssay= {}		# Assay ID/Assay keys
assaySpecimen = {}	# Assay ID/Specimen ID and Specimen keys

ASSAY_NOTE_LENGTH = 255

loaddate = loadlib.loaddate

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
 
    try:
        diagFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        errorFile.write('\n\nEnd Date/Time: %s\n' % (mgi_utils.date()))
        diagFile.close()
        errorFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
# Purpose: process command line options
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global diagFile, errorFile, errorFileName, diagFileName
    global outAccFile, outPrepFile, outAssayFile, outAssayNoteFile
    global outSpecimenFile, outResultStFile, outResultFile, outResultImageFile
    global inPrepFile, inAssayFile, inSpecimenFile, inResultsFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    diagFileName = datadir + '/immunoload.diagnostics'
    errorFileName = datadir + '/immunoload.error'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
		
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
		
    # Input Files

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
        inResultsFile = open(inResultsFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResultsFileName)

    # Output Files

    try:
        outPrepFile = open(outPrepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outPrepFileName)

    try:
        outAssayFile = open(outAssayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outAssayFileName)

    try:
        outAssayNoteFile = open(outAssayNoteFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outAssayNoteFileName)

    try:
        outSpecimenFile = open(outSpecimenFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outSpecimenFileName)

    try:
        outResultStFile = open(outResultStFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outResultStFileName)

    try:
        outResultFile = open(outResultFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outResultFileName)

    try:
        outAccFile = open(outAccFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outAccFileName)

    try:
        outResultImageFile = open(outResultImageFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outResultImageFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlDatabase()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    return

# Purpose: verify processing mode
# Returns: nothing
# Assumes: nothing
# Effects: if the processing mode is not valid, exits.
#	   else, sets global variables
# Throws:  nothing

def verifyMode():

    global DEBUG

    if mode == 'preview':
        DEBUG = 1
        bcpon = 0
    elif mode != 'load':
        exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global accKey, mgiKey, antibodyPrepKey, assayKey
    global specimenKey, resultKey

    results = db.sql('select maxKey = max(_AntibodyPrep_key) + 1 from GXD_AntibodyPrep', 'auto')
    antibodyPrepKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Assay_key) + 1 from GXD_Assay', 'auto')
    assayKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Specimen_key) + 1 from GXD_Specimen', 'auto')
    specimenKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Result_key) from GXD_InSituResult', 'auto')
    resultKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Accession_key) + 1 from ACC_Accession', 'auto')
    accKey = results[0]['maxKey']

    results = db.sql('select maxKey = maxNumericPart + 1 from ACC_AccessionMax ' + \
        'where prefixPart = "%s"' % (mgiPrefix), 'auto')
    mgiKey = results[0]['maxKey']

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles(
   recordsProcessed	# number of records processed (integer)
   ):

    if DEBUG or not bcpon:
        return

    outPrepFile.close()
    outAssayFile.close()
    outAssayNoteFile.close()
    outSpecimenFile.close()
    outResultStFile.close()
    outResultFile.close()
    outAccFile.close()
    outResultImageFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"%s' % (bcpdelim) + '" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, prepTable, outPrepFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, assayTable, outAssayFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, assaynoteTable, outAssayNoteFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, specimenTable, outSpecimenFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, resultStTable, outResultStFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, resultTable, outResultFileName, bcpII)
    bcp7 = '%s%s in %s %s' % (bcpI, accTable, outAccFileName, bcpII)
    bcp8 = '%s%s in %s %s' % (bcpI, resultImageTable, outResultImageFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    # load the cache tables for the records processed (by assay Key)

    for i in range(assayKey - recordsProcessed, assayKey + 1):
	# run cache load by assayKey
        gxdexpression.process(i)

    # update the max Accession ID value
    db.sql('exec ACC_setMax %d' % (recordsProcessed), None)

    # update statistics
    db.sql('update statistics %s' % (prepTable), None)
    db.sql('update statistics %s' % (assayTable), None)
    db.sql('update statistics %s' % (assaynoteTable), None)
    db.sql('update statistics %s' % (specimenTable), None)
    db.sql('update statistics %s' % (resultStTable), None)
    db.sql('update statistics %s' % (resultTable), None)
    db.sql('update statistics %s' % (resultImageTable), None)

    return

# Purpose:  processes antibody prep data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processPrepFile():

    global assayPrep, antibodyPrepKey

    # This dictionary is used to keep track of each combination of antibody key,
    # secondary key, label key that are added.
    # If the combination exists on multiple records in the antibody prep input
    # file, only one antibody prep record will be created in the database and
    # it will be shared by multiple assays.
    prepLookup = {}

    lineNum = 0
    # For each line in the input file

    for line in inPrepFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    prepID = tokens[1]
	    secondary = tokens[2]
	    labelledWith = tokens[3]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	antibodyKey = gxdloadlib.verifyAntibody(prepID, lineNum, errorFile)
	secondaryKey = gxdloadlib.verifyPrepSecondary(secondary, lineNum, errorFile)
	labelKey = gxdloadlib.verifyPrepLabel(labelledWith, lineNum, errorFile)

        if antibodyKey == 0:
            # set error flag to true
	    errorFile.write('\ngxdloadlib.verifyAntibody  %s' % (prepID))
            error = 1

	if secondaryKey == 0: 
            # set error flag to true
	    errorFile.write('\ngxdloadlib.verifyPrepSecondary:  %s' % (secondary))
            error = 1

	if labelKey == 0:
            # set error flag to true
	    errorFile.write('\ngxdloadlib.verifyPrepLabel:  %s' % (labelledWith))
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

	# Determine if the current combination of antibody key, secondary key,
        # label key has already been added
        # to the output file.
        #
	key = '%s:%s:%s' % (str(antibodyKey),
                                  str(secondaryKey),
                                  str(labelKey))

        # If a probe prep record has already been created, add the existing
        # probe prep key to the lookup for the current assayID.
        #
	if prepLookup.has_key(key):
	    assayPrep[assayID] = prepLookup[key]

        # Otherwise, add a new probe prep key to the lookup for the current
        # assayID and also add a new entry to the dictionary for this
        # combination of probe key, secondary key, label key. 
        #
	else:
	    outPrepFile.write(str(antibodyPrepKey) + TAB + \
	        str(antibodyKey) + TAB + \
	        str(secondaryKey) + TAB + \
	        str(labelKey) + TAB + \
	        loaddate + TAB + loaddate + CRT)

	    assayPrep[assayID] = antibodyPrepKey
	    prepLookup[key] = antibodyPrepKey
	    antibodyPrepKey = antibodyPrepKey + 1

    #	end of "for line in inPrepFile.readlines():"

    return

# Purpose:  processes assay data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processAssayFile():

    global assayAssay, assayKey, accKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inAssayFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    markerID = tokens[1]
	    jnum = tokens[2]
	    assayType = tokens[3]
	    reporterGene = tokens[4]
	    note = tokens[5]
	    createdBy = tokens[6]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)
        referenceKey = loadlib.verifyReference(jnum, lineNum, errorFile)
	assayTypeKey = gxdloadlib.verifyAssayType(assayType, lineNum, errorFile)
	createdByKey = loadlib.verifyUser(createdBy, lineNum, errorFile)

        if markerKey == 0 or referenceKey == 0 or assayTypeKey == 0:
            # set error flag to true
            error = 1

        if len(reporterGene) > 0:
            reporterGeneKey = gxdloadlib.verifyReporterGene(reporterGene, lineNum, errorFile)
	    if reporterGeneKey == 0:
                error = 1
        else:
            reporterGeneKey = ''

        # if errors, continue to next record
        if error:
            continue

	if assayPrep.has_key(assayID):
	    antibodyPrepKey = assayPrep[assayID]
	else:
	    antibodyPrepKey = ''

        # if no errors, process

        outAssayFile.write(str(assayKey) + TAB + \
	    str(assayTypeKey) + TAB + \
	    str(referenceKey) + TAB + \
	    str(markerKey) + TAB + \
	    TAB + \
	    str(antibodyPrepKey) + TAB + \
	    TAB + \
            str(reporterGeneKey) + TAB + \
            str(createdByKey) + TAB + \
            str(createdByKey) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	if len(note) > 0:
	    i = 0
	    sequenceNum = 1
	    while i < len(note):
		outAssayNoteFile.write(str(assayKey) + TAB + \
		    str(sequenceNum) + TAB + \
		    note[i:i+ASSAY_NOTE_LENGTH] + TAB + \
		    loaddate + TAB + loaddate + CRT)
		i = i + ASSAY_NOTE_LENGTH
		sequenceNum = sequenceNum + 1

        # MGI Accession ID for the assay

	outAccFile.write(str(accKey) + TAB + \
	    mgiPrefix + str(mgiKey) + TAB + \
	    mgiPrefix + TAB + \
	    str(mgiKey) + TAB + \
	    accLogicalDBKey + TAB + \
	    str(assayKey) + TAB + \
	    assayMgiTypeKey + TAB + \
	    accPrivate + TAB + \
	    accPreferred + TAB + \
            str(createdByKey) + TAB + \
            str(createdByKey) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	assayAssay[assayID] = assayKey
	accKey = accKey + 1
	mgiKey = mgiKey + 1
        assayKey = assayKey + 1

    #	end of "for line in inAssayFile.readlines():"

    return lineNum

# Purpose:  processes specimen data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processSpecimenFile():

    global assaySpecimen, specimenKey

    lineNum = 0
    # For each line in the input file

    for line in inSpecimenFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    specimenID = tokens[1]
	    specimenLabel = tokens[2]
	    genotypeID = tokens[3]
	    age = tokens[4]
	    ageNote = tokens[5]
	    gender = tokens[6]
	    fixation = tokens[7]
	    embedding = tokens[8]
	    hybridization = tokens[9]
	    specimenNote = tokens[10]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	if gxdloadlib.verifyHybridization(hybridization, lineNum, errorFile) == 0:
	    error = 1

	genotypeKey = gxdloadlib.verifyGenotype(genotypeID, lineNum, errorFile)
	fixationKey = gxdloadlib.verifyFixationMethod(fixation, lineNum, errorFile)
	embeddingKey = gxdloadlib.verifyEmbeddingMethod(embedding, lineNum, errorFile)
	ageMin, ageMax = agelib.ageMinMax(age)

        if genotypeKey == 0 or ageMin < 0 or ageMax < 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

        outSpecimenFile.write(
	    str(specimenKey) + TAB + \
	    str(assayAssay[assayID]) + TAB + \
	    str(embeddingKey) + TAB + \
	    str(fixationKey) + TAB + \
	    str(genotypeKey) + TAB + \
	    specimenID + TAB + \
	    specimenLabel + TAB + \
	    gender + TAB + \
	    age + TAB + \
	    str(ageMin) + TAB + \
	    str(ageMax) + TAB + \
	    mgi_utils.prvalue(ageNote) + TAB + \
	    hybridization + TAB + \
	    mgi_utils.prvalue(specimenNote) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	key = '%s:%s' % (assayID, specimenID)
	assaySpecimen[key] = specimenKey
        specimenKey = specimenKey + 1

    #	end of "for line in inSpecimenFile.readlines():"

    return

# Purpose:  processes results data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processResultsFile():

    global resultKey

    prevAssay = 0
    prevSpecimen = 0
    prevResult = 0
    lineNum = 0
    # For each line in the input file

    for line in inResultsFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    specimenID = tokens[1]
	    resultID = tokens[2]
	    strength = tokens[3]
	    pattern = tokens[4]
	    structureName = tokens[5]
	    structureTS = tokens[6]
	    resultNote = tokens[7]
	    images = tokens[8]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	strengthKey = gxdloadlib.verifyStrength(strength, lineNum, errorFile)
	patternKey = gxdloadlib.verifyPattern(pattern, lineNum, errorFile)

	structureName = string.replace(structureName, '"', '""')
	structureKey = gxdloadlib.verifyStructure(structureName, structureTS, lineNum, errorFile)

        if strengthKey == 0 or patternKey == 0 or structureKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

	key = '%s:%s' % (assayID, specimenID)

	if not assaySpecimen.has_key(key):
	    errorFile.write('Cannot find Assay:Speciman key "%s"\n' % (key))
	    continue

	specimenKey = assaySpecimen[key]

	if prevAssay != assayID:
	    prevSpecimen = 0

        if prevSpecimen != specimenKey:
	    prevResult = 0

	if prevResult != resultID:

            resultKey = resultKey + 1

            outResultFile.write(
	        str(resultKey) + TAB + \
	        str(specimenKey) + TAB + \
	        str(strengthKey) + TAB + \
	        str(patternKey) + TAB + \
	        resultID + TAB + \
	        mgi_utils.prvalue(resultNote) + TAB + \
	        loaddate + TAB + loaddate + CRT)

            for image in string.split(images,','):
                if image != '':
                    outResultImageFile.write(str(resultKey) + TAB + \
					image + TAB + \
	        			loaddate + TAB + loaddate + CRT)

	outResultStFile.write(
	    str(resultKey) + TAB + \
	    str(structureKey) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	prevAssay = assayID
	prevSpecimen = specimenKey
	prevResult = resultID

    #	end of "for line in inResultsFile.readlines():"

    return

def process():

    processPrepFile()
    recordsProcessed = processAssayFile()
    processSpecimenFile()
    processResultsFile()
    bcpFiles(recordsProcessed)

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
process()
exit(0)
