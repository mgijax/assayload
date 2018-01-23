#!/usr/local/bin/python

#
# Program: gelload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Assays into Gel Structures
#
#	GXD_ProbePrep
#	GXD_Assay
#	GXD_GelLane
#	GXD_GelLaneStructure
#	GXD_GelRow
#	GXD_GelBand
#	ACC_Accession
#
# Requirements Satisfied by This Program:
#
# Usage:
#	gelload.py
#
# Envvars:
#
# Inputs:
#
#	Probe Prep file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Probe MGI ID
#		field 3: Probe Prep Type
#		field 4: Hybridization
#		field 5: Labelled With
#		field 6: Visualized With
#
#       Assay file, a tab-delimited file in the format:
#               field 1: Assay #
#               field 2: MGI Marker Accession ID
#               field 3: Reference (J:#####)
#               field 4: Assay Type
#               field 5: Reporter Gene
#               field 6: Assay Note
#               field 7: Created By
#
#	Gel Lane file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Lane #
#		field 3: Lane Label
#		field 4: Genotype
#		field 5: RNA type
#		field 6: Control
#		field 7: Sample Amount
#		field 8: Sex
#		field 9: Age
#		field 10: Age Note
#		field 11: Lane Note
#		field 12: MGI Structure name
#		field 13: MGI Structure Theiler Stage
#
#	Gel Row/Band file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Lane #
#		field 3: Row #
#		field 4: Size
#		field 5: Units
#		field 6: Strength
#		field 7: Row Note
#		field 8: Band Note
#
# Outputs:
#
#       BCP files:
#
#	GXD_ProbePrep.bcp		Probe Prep records
#	GXD_Assay.bcp			Assay records
#	GXD_GelLane.bcp			Gel Lanes
#	GXD_GelLaneStructure.bcp	Gel Lane Structures
#	GXD_GelRow.bcp			Gel Rows
#	GXD_GelBand.bcp			Gel Bands
#       ACC_Accession.bcp               Accession records
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
# 01/20/2010 lec
#	- TR9560/TR9782; remove verifyPrepCoverage
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
datadir = os.environ['ASSAYLOADDATADIR']	# file which contains the data files

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
inGelLaneFile = ''        # file descriptor
inGelBandFile = ''        # file descriptor

inPrepFileName = datadir + '/RT_PCR_probeprep.txt'
inAssayFileName = datadir + '/RT_PCR_assay.txt'
inGelLaneFileName = datadir + '/RT_PCR_gellane.txt'
inGelBandFileName = datadir + '/RT_PCR_gelband.txt'

# output files

outPrepFile = ''	# file descriptor
outAssayFile = ''	# file descriptor
outGelLaneFile = ''	# file descriptor
outGelLaneStFile = ''	# file descriptor
outGelRowFile = ''	# file descriptor
outGelBandFile = ''	# file descriptor
outAccFile = ''         # file descriptor
outAssayNoteFile = ''   # file descriptor

probeprepTable = 'GXD_ProbePrep'
assayTable = 'GXD_Assay'
assaynoteTable = 'GXD_AssayNote'
gelLaneTable = 'GXD_GelLane'
gelLaneStTable = 'GXD_GelLaneStructure'
gelRowTable = 'GXD_GelRow'
gelBandTable = 'GXD_GelBand'
accTable = 'ACC_Accession'

outPrepFileName = datadir + '/' + probeprepTable + '.bcp'
outAssayFileName = datadir + '/' + assayTable + '.bcp'
outGelLaneFileName = datadir + '/' + gelLaneTable + '.bcp'
outGelLaneStFileName = datadir + '/' + gelLaneStTable + '.bcp'
outGelRowFileName = datadir + '/' + gelRowTable + '.bcp'
outGelBandFileName = datadir + '/' + gelBandTable + '.bcp'
outAccFileName = datadir + '/' + accTable + '.bcp'
outAssayNoteFileName = datadir + '/' + assaynoteTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

# primary keys

prepKey = 0		# GXD_ProbePrep._ProbePrep_key
assayKey = 0		# GXD_Assay._Assay_key
gelLaneKey = 0		# GXD_GelLane._GelLane_key
gelRowKey = 0		# GXD_GelRow._GelRow_key
gelBandKey = 0		# GXD_GelBand._GelBand_key
accKey = 0              # ACC_Accession._Accession_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart

# accession constants

assayMgiTypeKey = '8'	# Assay
mgiPrefix = "MGI:"	# Prefix for MGI accession ID
accLogicalDBKey = '1'	# Logical DB Key for MGI accession ID
accPrivate = '0'	# Private status for MGI accession ID (false)
accPreferred = '1'	# Preferred status MGI accession ID (true)

assayProbePrep = {}	# Assay ID/Probe Prep keys
assayAssay = {}		# Assay ID/Assay keys
assayGelLane = {}	# Assay ID/Lane ID and Lane keys

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
    global outGelLaneFile, outGelLaneStFile, outGelRowFile, outGelBandFile
    global inPrimerFile, inPrepFile, inAssayFile, inGelLaneFile, inGelBandFile
 
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
 
    fdate = mgi_utils.date('%m%d%Y')	# current date
    diagFileName = sys.argv[0] + '.' + fdate + '.diagnostics'
    errorFileName = sys.argv[0] + '.' + fdate + '.error'

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
        inGelLaneFile = open(inGelLaneFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inGelLaneFileName)

    try:
        inGelBandFile = open(inGelBandFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inGelBandFileName)

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
        outGelLaneFile = open(outGelLaneFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outGelLaneFileName)

    try:
        outGelLaneStFile = open(outGelLaneStFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outGelLaneStFileName)

    try:
        outGelRowFile = open(outGelRowFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outGelRowFileName)

    try:
        outGelBandFile = open(outGelBandFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outGelBandFileName)

    try:
        outAccFile = open(outAccFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outAccFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

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

    global accKey, mgiKey, prepKey, assayKey
    global gelLaneKey, gelRowKey, gelBandKey

    results = db.sql('select maxKey = max(_ProbePrep_key) + 1 from GXD_ProbePrep', 'auto')
    prepKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Assay_key) + 1 from GXD_Assay', 'auto')
    assayKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_GelLane_key) + 1 from GXD_GelLane', 'auto')
    gelLaneKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_GelRow_key) from GXD_GelRow', 'auto')
    gelRowKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_GelBand_key) + 1 from GXD_GelBand', 'auto')
    gelBandKey = results[0]['maxKey']

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

    outPrepFile.close()
    outAssayFile.close()
    outAssayNoteFile.close()
    outGelLaneFile.close()
    outGelLaneStFile.close()
    outGelRowFile.close()
    outGelBandFile.close()
    outAccFile.close()

    # update the max Accession ID value
    db.sql('select * from ACC_setMax (%d)' % (recordsProcessed), None)

    db.commit()
    db.useOneConnection(0)

    if DEBUG or not bcpon:
        return

    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'
    currentDir = os.getcwd()

    bcp1 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), probeprepTable, currentDir, outPrepFileName)
    bcp2 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), assayTable, currentDir, outAssayFileName)
    bcp3 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), assaynoteTable, currentDir, outAssayNoteFileName)
    bcp4 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), gelLaneTable, currentDir, outGelLaneFileName)
    bcp5 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), gelLaneStTable, currentDir, outGelLaneStFileName)
    bcp6 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), gelRowTable, currentDir, outGelRowFileName)
    bcp7 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), gelBandTable, currentDir, outGelBandFileName)
    bcp8 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
	% (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), accTable, currentDir, outAccFileName)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    return

# Purpose:  processes probe prep data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processPrepFile():

    global assayProbePrep, prepKey

    lineNum = 0
    # For each line in the input file

    for line in inPrepFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    probeID = tokens[1]
	    prepType = tokens[2]
	    hybridization = tokens[3]
	    labelledWith = tokens[4]
	    visualization = tokens[5]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	if gxdloadlib.verifyPrepType(prepType, lineNum, errorFile) == 0:
	    error = 1

	probeKey = loadlib.verifyProbe(probeID, lineNum, errorFile)
	senseKey = gxdloadlib.verifyPrepSense(hybridization, lineNum, errorFile)
	labelKey = gxdloadlib.verifyPrepLabel(labelledWith, lineNum, errorFile)
	visualizationKey = gxdloadlib.verifyPrepVisualization(visualization, lineNum, errorFile)

        if probeKey == 0 or senseKey == 0 or labelKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

        outPrepFile.write(str(prepKey) + TAB + \
	    str(probeKey) + TAB + \
	    str(senseKey) + TAB + \
	    str(labelKey) + TAB + \
	    str(visualizationKey) + TAB + \
	    prepType + TAB + \
	    loaddate + TAB + loaddate + CRT)

	assayProbePrep[assayID] = prepKey
        prepKey = prepKey + 1

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

	if assayProbePrep.has_key(assayID):
	    probePrepKey = assayProbePrep[assayID]
	else:
	    probePrepKey = ''

        # if no errors, process

        outAssayFile.write(str(assayKey) + TAB + \
	    str(assayTypeKey) + TAB + \
	    str(referenceKey) + TAB + \
	    str(markerKey) + TAB + \
	    str(probePrepKey) + TAB + \
	    TAB + \
	    TAB + \
            str(reporterGeneKey) + TAB + \
            str(createdByKey) + TAB + \
            str(createdByKey) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	if len(note) > 0:
	    i = 0
	    while i < len(note):
		outAssayNoteFile.write(str(assayKey) + TAB + \
		    note[i:i+ASSAY_NOTE_LENGTH] + TAB + \
		    loaddate + TAB + loaddate + CRT)
		i = i + ASSAY_NOTE_LENGTH

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

# Purpose:  processes gel lane data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processGelLaneFile():

    global assayGelLane, gelLaneKey

    lineNum = 0

    # For each line in the input file

    for line in inGelLaneFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    laneID = tokens[1]
	    laneLabel = tokens[2]
	    genotypeID = tokens[3]
	    rnaType = tokens[4]
	    control = tokens[5]
	    sampleAmount = tokens[6]
	    gender = tokens[7]
	    age = tokens[8]
	    ageNote = tokens[9]
	    laneNote = tokens[10]
	    emapaID = tokens[11]
	    structureTS = tokens[12]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	# if control is set to "No", then there *is* a structure
	# else there are no structures

	hasStructure = 0
	if control == "No":
	    hasStructure = 1

	genotypeKey = gxdloadlib.verifyGenotype(genotypeID, lineNum, errorFile)
	rnaTypeKey = gxdloadlib.verifyGelRNAType(rnaType, lineNum, errorFile)
	controlKey = gxdloadlib.verifyGelControl(control, lineNum, errorFile)
	ageMin, ageMax = agelib.ageMinMax(age)

	if hasStructure:
	    structureKey = gxdloadlib.verifyTerm(emapaID, 90, '', lineNum, errorFile)
	    if structureKey == 0:
                error = 1

	#
	# if age = "Not Specified", then ageMin/ageMax = -1 which is < 0
	# so, removed this check:
	#	ageMin < 0 or ageMax < 0:
	#

        if genotypeKey == 0 or rnaTypeKey == 0 or controlKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

	key = '%s:%s' % (assayID, laneID)

	# if this is a lane that has not been added to the gel lane yet...

	if not assayGelLane.has_key(key):

            outGelLaneFile.write(
	        str(gelLaneKey) + TAB + \
	        str(assayAssay[assayID]) + TAB + \
	        str(genotypeKey) + TAB + \
	        str(rnaTypeKey) + TAB + \
	        str(controlKey) + TAB + \
	        str(laneID) + TAB + \
	        laneLabel + TAB + \
	        mgi_utils.prvalue(sampleAmount) + TAB + \
	        gender + TAB + \
	        age + TAB + \
	        str(ageMin) + TAB + \
	        str(ageMax) + TAB + \
	        mgi_utils.prvalue(ageNote) + TAB + \
	        mgi_utils.prvalue(laneNote) + TAB + \
	        loaddate + TAB + loaddate + CRT)

	    if hasStructure:
	        outGelLaneStFile.write(
	            str(gelLaneKey) + TAB + \
	            str(structureKey) + TAB + \
	            loaddate + TAB + loaddate + CRT)

	    assayGelLane[key] = gelLaneKey
            gelLaneKey = gelLaneKey + 1

	# else if gel lanes has more than one structure...

	else:
	    if hasStructure:
	        outGelLaneStFile.write(
	            str(assayGelLane[key]) + TAB + \
	            str(structureKey) + TAB + \
	            loaddate + TAB + loaddate + CRT)

    #	end of "for line in inGelLaneFile.readlines():"

    #print assayGelLane

    return

# Purpose:  processes gel row/band data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processGelBandFile():

    global gelRowKey, gelBandKey

    lineNum = 0
    prevAssay = 0
    prevLane = 0
    prevRow = 0

    # For each line in the input file

    for line in inGelBandFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    laneID = tokens[1]
	    rowID = tokens[2]
	    bandSize = tokens[3]
	    bandUnits = tokens[4]
	    bandStrength = tokens[5]
	    rowNote = tokens[6]
	    bandNote = tokens[7]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	unitsKey = gxdloadlib.verifyGelUnits(bandUnits, lineNum, errorFile)
	strengthKey = gxdloadlib.verifyGelStrength(bandStrength, lineNum, errorFile)

        if unitsKey == 0 or strengthKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

	# new Assay means new Row

	if prevAssay != assayID:

          gelRowKey = gelRowKey + 1

          outGelRowFile.write(
	      str(gelRowKey) + TAB + \
	      str(assayAssay[assayID]) + TAB + \
	      str(unitsKey) + TAB + \
	      str(rowID) + TAB + \
	      mgi_utils.prvalue(bandSize) + TAB + \
	      mgi_utils.prvalue(rowNote) + TAB + \
	      loaddate + TAB + loaddate + CRT)

	  prevAssay = assayID

	# determine the lane key based on assayID and laneID
	key = '%s:%s' % (assayID, laneID)
	laneKey = assayGelLane[key]

	outGelBandFile.write(
	    str(gelBandKey) + TAB + \
	    str(laneKey) + TAB + \
	    str(gelRowKey) + TAB + \
	    str(strengthKey) + TAB + \
	    mgi_utils.prvalue(bandNote) + TAB + \
	    loaddate + TAB + loaddate + CRT)

        gelBandKey = gelBandKey + 1

    #	end of "for line in inGelLaneFile.readlines():"

    return

def process():

    processPrepFile()
    recordsProcessed = processAssayFile()
    processGelLaneFile()
    processGelBandFile()
    bcpFiles(recordsProcessed)

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
process()
exit(0)

