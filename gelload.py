#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: gelload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Assays into Gel Structures
#
#	PRB_Probe (primers)
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
#	program.py
#	-S = database server
#	-D = database
#	-U = user
#	-P = password file
#	-M = mode
#
# Envvars:
#
# Inputs:
#
#       Primer file, a tab-delimited file in the format:
#		field 1: Assay #
#               field 2: MGI Marker Accession ID
#               field 3: Primer Name
#               field 4: Reference (J:#####)
#               field 5: Region Covered
#               field 6: Sequence 1
#               field 7: Sequence 2
#               field 8: Repeat Unit
#               field 9: More Than One Product (y/n)
#               field 10: Product Size
#
#	Probe Prep file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Probe Prep Type
#		field 3: Hybridization
#		field 4: Labelled With
#		field 5: Label Coverage
#		field 6: Visualized With
#
#	Assay file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: MGI Marker Accession ID
#		field 3: Reference (J:#####)
#		field 4: Assay Type
#		field 5: Created By
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
#       PRB_Probe.bcp                   master Primer records
#	PRB_Marker.bcp			Primer/Marker records
#       PRB_Reference.bcp         	Primer Reference records
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

import sys
import os
import string
import getopt
import db
import mgi_utils
import agelib
import loadlib
import gxdloadlib

#globals

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
bcpdelim = TAB		# bcp file delimiter

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

datadir = os.environ['RTPCRDATADIR']	# file which contains the data files

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor

# input files

inPrimerFile = ''         # file descriptor
inPrepFile = ''           # file descriptor
inAssayFile = ''          # file descriptor
inGelLaneFile = ''        # file descriptor
inGelBandFile = ''        # file descriptor

inPrimerFileName = datadir + '/RT_PCR_primer.txt'
inPrepFileName = datadir + '/RT_PCR_probeprep.txt'
inAssayFileName = datadir + '/RT_PCR_assay.txt'
inGelLaneFileName = datadir + '/RT_PCR_gellane.txt'
inGelBandFileName = datadir + '/RT_PCR_gelband.txt'

# output files

outPrimerFile = ''      # file descriptor
outMarkerFile = ''	# file descriptor
outRefFile = ''         # file descriptor
outPrepFile = ''	# file descriptor
outAssayFile = ''	# file descriptor
outGelLaneFile = ''	# file descriptor
outGelLaneStFile = ''	# file descriptor
outGelRowFile = ''	# file descriptor
outGelBandFile = ''	# file descriptor
outAccFile = ''         # file descriptor

probeTable = 'PRB_Probe'
markerTable = 'PRB_Marker'
refTable = 'PRB_Reference'
probeprepTable = 'GXD_ProbePrep'
assayTable = 'GXD_Assay'
gelLaneTable = 'GXD_GelLane'
gelLaneStTable = 'GXD_GelLaneStructure'
gelRowTable = 'GXD_GelRow'
gelBandTable = 'GXD_GelBand'
accTable = 'ACC_Accession'

outPrimerFileName = datadir + '/' + probeTable + '.bcp'
outMarkerFileName = datadir + '/' + markerTable + '.bcp'
outRefFileName = datadir + '/' + refTable + '.bcp'
outPrepFileName = datadir + '/' + probeprepTable + '.bcp'
outAssayFileName = datadir + '/' + assayTable + '.bcp'
outGelLaneFileName = datadir + '/' + gelLaneTable + '.bcp'
outGelLaneStFileName = datadir + '/' + gelLaneStTable + '.bcp'
outGelRowFileName = datadir + '/' + gelRowTable + '.bcp'
outGelBandFileName = datadir + '/' + gelBandTable + '.bcp'
outAccFileName = datadir + '/' + accTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
passwordFileName = ''	# password file name

mode = ''		# processing mode (load, preview)

createdBy = os.environ['CREATEDBY']

# primary keys

probeKey = 0            # PRB_Probe._Probe_key
refKey = 0		# PRB_Reference._Reference_key
prepKey = 0		# GXD_ProbePrep._ProbePrep_key
assayKey = 0		# GXD_Assay._Assay_key
gelLaneKey = 0		# GXD_GelLane._GelLane_key
gelRowKey = 0		# GXD_GelRow._GelRow_key
gelBandKey = 0		# GXD_GelBand._GelBand_key
accKey = 0              # ACC_Accession._Accession_key
mgiKey = 0              # ACC_AccessionMax.maxNumericPart

# primer constants

dnaType = 'primer'	# PRB_Probe.DNAtype
relationship = 'A'	# PRB_Marker.relationship
NA = '-2'		# for Not Applicable fields
primerMgiTypeKey = '3'	# Molecular Segment

# accession constants

assayMgiTypeKey = '8'	# Assay
mgiPrefix = "MGI:"	# Prefix for MGI accession ID
accLogicalDBKey = '1'	# Logical DB Key for MGI accession ID
accPrivate = '0'	# Private status for MGI accession ID (false)
accPreferred = '1'	# Preferred status MGI accession ID (true)

assayPrimer = {}	# Assay ID/Primer keys
assayProbePrep = {}	# Assay ID/Probe Prep keys
assayAssay= {}		# Assay ID/Assay keys
assayGelLane = {}	# Assay ID/Lane ID and Lane keys

loaddate = loadlib.loaddate

# Purpose: displays correct usage of this program
# Returns: nothing
# Assumes: nothing
# Effects: exits with status of 1
# Throws: nothing
 
def showUsage():
    usage = 'usage: %s -S server\n' % sys.argv[0] + \
        '-D database\n' + \
        '-U user\n' + \
        '-P password file\n' + \
        '-M mode\n'

    exit(1, usage)
 
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
#          calls showUsage() if usage error
#          exits if files cannot be opened
# Throws: nothing

def init():
    global diagFile, errorFile, errorFileName, diagFileName, passwordFileName
    global mode
    global outPrimerFile, outMarkerFile, outRefFile, outAccFile, outPrepFile, outAssayFile
    global outGelLaneFile, outGelLaneStFile, outGelRowFile, outGelBandFile
    global inPrimerFile, inPrepFile, inAssayFile, inGelLaneFile, inGelBandFile
 
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:M:')
    except:
        showUsage()
 
    #
    # Set server, database, user, passwords depending on options specified
    #
 
    server = ''
    database = ''
    user = ''
    password = ''
 
    for opt in optlist:
        if opt[0] == '-S':
            server = opt[1]
        elif opt[0] == '-D':
            database = opt[1]
        elif opt[0] == '-U':
            user = opt[1]
        elif opt[0] == '-P':
            passwordFileName = opt[1]
        elif opt[0] == '-M':
            mode = opt[1]
        else:
            showUsage()

    # User must specify Server, Database, User and Password
    password = string.strip(open(passwordFileName, 'r').readline())
    if server == '' or database == '' or user == '' or password == '' \
	or mode == '':
        showUsage()

    # Initialize db.py DBMS parameters
    db.set_sqlLogin(user, password, server, database)
    db.useOneConnection(1)
 
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
        inPrimerFile = open(inPrimerFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inPrimerFileName)

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
        outPrimerFile = open(outPrimerFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outPrimerFileName)

    try:
        outMarkerFile = open(outMarkerFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outMarkerFileName)

    try:
        outRefFile = open(outRefFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outRefFileName)

    try:
        outPrepFile = open(outPrepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outPrepFileName)

    try:
        outAssayFile = open(outAssayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outAssayFileName)

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

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (server))
    diagFile.write('Database: %s\n' % (database))
    diagFile.write('User: %s\n' % (user))

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

    global probeKey, refKey, accKey, mgiKey, prepKey, assayKey
    global gelLaneKey, gelRowKey, gelBandKey

    results = db.sql('select maxKey = max(_Probe_key) + 1 from PRB_Probe', 'auto')
    probeKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Reference_key) + 1 from PRB_Reference', 'auto')
    refKey = results[0]['maxKey']

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

    if DEBUG or not bcpon:
        return

    outPrimerFile.close()
    outMarkerFile.close()
    outRefFile.close()
    outPrepFile.close()
    outAssayFile.close()
    outGelLaneFile.close()
    outGelLaneStFile.close()
    outGelRowFile.close()
    outGelBandFile.close()
    outAccFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"%s' % (bcpdelim) + '" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())

    bcp1 = '%s%s in %s %s' % (bcpI, probeTable, outPrimerFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, markerTable, outMarkerFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, refTable, outRefFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, probeprepTable, outPrepFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, assayTable, outAssayFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, gelLaneTable, outGelLaneFileName, bcpII)
    bcp7 = '%s%s in %s %s' % (bcpI, gelLaneStTable, outGelLaneStFileName, bcpII)
    bcp8 = '%s%s in %s %s' % (bcpI, gelRowTable, outGelRowFileName, bcpII)
    bcp9 = '%s%s in %s %s' % (bcpI, gelBandTable, outGelBandFileName, bcpII)
    bcp10 = '%s%s in %s %s' % (bcpI, accTable, outAccFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6, bcp7, bcp8, bcp9, bcp10]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    # load the cache tables for the records processed (by assay Key)

    for i in range(assayKey - recordsProcessed, assayKey + 1):
	db.sql('exec GXD_loadCacheByAssay %d' % (i), None)

    # update the max Accession ID value
    db.sql('exec ACC_setMax %d' % (recordsProcessed), None)

    # update statistics
    db.sql('update statistics %s' % (probeTable), None)
    db.sql('update statistics %s' % (markerTable), None)
    db.sql('update statistics %s' % (refTable), None)
    db.sql('update statistics %s' % (probeprepTable), None)
    db.sql('update statistics %s' % (assayTable), None)
    db.sql('update statistics %s' % (gelLaneTable), None)
    db.sql('update statistics %s' % (gelLaneStTable), None)
    db.sql('update statistics %s' % (gelRowTable), None)
    db.sql('update statistics %s' % (gelBandTable), None)
    db.sql('update statistics %s' % (accTable), None)

    return

# Purpose:  processes primer data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def processPrimerFile():

    global assayPrimer
    global probeKey, refKey, accKey, mgiKey

    lineNum = 0
    # For each line in the input file

    for line in inPrimerFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

        try:
	    assayID = tokens[0]
	    markerID = tokens[1]
	    name = tokens[2]
	    jnum = tokens[3]
	    regionCovered = tokens[4]
	    sequence1 = tokens[5]
	    sequence2 = tokens[6]
	    repeatUnit = tokens[7]
	    moreProduct = tokens[8]
	    productSize = tokens[9]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)
        referenceKey = loadlib.verifyReference(jnum, lineNum, errorFile)

        if markerKey == 0 or referenceKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

        outPrimerFile.write(str(probeKey) + TAB + \
	    name + TAB + \
	    TAB + \
	    NA + TAB + \
	    NA + TAB + \
	    sequence1 + TAB + \
	    sequence2 + TAB + \
	    mgi_utils.prvalue(regionCovered) + TAB + \
	    TAB + \
	    TAB + \
	    TAB + \
	    dnaType + TAB + \
	    mgi_utils.prvalue(repeatUnit) + TAB + \
	    mgi_utils.prvalue(productSize) + TAB + \
	    moreProduct + TAB + \
	    loaddate + TAB + loaddate + CRT)

	outMarkerFile.write(str(probeKey) + TAB + \
	    str(markerKey) + TAB + \
	    relationship + TAB + \
	    loaddate + TAB + loaddate + CRT)

        outRefFile.write(str(refKey) + TAB + str(probeKey) + TAB + str(referenceKey) + TAB + \
	    TAB + '0' + TAB + '0' + TAB + loaddate + TAB + loaddate + CRT)

        # MGI Accession ID for the primer

	outAccFile.write(str(accKey) + TAB + \
	    mgiPrefix + str(mgiKey) + TAB + \
	    mgiPrefix + TAB + \
	    str(mgiKey) + TAB + \
	    accLogicalDBKey + TAB + \
	    str(probeKey) + TAB + \
	    primerMgiTypeKey + TAB + \
	    accPrivate + TAB + \
	    accPreferred + TAB + \
	    loaddate + TAB + loaddate + TAB + loaddate + CRT)

	assayPrimer[assayID] = probeKey

        accKey = accKey + 1
        mgiKey = mgiKey + 1
	refKey = refKey + 1
        probeKey = probeKey + 1

    #	end of "for line in inPrimerFile.readlines():"

    return lineNum

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
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    prepType = tokens[1]
	    hybridization = tokens[2]
	    labelledWith = tokens[3]
	    labelCoverage = tokens[4]
	    visualization = tokens[5]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	if gxdloadlib.verifyPrepType(prepType, lineNum, errorFile) == 0:
	    error = 1

	senseKey = gxdloadlib.verifyPrepSense(hybridization, lineNum, errorFile)
	labelKey = gxdloadlib.verifyPrepLabel(labelledWith, lineNum, errorFile)
	coverageKey = gxdloadlib.verifyPrepCoverage(labelCoverage, lineNum, errorFile)
	visualizationKey = gxdloadlib.verifyPrepVisualization(visualization, lineNum, errorFile)

        if senseKey == 0 or labelKey == 0 or coverageKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

        outPrepFile.write(str(prepKey) + TAB + \
	    str(assayPrimer[assayID]) + TAB + \
	    str(senseKey) + TAB + \
	    str(labelKey) + TAB + \
	    str(coverageKey) + TAB + \
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
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    markerID = tokens[1]
	    jnum = tokens[2]
	    assayType = tokens[3]
	    createdBy = tokens[4]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	markerKey = loadlib.verifyMarker(markerID, lineNum, errorFile)
        referenceKey = loadlib.verifyReference(jnum, lineNum, errorFile)
	assayTypeKey = gxdloadlib.verifyAssayType(assayType, lineNum, errorFile)

        if markerKey == 0 or referenceKey == 0 or assayTypeKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

        outAssayFile.write(str(assayKey) + TAB + \
	    str(assayTypeKey) + TAB + \
	    str(referenceKey) + TAB + \
	    str(markerKey) + TAB + \
	    str(assayProbePrep[assayID]) + TAB + \
	    TAB + \
	    TAB + \
	    TAB + \
	    createdBy + TAB + \
	    createdBy + TAB + \
	    loaddate + TAB + loaddate + CRT)

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
	    loaddate + TAB + loaddate + TAB + loaddate + CRT)

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
	    structureName = tokens[11]
	    structureTS = tokens[12]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	genotypeKey = gxdloadlib.verifyGenotype(genotypeID, lineNum, errorFile)
	rnaTypeKey = gxdloadlib.verifyGelRNAType(rnaType, lineNum, errorFile)
	controlKey = gxdloadlib.verifyGelControl(control, lineNum, errorFile)
	ageMin, ageMax = agelib.ageMinMax(age)
	structureKey = gxdloadlib.verifyStructure(structureName, structureTS, lineNum, errorFile)

        if genotypeKey == 0 or rnaTypeKey == 0 or controlKey == 0 or \
		ageMin < 0 or ageMax < 0 or structureKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

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

	outGelLaneStFile.write(
	    str(gelLaneKey) + TAB + \
	    str(structureKey) + TAB + \
	    loaddate + TAB + loaddate + CRT)

	key = '%s:%s' % (assayID, laneID)
	assayGelLane[key] = gelLaneKey
        gelLaneKey = gelLaneKey + 1

    #	end of "for line in inGelLaneFile.readlines():"

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

    recordsProcessed = processPrimerFile()
    processPrepFile()
    recordsProcessed = recordsProcessed + processAssayFile()
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

# $Log$
# Revision 1.12  2003/10/01 17:53:08  lec
# removed unnecessary imports
#
# Revision 1.11  2003/10/01 17:48:27  lec
# removed unnecessary imports
#
# Revision 1.10  2003/09/26 16:23:55  lec
# MGI 2.97
#
# Revision 1.9  2003/09/24 12:29:58  lec
# TR 5154
#
# Revision 1.8  2003/07/18 15:44:09  lec
# rtpcr.py
#
# Revision 1.7  2003/07/11 16:24:14  lec
# TR 4800
#
# Revision 1.6  2003/07/01 16:48:08  lec
# TR 4800
#
# Revision 1.5  2003/06/23 17:20:44  lec
# TR4800
#
# Revision 1.4  2003/06/18 15:56:16  lec
# TR 4800
#
# Revision 1.3  2003/06/18 13:19:32  lec
# TR 4800
#
# Revision 1.2  2003/06/17 12:17:49  lec
# TR 4800
#
