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
import accessionlib
import agelib

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

# dictionaries to cache data for quicker lookup

referenceDict = {}      # references
markerDict = {}      	# markers
senseDict = {}		# probe sense
labelDict = {}		# probe label
coverageDict = {}	# probe coverage
visualDict = {}		# probe visualization
assayTypeDict = {}	# assay type
gelRNATypeDict = {}	# gel rna types
gelControlDict = {}	# gel control
genotypeDict = {}       # genotypes
structureDict = {}	# anatomical structures
gelUnitsDict = {}	# gel units
gelStrengthDict = {}	# gel strength
prepTypeList = ['DNA', 'RNA', 'Not Specified'] 	# lookup of probe prep types

assayPrimer = {}	# Assay ID/Primer keys
assayProbePrep = {}	# Assay ID/Probe Prep keys
assayAssay= {}		# Assay ID/Assay keys
assayGelLane = {}	# Assay ID/Lane ID and Lane keys

cdate = mgi_utils.date('%m/%d/%Y')	# current date

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
    global diagFile, errorFile, inputFile, errorFileName, diagFileName, passwordFileName
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


# Purpose:  verify Marker Accession ID
# Returns:  Marker Key if Marker is valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Marker exists either in the marker dictionary or the database
#	writes to the error file if the Marker is invalid
#	adds the marker id and key to the marker dictionary if the Marker is valid
# Throws:  nothing

def verifyMarker(
    markerID, 	# Accession ID of the Marker (string)
    lineNum	# line number (integer)
    ):

    global markerDict

    markerKey = 0

    if markerDict.has_key(markerID):
        return markerDict[markerID]
    else:
        results = db.sql('select _Object_key from MRK_Acc_View where accID = "%s" ' % (markerID), 'auto')

        for r in results:
            if r['_Object_key'] is None:
                errorFile.write('Invalid Mouse Marker (%d) %s\n' % (lineNum, markerID))
                markerKey = 0
            else:
                markerKey = r['_Object_key']
                markerDict[markerID] = markerKey

    return markerKey

# Purpose:  verifies the input reference (J:)
# Returns:  the primary key of the reference or 0 if invalid
# Assumes:  nothing
# Effects:  verifies that the Reference exists by checking the referenceDict
#	dictionary for the reference ID or the database.
#	writes to the error file if the Reference is invalid.
#	adds the Reference ID/Key to the global referenceDict dictionary if the
#	reference is valid.
# Throws:

def verifyReference(
    referenceID,          # reference accession ID; J:#### (string)
    lineNum		  # line number (integer)
    ):

    global referenceDict

    if referenceDict.has_key(referenceID):
        referenceKey = referenceDict[referenceID]
    else:
        referenceKey = accessionlib.get_Object_key(referenceID, 'Reference')
        if referenceKey is None:
            errorFile.write('Invalid Reference (%d): %s\n' % (lineNum, referenceID))
            referenceKey = 0
        else:
            referenceDict[referenceID] = referenceKey

    return(referenceKey)

# Purpose:  verify Probe Prep Coverage
# Returns:  Probe Prep Coverage key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Coverage exists in the Coverage dictionary
#	writes to the error file if the Prep Coverage is invalid
# Throws:  nothing

def verifyPrepCoverage(
    coverage, 	# Coverage value (string)
    lineNum	# line number (integer)
    ):

    global coverageDict

    coverageKey = 0

    if len(coverageDict) == 0:
	results = db.sql('select _Coverage_key, coverage from GXD_LabelCoverage', 'auto')
	for r in results:
	    coverageDict[r['coverage']] = r['_Coverage_key']

    if coverageDict.has_key(coverage):
        coverageKey = coverageDict[coverage]
    else:
        errorFile.write('Invalid Prep Coverage (%d): %s\n' % (lineNum, coverage))

    return coverageKey

# Purpose:  verify Probe Prep Label
# Returns:  Probe Prep Label key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Label exists in the Label dictionary
#	writes to the error file if the Prep Label is invalid
# Throws:  nothing

def verifyPrepLabel(
    label, 	# Label value (string)
    lineNum	# line number (integer)
    ):

    global labelDict

    labelKey = 0

    if len(labelDict) == 0:
	results = db.sql('select _Label_key, label from GXD_Label', 'auto')
	for r in results:
	    labelDict[r['label']] = r['_Label_key']

    if labelDict.has_key(label):
        labelKey = labelDict[label]
    else:
        errorFile.write('Invalid Prep Label (%d): %s\n' % (lineNum, label))

    return labelKey

# Purpose:  verify Probe Prep Sense
# Returns:  Probe Prep Sense key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Sense exists in the sense dictionary
#	writes to the error file if the Prep Sense is invalid
# Throws:  nothing

def verifyPrepSense(
    sense, 	# Sense value (string)
    lineNum	# line number (integer)
    ):

    global senseDict

    senseKey = 0

    if len(senseDict) == 0:
	results = db.sql('select _Sense_key, sense from GXD_ProbeSense', 'auto')
	for r in results:
	    senseDict[r['sense']] = r['_Sense_key']

    if senseDict.has_key(sense):
        senseKey = senseDict[sense]
    else:
        errorFile.write('Invalid Prep Sense (%d): %s\n' % (lineNum, sense))

    return senseKey

# Purpose:  verify Probe Prep Type
# Returns:  1 if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Type exists in the probeType list
#	writes to the error file if the Prep Type is invalid
# Throws:  nothing

def verifyPrepType(
    prepType, 	# Type value (string)
    lineNum	# line number (integer)
    ):

    if prepType in prepTypeList:
	return 1
    else:
        errorFile.write('Invalid Prep Type (%d): %s\n' % (lineNum, prepType))
        return 0

# Purpose:  verify Probe Prep Visualization
# Returns:  Probe Prep Visualization key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Visualization exists in the Visualization dictionary
#	writes to the error file if the Prep Visualization is invalid
# Throws:  nothing

def verifyPrepVisualization(
    visualization, 	# Visualization value (string)
    lineNum		# line number (integer)
    ):

    global visualDict

    visualKey = 0

    if len(visualDict) == 0:
	results = db.sql('select _Visualization_key, visualization from GXD_VisualizationMethod', 'auto')
	for r in results:
	    visualDict[r['visualization']] = r['_Visualization_key']

    if visualDict.has_key(visualization):
        visualKey = visualDict[visualization]
    else:
        errorFile.write('Invalid Prep Visualization (%d): %s\n' % (lineNum, visualization))

    return visualKey

# Purpose:  verify Assay Type
# Returns:  Assay Type key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Assay Type exists in the Assay Type dictionary
#	writes to the error file if the Assay Type is invalid
# Throws:  nothing

def verifyAssayType(
    assayType, 	# Assay Type value (string)
    lineNum	# line number (integer)
    ):

    global assayTypeDict

    assayTypeKey = 0

    if len(assayTypeDict) == 0:
	results = db.sql('select _AssayType_key, assayType from GXD_AssayType', 'auto')
	for r in results:
	    assayTypeDict[r['assayType']] = r['_AssayType_key']

    if assayTypeDict.has_key(assayType):
        assayTypeKey = assayTypeDict[assayType]
    else:
        errorFile.write('Invalid Assay Type (%d): %s\n' % (lineNum, assayType))

    return assayTypeKey

# Purpose:  verify Gel RNA Type
# Returns:  Gel RNA Type key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel RNA Type exists in the Gel RNA Type dictionary
#	writes to the error file if the Gel RNA Type is invalid
# Throws:  nothing

def verifyGelRNAType(
    gelRNAType, 	# Gel RNA Type value (string)
    lineNum	# line number (integer)
    ):

    global gelRNATypeDict

    gelRNATypeKey = 0

    if len(gelRNATypeDict) == 0:
	results = db.sql('select _GelRNAType_key, rnaType from GXD_GelRNAType', 'auto')
	for r in results:
	    gelRNATypeDict[r['rnaType']] = r['_GelRNAType_key']

    if gelRNATypeDict.has_key(gelRNAType):
        gelRNATypeKey = gelRNATypeDict[gelRNAType]
    else:
        errorFile.write('Invalid Gel RNA Type (%d): %s\n' % (lineNum, gelRNAType))

    return gelRNATypeKey

# Purpose:  verify Gel Control
# Returns:  Gel Control key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Control exists in the Gel Control dictionary
#	writes to the error file if the Gel Control is invalid
# Throws:  nothing

def verifyGelControl(
    gelControl, 	# Gel Control value (string)
    lineNum	# line number (integer)
    ):

    global gelControlDict

    gelControlKey = 0

    if len(gelControlDict) == 0:
	results = db.sql('select _GelControl_key, gelLaneContent from GXD_GelControl', 'auto')
	for r in results:
	    gelControlDict[r['gelLaneContent']] = r['_GelControl_key']

    if gelControlDict.has_key(gelControl):
        gelControlKey = gelControlDict[gelControl]
    else:
        errorFile.write('Invalid Gel Control (%d): %s\n' % (lineNum, gelControl))

    return gelControlKey

# Purpose:  verifies the genotype
# Returns:  the primary key of the genotype or 0 if invalid
# Assumes:  nothing
# Effects:  verifies that the Genotype exists by checking the genotypeDict
#	dictionary for the Genotype ID or the database.
#	writes to the error file if the Genotype is invalid.
#	adds the Genotype ID/Key to the global genotypeDict dictionary if the
#	genotype is valid.
# Throws:

def verifyGenotype(
    genotypeID,          # genotype accession ID; MGI:#### (string)
    lineNum		 # line number (integer)
    ):

    global genotypeDict

    if genotypeDict.has_key(genotypeID):
        genotypeKey = genotypeDict[genotypeID]
    else:
        genotypeKey = accessionlib.get_Object_key(genotypeID, 'Genotype')
        if genotypeKey is None:
            errorFile.write('Invalid Genotype (%d): %s\n' % (lineNum, genotypeID))
            genotypeKey = 0
        else:
            genotypeDict[genotypeID] = genotypeKey

    return(genotypeKey)

# Purpose:  verifies the anatomical structure
# Returns:  the primary key of the anatomical structure or 0 if invalid
# Assumes:  nothing
# Effects:  verifies that the Anatomical Structure exists by checking the structureDict
#	dictionary for the Structure Name or the database.
#	writes to the error file if the Anatomical Structure is invalid.
#	adds the Structure Name/TS/Key to the global structureDict dictionary if the
#	structure is valid.
# Throws:

def verifyStructure(
    structureName,       # structure name (string)
    theilerStage,	 # theiler stage (integer)
    lineNum		 # line number (integer)
    ):

    global structureDict

    key = '%s:%s' % (structureName, theilerStage)

    if structureDict.has_key(key):
        structureKey = structureDict[key]
    else:
        results = db.sql('select s._Structure_key ' + \
                'from GXD_Structure s, GXD_TheilerStage t ' + \
                'where s._Stage_key = t._Stage_key ' + \
                'and t.stage = %s ' % (str(theilerStage)) + \
                'and s.printName = "%s" ' % (structureName), 'auto')
        if len(results) == 0:
            errorFile.write('Invalid Structure (%d): %s:%d\n' % (lineNum, structureName, theilerStage))
            structureKey = 0
        else:
	    for r in results:
                structureKey = r['_Structure_key']
                structureDict[key] = structureKey

    return(structureKey)

# Purpose:  verify Gel Units
# Returns:  Gel Units key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Units exists in the Gel Units dictionary
#	writes to the error file if the Gel Units is invalid
# Throws:  nothing

def verifyGelUnits(
    gelUnits, 	# Gel Units value (string)
    lineNum	# line number (integer)
    ):

    global gelUnitsDict

    gelUnitsKey = 0

    if len(gelUnitsDict) == 0:
	results = db.sql('select _GelUnits_key, units from GXD_GelUnits', 'auto')
	for r in results:
	    gelUnitsDict[r['units']] = r['_GelUnits_key']

    if gelUnitsDict.has_key(gelUnits):
        gelUnitsKey = gelUnitsDict[gelUnits]
    else:
        errorFile.write('Invalid Gel Units (%d): %s\n' % (lineNum, gelUnits))

    return gelUnitsKey

# Purpose:  verify Gel Strength
# Returns:  Gel Strength key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Strength exists in the Gel Strength dictionary
#	writes to the error file if the Gel Strength is invalid
# Throws:  nothing

def verifyGelStrength(
    gelStrength, 	# Gel Strength value (string)
    lineNum	# line number (integer)
    ):

    global gelStrengthDict

    gelStrengthKey = 0

    if len(gelStrengthDict) == 0:
	results = db.sql('select _Strength_key, strength from GXD_Strength', 'auto')
	for r in results:
	    gelStrengthDict[r['strength']] = r['_Strength_key']

    if gelStrengthDict.has_key(gelStrength):
        gelStrengthKey = gelStrengthDict[gelStrength]
    else:
        errorFile.write('Invalid Gel Strength (%d): %s\n' % (lineNum, gelStrength))

    return gelStrengthKey

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
    truncateDB = 'dump transaction %s with truncate_only' % (db.get_sqlDatabase())

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
	db.sql(truncateDB, None)

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
        tokens = string.split(line[:-1], '\t')

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

	markerKey = verifyMarker(markerID, lineNum)
        referenceKey = verifyReference(jnum, lineNum)

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
	    cdate + TAB + cdate + CRT)

	outMarkerFile.write(str(probeKey) + TAB + \
	    str(markerKey) + TAB + \
	    relationship + TAB + \
	    cdate + TAB + cdate + CRT)

        outRefFile.write(str(refKey) + TAB + str(probeKey) + TAB + str(referenceKey) + TAB + \
	    TAB + '0' + TAB + '0' + TAB + cdate + TAB + cdate + CRT)

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
	    cdate + TAB + cdate + TAB + cdate + CRT)

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

	if verifyPrepType(prepType, lineNum) == 0:
	    error = 1

	senseKey = verifyPrepSense(hybridization, lineNum)
	labelKey = verifyPrepLabel(labelledWith, lineNum)
	coverageKey = verifyPrepCoverage(labelCoverage, lineNum)
	visualizationKey = verifyPrepVisualization(visualization, lineNum)

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
	    cdate + TAB + cdate + CRT)

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

	markerKey = verifyMarker(markerID, lineNum)
        referenceKey = verifyReference(jnum, lineNum)
	assayTypeKey = verifyAssayType(assayType, lineNum)

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
#	    TAB + \
#	    createdBy + TAB + \
#	    createdBy + TAB + \
	    cdate + TAB + cdate + CRT)

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
	    cdate + TAB + cdate + TAB + cdate + CRT)

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

	genotypeKey = verifyGenotype(genotypeID, lineNum)
	rnaTypeKey = verifyGelRNAType(rnaType, lineNum)
	controlKey = verifyGelControl(control, lineNum)
	ageMin, ageMax = agelib.ageMinMax(age)
	structureKey = verifyStructure(structureName, structureTS, lineNum)

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
	    cdate + TAB + cdate + CRT)

	outGelLaneStFile.write(
	    str(gelLaneKey) + TAB + \
	    str(structureKey) + TAB + \
	    cdate + TAB + cdate + CRT)

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

	unitsKey = verifyGelUnits(bandUnits, lineNum)
	strengthKey = verifyGelStrength(bandStrength, lineNum)

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
	      cdate + TAB + cdate + CRT)

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
	    cdate + TAB + cdate + CRT)

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
