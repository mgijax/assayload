#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: insituload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Assays into InSitu Structures
#
#	GXD_ProbePrep
#	GXD_Assay
#	GXD_Specimen
#	GXD_InSituResults
#	GXD_ISResultStructure
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
#	Probe Prep file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: Probe MGI ID
#		field 3: Probe Prep Type
#		field 4: Hybridization
#		field 5: Labelled With
#		field 6: Label Coverage
#		field 7: Visualized With
#
#	Assay file, a tab-delimited file in the format:
#		field 1: Assay #
#		field 2: MGI Marker Accession ID
#		field 3: Reference (J:#####)
#		field 4: Assay Type
#		field 5: Reporter Gene
#		field 6: Created By
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
#
# Outputs:
#
#       BCP files:
#
#	GXD_ProbePrep.bcp		Probe Prep records
#	GXD_Assay.bcp			Assay records
#	GXD_Specimen.bcp		Specimens
#	GXD_InSituResult.bcp		InSitu Results
#	GXD_ISResultStructure.bcp	InSitu Result Structures
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

datadir = os.environ['INSITUDATADIR']	# file which contains the data files

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor

# input files

inPrepFile = ''           # file descriptor
inAssayFile = ''          # file descriptor
inSpecimenFile = ''        # file descriptor
inResultsFile = ''        # file descriptor

inPrepFileName = datadir + '/In_Situ_probeprep.txt'
inAssayFileName = datadir + '/In_Situ_assay.txt'
inSpecimenFileName = datadir + '/In_Situ_specimen.txt'
inResultsFileName = datadir + '/In_Situ_results.txt'

# output files

outPrepFile = ''	# file descriptor
outAssayFile = ''	# file descriptor
outSpecimenFile = ''	# file descriptor
outResultStFile = ''	# file descriptor
outResultFile = ''	# file descriptor
outAccFile = ''         # file descriptor

probeprepTable = 'GXD_ProbePrep'
assayTable = 'GXD_Assay'
specimenTable = 'GXD_Specimen'
resultTable = 'GXD_InSituResult'
resultStTable = 'GXD_ISResultStructure'
accTable = 'ACC_Accession'

outPrepFileName = datadir + '/' + probeprepTable + '.bcp'
outAssayFileName = datadir + '/' + assayTable + '.bcp'
outSpecimenFileName = datadir + '/' + specimenTable + '.bcp'
outResultFileName = datadir + '/' + resultTable + '.bcp'
outResultStFileName = datadir + '/' + resultStTable + '.bcp'
outAccFileName = datadir + '/' + accTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
passwordFileName = ''	# password file name

mode = ''		# processing mode (load, preview)

# primary keys

prepKey = 0		# GXD_ProbePrep._ProbePrep_key
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

# dictionaries to cache data for quicker lookup

referenceDict = {}      # references
markerDict = {}      	# markers
probeDict = {}		# probes
assayTypeDict = {}	# assay type
senseDict = {}          # probe sense
labelDict = {}          # probe label
coverageDict = {}       # probe coverage
visualDict = {}         # probe visualization
embeddingDict = {}	# embedding method
fixationDict = {}	# fixation
genotypeDict = {}	# genotype
strengthDict = {}	# strength
patternDict = {}	# pattern
structureDict = {}      # anatomical structures
reporterGeneDict = {}	# reporter gene
prepTypeList = ['DNA', 'RNA', 'Not Specified'] 	# lookup of probe prep types
hybridizationList = ['section', 'whole mount', 'section from whole mount']

assayProbePrep = {}	# dictionary of Assay ID/Probe Prep keys
assayAssay= {}		# dictionary of Assay ID/Assay keys
assaySpecimen = {}	# dictionary of Assay ID/Specimen ID and Specimen keys

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
    global outAccFile, outPrepFile, outAssayFile
    global outSpecimenFile, outResultStFile, outResultFile
    global inPrepFile, inAssayFile, inSpecimenFile, inResultsFile
 
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

# Purpose:  verify Probe Accession ID
# Returns:  Probe Key if Probe is valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Probe exists either in the probe dictionary or the database
#	writes to the error file if the Probe is invalid
#	adds the probe id and key to the probe dictionary if the Probe is valid
# Throws:  nothing

def verifyProbe(
    probeID, 	# Accession ID of the Probe (string)
    lineNum	# line number (integer)
    ):

    global probeDict

    probeKey = 0

    if probeDict.has_key(probeID):
        return probeDict[probeID]
    else:
        results = db.sql('select _Object_key from PRB_Acc_View where accID = "%s" ' % (probeID), 'auto')

        for r in results:
            if r['_Object_key'] is None:
                errorFile.write('Invalid Mouse Probe (%d) %s\n' % (lineNum, probeID))
                probeKey = 0
            else:
                probeKey = r['_Object_key']
                probeDict[probeID] = probeKey

    return probeKey

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

# Purpose:  verify Embedding Method
# Returns:  Embedding Method key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Embedding Method exists in the Embedding Method dictionary
#	writes to the error file if the Embedding Method is invalid
# Throws:  nothing

def verifyEmbeddingMethod(
    embedding, 	# Embedding Method value (string)
    lineNum	# line number (integer)
    ):

    global embeddingDict

    embeddingKey = 0

    if len(embeddingDict) == 0:
	results = db.sql('select _Embedding_key, embeddingMethod from GXD_EmbeddingMethod', 'auto')
	for r in results:
	    embeddingDict[r['embeddingMethod']] = r['_Embedding_key']

    if embeddingDict.has_key(embedding):
        embeddingKey = embeddingDict[embedding]
    else:
        errorFile.write('Invalid Embedding Method (%d): %s\n' % (lineNum, embedding))

    return embeddingKey

# Purpose:  verify Fixation Method
# Returns:  Fixation key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Fixation Method exists in the Fixation dictionary
#	writes to the error file if the Fixation Method is invalid
# Throws:  nothing

def verifyFixationMethod(
    fixation, 	# Fixation Method value (string)
    lineNum	# line number (integer)
    ):

    global fixationDict

    fixationKey = 0

    if len(fixationDict) == 0:
	results = db.sql('select _Fixation_key, fixation from GXD_FixationMethod', 'auto')
	for r in results:
	    fixationDict[r['fixation']] = r['_Fixation_key']

    if fixationDict.has_key(fixation):
        fixationKey = fixationDict[fixation]
    else:
        errorFile.write('Invalid Fixation (%d): %s\n' % (lineNum, fixation))

    return fixationKey

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
    theilerStage,	 # theiler stage (string)
    lineNum		 # line number (integer)
    ):

    global structureDict

    key = '%s:%s' % (structureName, theilerStage)

    if structureDict.has_key(key):
        structureKey = structureDict[key]
    else:
	results = db.sql('select s._Structure_key ' + \
		'from GXD_Structure s, GXD_StructureName n, GXD_TheilerStage t ' + \
		'where s._Structure_key = n._Structure_key ' + \
		'and s._Stage_key = t._Stage_key ' + \
		'and t.stage = %s ' % (str(theilerStage)) + \
		'and n.structure = "%s" ' % (structureName) + \
		'union ' + \
		'select s._Structure_key ' + \
		'from GXD_Structure s, GXD_TheilerStage t ' + \
		'where s._Stage_key = t._Stage_key ' + \
		'and t.stage = %s ' % (str(theilerStage)) + \
		'and s.printName = "%s" ' % (structureName), 'auto')
        if len(results) == 0:
            errorFile.write('Invalid Structure (%d): %s:%s\n' % (lineNum, structureName, theilerStage))
            structureKey = 0
        else:
	    for r in results:
                structureKey = r['_Structure_key']
                structureDict[key] = structureKey

    return(structureKey)

# Purpose:  verify Strength
# Returns:  Strength key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Strength exists in the Strength dictionary
#	writes to the error file if the Strength is invalid
# Throws:  nothing

def verifyStrength(
    strength, 	# Strength value (string)
    lineNum	# line number (integer)
    ):

    global strengthDict

    strengthKey = 0

    if len(strengthDict) == 0:
	results = db.sql('select _Strength_key, strength from GXD_Strength', 'auto')
	for r in results:
	    strengthDict[r['strength']] = r['_Strength_key']

    if strengthDict.has_key(strength):
        strengthKey = strengthDict[strength]
    else:
        errorFile.write('Invalid Strength (%d): %s\n' % (lineNum, strength))

    return strengthKey

# Purpose:  verify Pattern
# Returns:  Pattern key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Pattern exists in the Pattern dictionary
#	writes to the error file if the Pattern is invalid
# Throws:  nothing

def verifyPattern(
    pattern, 	# Pattern value (string)
    lineNum	# line number (integer)
    ):

    global patternDict

    patternKey = 0

    if len(patternDict) == 0:
	results = db.sql('select _Pattern_key, pattern from GXD_Pattern', 'auto')
	for r in results:
	    patternDict[r['pattern']] = r['_Pattern_key']

    if patternDict.has_key(pattern):
        patternKey = patternDict[pattern]
    else:
        errorFile.write('Invalid Pattern (%d): %s\n' % (lineNum, pattern))

    return patternKey

# Purpose:  verify Hybridization Type
# Returns:  1 if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Hybridization exists in the hybridization list
#	writes to the error file if the Hybridization is invalid
# Throws:  nothing

def verifyHybridization(
    hybridization, # value (string)
    lineNum	# line number (integer)
    ):

    if hybridization in hybridizationList:
	return 1
    else:
        errorFile.write('Invalid Prep Type (%d): %s\n' % (lineNum, hybridization))
        return 0

# Purpose:  verify Reporter Gene
# Returns:  Reporter Gene key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Reporter Gene exists in the Reporter Gene dictionary
#	writes to the error file if the Reporter Gene is invalid
# Throws:  nothing

def verifyReporterGene(
    reporterGene, # Reporter Gene value (string)
    lineNum	# line number (integer)
    ):

    global reporterGeneDict

    reporterGeneKey = 0

    if len(reporterGeneDict) == 0:
	results = db.sql('select t._Term_key, t.term ' + \
	    'from VOC_Vocab v, VOC_Term t ' + \
	    'where v.name = "GXD Reporter Gene" ' + \
	    'and v._Vocab_key = t._Vocab_key', 'auto')
	for r in results:
	    reporterGeneDict[r['term']] = r['_Term_key']

    if reporterGeneDict.has_key(reporterGene):
        reporterGeneKey = reporterGeneDict[reporterGene]
    else:
        errorFile.write('Invalid Reporter Gene (%d): %s\n' % (lineNum, reporterGene))

    return reporterGeneKey

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global accKey, mgiKey, prepKey, assayKey
    global specimenKey, resultKey

    results = db.sql('select maxKey = max(_ProbePrep_key) + 1 from GXD_ProbePrep', 'auto')
    prepKey = results[0]['maxKey']

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
    outSpecimenFile.close()
    outResultStFile.close()
    outResultFile.close()
    outAccFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"%s' % (bcpdelim) + '" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())
    truncateDB = 'dump transaction %s with truncate_only' % (db.get_sqlDatabase())

    bcp1 = '%s%s in %s %s' % (bcpI, probeprepTable, outPrepFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, assayTable, outAssayFileName, bcpII)
    bcp3 = '%s%s in %s %s' % (bcpI, specimenTable, outSpecimenFileName, bcpII)
    bcp4 = '%s%s in %s %s' % (bcpI, resultStTable, outResultStFileName, bcpII)
    bcp5 = '%s%s in %s %s' % (bcpI, resultTable, outResultFileName, bcpII)
    bcp6 = '%s%s in %s %s' % (bcpI, accTable, outAccFileName, bcpII)

    for bcpCmd in [bcp1, bcp2, bcp3, bcp4, bcp5, bcp6]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)
	db.sql(truncateDB, None)

    # load the cache tables for the records processed (by assay Key)

    for i in range(assayKey - recordsProcessed, assayKey + 1):
	db.sql('exec GXD_loadCacheByAssay %d' % (i), None)

    # update the max Accession ID value
    db.sql('exec ACC_setMax %d' % (recordsProcessed), None)

    # update statistics
    db.sql('update statistics %s' % (probeprepTable), None)
    db.sql('update statistics %s' % (assayTable), None)
    db.sql('update statistics %s' % (specimenTable), None)
    db.sql('update statistics %s' % (resultStTable), None)
    db.sql('update statistics %s' % (resultTable), None)
    db.sql('update statistics %s' % (accTable), None)

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
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    probeID = tokens[1]
	    prepType = tokens[2]
	    hybridization = tokens[3]
	    labelledWith = tokens[4]
	    labelCoverage = tokens[5]
	    visualization = tokens[6]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	if verifyPrepType(prepType, lineNum) == 0:
	    error = 1

	probeKey = verifyProbe(probeID, lineNum)
	senseKey = verifyPrepSense(hybridization, lineNum)
	labelKey = verifyPrepLabel(labelledWith, lineNum)
	coverageKey = verifyPrepCoverage(labelCoverage, lineNum)
	visualizationKey = verifyPrepVisualization(visualization, lineNum)

        if probeKey == 0 or senseKey == 0 or labelKey == 0 or coverageKey == 0:
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
	    reporterGene = tokens[4]
	    createdBy = tokens[5]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	markerKey = verifyMarker(markerID, lineNum)
        referenceKey = verifyReference(jnum, lineNum)
	assayTypeKey = verifyAssayType(assayType, lineNum)

        if markerKey == 0 or referenceKey == 0 or assayTypeKey == 0:
            # set error flag to true
            error = 1

        if len(reporterGene) > 0:
            reporterGeneKey = verifyReporterGene(reporterGene, lineNum)
	    if reporterGeneKey == 0:
                error = 1
        else:
            reporterGeneKey = ''

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
#            str(reporterGeneKey) + TAB + \
#            createdBy + TAB + \
#            createdBy + TAB + \
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
        tokens = string.split(line[:-1], '\t')

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

	if verifyHybridization(hybridization, lineNum) == 0:
	    error = 1

	genotypeKey = verifyGenotype(genotypeID, lineNum)
	fixationKey = verifyFixationMethod(fixation, lineNum)
	embeddingKey = verifyEmbeddingMethod(embedding, lineNum)
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
	    cdate + TAB + cdate + CRT)

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
    prevResult = 0
    lineNum = 0
    # For each line in the input file

    for line in inResultsFile.readlines():

        error = 0
        lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], '\t')

        try:
	    assayID = tokens[0]
	    specimenID = tokens[1]
	    resultID = tokens[2]
	    strength = tokens[3]
	    pattern = tokens[4]
	    structureName = tokens[5]
	    structureTS = tokens[6]
	    resultNote = tokens[7]
        except:
            exit(1, 'Invalid Line (%d): %s\n' % (lineNum, line))

	strengthKey = verifyStrength(strength, lineNum)
	patternKey = verifyPattern(pattern, lineNum)
	structureKey = verifyStructure(structureName, structureTS, lineNum)

        if strengthKey == 0 or patternKey == 0 or structureKey == 0:
            # set error flag to true
            error = 1

        # if errors, continue to next record
        if error:
            continue

        # if no errors, process

	key = '%s:%s' % (assayID, specimenID)
	specimenKey = assaySpecimen[key]

	if prevAssay != assayID:
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
	      cdate + TAB + cdate + CRT)

	outResultStFile.write(
	    str(resultKey) + TAB + \
	    str(structureKey) + TAB + \
	    cdate + TAB + cdate + CRT)

	prevAssay = assayID
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

# $Log$
# Revision 1.5  2003/06/20 15:16:11  lec
# TR 4800
#
# Revision 1.4  2003/06/18 18:21:50  lec
# TR 4800
#
# Revision 1.3  2003/06/18 17:57:14  lec
# TR 4800
#
# Revision 1.2  2003/06/18 15:56:16  lec
# TR 4800
#
# Revision 1.1  2003/06/18 13:19:23  lec
# new
#
#
