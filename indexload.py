#!/usr/local/bin/python

#
# Program: indexload.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Assays into Index Structures
#
#	GXD_Index
#	GXD_Index_Stages
#
# Requirements Satisfied by This Program:
#
# Usage:
#	indexload.py
#
# Envvars:
#
#	LOADFILE5	marker notes; leave blank if not individual notes are needed
#
# Inputs:
#
#       LoadFile5.txt, a tab-delimited file in the format:
#
#               field 0: J number
#		field 1: Marker Symbol
#		field 2: Marker MGI ID
#		field 3: Comment/Note
#
# Outputs:
#
#       BCP files:
#
#	GXD_Index.bcp		Index records
#	GXD_Index_Stages.bcp	Stage records
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# Exit Codes:
#
# Assumes:
#
#	That no one else is adding Index records to the database.
#
# Bugs:
#
# Implementation:
#

import sys
import os
import string
import db
import mgi_utils
import loadlib
import gxdloadlib

#
# from configuration file
#
mode = os.environ['ASSAYLOADMODE']
createdBy = os.environ['CREATEDBY']
reference = os.environ['REFERENCE']
indexpriority = os.environ['IDXPRIORITY']
indexComments = os.environ['IDXCOMMENTS']

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor

# input files

inCommentsFile = ''	# file description

# output files

outIndexFile = ''	# file descriptor
outStagesFile = ''	# file descriptor

indexTable = 'GXD_Index'
stagesTable = 'GXD_Index_Stages'

outIndexFileName = indexTable + '.bcp'
outStagesFileName = stagesTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name

referenceKey = ''	# reference key
priorityKey = ''	# priority key
conditionalKey = 4834242 # conditional key defaults to "not applicable"
createdByKey = ''	# created by key

# comments lookup
commentsLookup = {}

# constants

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
    global inCommentsFile, outIndexFile, outStagesFile
    global referenceKey, priorityKey, createdByKey, indexComments
 
    db.useOneConnection(1)
 
    diagFileName = 'indexload.diagnostics'
    errorFileName = 'indexload.error'

    try:
        diagFile = open(diagFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % diagFileName)
		
    try:
        errorFile = open(errorFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % errorFileName)
		
    # Input Files

    #try:
#	inCommentsFileName = os.environ['LOADFILE5']
#        inCommentsFile = open(inCommentsFileName, 'r')
#    except:
#	pass
        #exit(1, 'Could not open file %s\n' % inCommentsFileName)

    # Output Files

    try:
        outIndexFile = open(outIndexFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outIndexFileName)

    try:
        outStagesFile = open(outStagesFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outStagesFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (db.get_sqlServer()))
    diagFile.write('Database: %s\n' % (db.get_sqlServer()))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    referenceKey = loadlib.verifyReference(reference, 0, errorFile)
    priorityKey = gxdloadlib.verifyIdxPriority(indexpriority, 0, errorFile)
    createdByKey = loadlib.verifyUser(createdBy, 0, errorFile)

    # comments lookup
    try:
    	for line in inCommentsFile.readlines():
		tokens = string.split(line[:-1], TAB)
		markerID = tokens[2]
		comments = tokens[3]
		commentsLookup[markerID] = []
		commentsLookup[markerID] = comments
    	inCommentsFile.close()
    except:
	pass

    return

# Purpose: verify processing mode
# Returns: nothing
# Assumes: nothing
# Effects: if the processing mode is not valid, exits.
#	   else, sets global variables DEBUG and bcpon
# Throws:  nothing

def verifyMode():

    global DEBUG

    if mode == 'preview':
        DEBUG = 1
        bcpon = 0
    elif mode != 'load':
        exit(1, 'Invalid Processing Mode:  %s\n' % (mode))

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():

    outIndexFile.close()
    outStagesFile.close()

    db.commit()
    db.useOneConnection(0)

    if DEBUG or not bcpon:
        return

    bcpCommand = os.environ['PG_DBUTILS'] + '/bin/bcpin.csh'
    currentDir = os.getcwd()

    bcp1 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
        % (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), indexTable, currentDir, outIndexFileName)

    bcp2 =  '%s %s %s %s %s %s "\\t" "\\n" mgd' \
        % (bcpCommand, db.get_sqlServer(), db.get_sqlDatabase(), stagesTable, currentDir, outStagesFileName)

    for bcpCmd in [bcp1, bcp2]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)

    return

# Purpose:  processes assay data
# Returns:  nothing
# Assumes:  nothing
# Effects:  reads in the appropriate assay data to create the output files
# Throws:   nothing

def processAssay():

    # currently existing indexes

    db.sql('''
	create temporary table indexExist on commit drop as
        select distinct a._Refs_key, a._Marker_key, i._Index_key, s._IndexAssay_key, s._StageID_key
	from GXD_Assay a, GXD_Index i, GXD_Index_Stages s
	where a._Refs_key = %s 
	and a._Refs_key = i._Refs_key
	and a._Marker_key = i._Marker_key
	and i._Index_key = s._Index_key''' % (referenceKey), None)

    results = db.sql('select * from indexExist', 'auto')

    # new indexes

    db.sql('''
	create temporary table indexToAdd on commit drop as
        select distinct a._Refs_key, a._Marker_key, aa.accID
	from GXD_Assay a, ACC_Accession aa
	where a._Refs_key = %s
	and a._Marker_key = aa._Object_key
	and aa._MGIType_key = 2
	and aa._LogicalDB_key = 1
	and aa.prefixPart = 'MGI:'
	and aa.preferred = 1
	''' % (referenceKey), None)

    # store new assyas
    results = db.sql('select max(_Index_key) + 1 as maxKey from GXD_Index', 'auto')
    indexKey = results[0]['maxKey']

    indexAssay = {}
    results = db.sql('select distinct _Marker_key from indexToAdd', 'auto')
    for r in results:
	indexAssay[r['_Marker_key']] = indexKey
	indexKey = indexKey + 1

    # store current assyas
    results = db.sql('select _Marker_key, _Index_key from indexExist', 'auto')
    for r in results:
	indexAssay[r['_Marker_key']] = r['_Index_key']

    # store current stages
    indexedAlready = []
    results = db.sql('select _Index_key, _IndexAssay_key, _StageID_key from indexExist', 'auto')
    for r in results:
	indexedTuple = (r['_Index_key'], r['_IndexAssay_key'], r['_StageID_key'])
	indexedAlready.append(indexedTuple)

    # select new indexes (those that do NOT exist)

    results = db.sql('''select a.* from indexToAdd a
	where not exists 
	(select 1 from indexExist e 
	  where a._Refs_key = e._Refs_key
	  and a._Marker_key = e._Marker_key)
	''', 'auto')

    for r in results:

	indexKey = indexAssay[r['_Marker_key']]

	if commentsLookup.has_key(r['accID']):
	    writeIndexComments = commentsLookup[r['accID']]
	else:
            writeIndexComments = indexComments

	outIndexFile.write(str(indexKey) + TAB + \
	     str(referenceKey) + TAB + \
	     str(r['_Marker_key']) + TAB + \
	     str(priorityKey) + TAB + \
	     str(conditionalKey) + TAB + \
	     str(writeIndexComments) + TAB + \
	     str(createdByKey) + TAB + str(createdByKey) + TAB + \
	     loaddate + TAB + loaddate + CRT)

    #
    # select stages
    #

    results = db.sql('''(select distinct i._Marker_key, a._AssayType_key, s.age, s.hybridization 
	from indexToAdd i, GXD_Assay a, GXD_Specimen s 
	where i._Refs_key = a._Refs_key 
	and i._Marker_key = a._Marker_key 
	and a._Assay_key = s._Assay_key 
	union 
	select distinct i._Marker_key, a._AssayType_key, s.age, 'NA' as hybridization
	from indexToAdd i, GXD_Assay a, GXD_GelLane s 
	where i._Refs_key = a._Refs_key 
	and i._Marker_key = a._Marker_key 
	and a._Assay_key = s._Assay_key 
        )
	order by _Marker_key, _AssayType_key, age''', 'auto')

    for r in results:

	indexKey = indexAssay[r['_Marker_key']]

	if r['_AssayType_key'] == 1 and r['hybridization'] == 'whole mount':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RNA-WM', 0, errorFile)
	elif r['_AssayType_key'] == 1 and r['hybridization'] == 'section':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RNA-sxn', 0, errorFile)
	elif r['_AssayType_key'] == 2:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Northern', 0, errorFile)
	elif r['_AssayType_key'] == 3:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('S1 nuc', 0, errorFile)
	elif r['_AssayType_key'] == 4:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RNAse prot', 0, errorFile)
	elif r['_AssayType_key'] == 5:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RT-PCR', 0, errorFile)
	elif r['_AssayType_key'] == 6 and r['hybridization'] == 'whole mount':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Prot-WM', 0, errorFile)
	elif r['_AssayType_key'] == 6 and r['hybridization'] == 'section':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Prot-sxn', 0, errorFile)
	elif r['_AssayType_key'] == 7 and r['hybridization'] == 'whole mount':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Prot-WM', 0, errorFile)
	elif r['_AssayType_key'] == 7 and r['hybridization'] == 'section':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Prot-sxn', 0, errorFile)
	elif r['_AssayType_key'] == 8:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Western', 0, errorFile)
	elif r['_AssayType_key'] == 9:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('Knock in', 0, errorFile)

	#
	# if "stages" == 'embryonic day 15-16', then
	#	stages[0] = 'embryonic day 15'
	#	stages[1] = 'embryonic day 16'
	#

        stages = []

	if string.find(r['age'], 'embryonic day') >= 0:
	    i = string.find(r['age'], 'embryonic day')
	    # all embryonic day stages sorted by '-'
	    allstages = string.split(r['age'][i + 14:], '-')
	    for i in allstages:
		i = string.replace(i, '7.25', '7.5')
		i = string.replace(i, '7.75', '8')
		i = string.replace(i, '10.25', '10.5')
		i = string.replace(i, '10.75', '11.0')
		i = string.replace(i, '11.25', '11.5')
	        stages.append(i)
	elif string.find(r['age'], 'postnatal') >= 0:
	    stages.append('A')

	for s in stages:

	    # the age may have a '.0'...the vocbaulary does not
	    # modify the age so that it can be found in the vocabulary
	    s = string.replace(s, '.0', '')

	    idxStageKey = gxdloadlib.verifyIdxStage(s, 0, errorFile)

	    indexedTuple = (indexKey, idxAssayKey, idxStageKey)
	    if indexedTuple in indexedAlready:
	        continue

            outStagesFile.write(str(indexKey) + TAB + \
                str(idxAssayKey) + TAB + \
                str(idxStageKey) + TAB + \
	        str(createdByKey) + TAB + str(createdByKey) + TAB + \
                loaddate + TAB + loaddate + CRT)

	    indexedAlready.append(indexedTuple)

    return

#
# Main
#

init()
verifyMode()
processAssay()
bcpFiles()
exit(0)

