#!/usr/local/bin/python

# $Header$
# $Name$

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
#	program.py
#	-S = database server
#	-D = database
#	-U = user
#	-P = password file
#	-M = mode
#	-R = reference (J:)
#
# Envvars:
#
# Inputs:
#
#	None
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
import getopt
import db
import mgi_utils
import accessionlib
import agelib
import loadlib
import gxdloadlib

#globals

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
bcpdelim = TAB		# bcp file delimiter

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

datadir = os.environ['ASSAYLOADDATADIR']	# file which contains the data files
createdBy = os.environ['CREATEDBY']

diagFile = ''		# diagnostic file descriptor
errorFile = ''		# error file descriptor

# output files

outIndexFile = ''	# file descriptor
outStagesFile = ''	# file descriptor

indexTable = 'GXD_Index'
stagesTable = 'GXD_Index_Stages'

outIndexFileName = datadir + '/' + indexTable + '.bcp'
outStagesFileName = datadir + '/' + stagesTable + '.bcp'

diagFileName = ''	# diagnostic file name
errorFileName = ''	# error file name
passwordFileName = ''	# password file name

mode = ''		# processing mode (load, preview)
reference = ''		# reference (J:)
referenceKey = ''	# reference key
priorityKey = ''	# priority key

# primary keys

indexKey = 0		# GXD_ProbePrep._ProbePrep_key

# constants
indexComments = 'Age of embryo at noon of plug day not specified in reference.'

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
        '-M mode\n' + \
	'-R reference (J:)\n'

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
    global mode, reference
    global outIndexFile, outStagesFile
    global referenceKey, priorityKey
 
    try:
        optlist, args = getopt.getopt(sys.argv[1:], 'S:D:U:P:M:R:')
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
        elif opt[0] == '-R':
            reference = opt[1]
        else:
            showUsage()

    # User must specify Server, Database, User and Password
    password = string.strip(open(passwordFileName, 'r').readline())
    if server == '' or database == '' or user == '' or password == '' \
	or mode == '' or reference == '':
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

    # Set Log File Descriptor
    db.set_sqlLogFD(diagFile)

    diagFile.write('Start Date/Time: %s\n' % (mgi_utils.date()))
    diagFile.write('Server: %s\n' % (server))
    diagFile.write('Database: %s\n' % (database))
    diagFile.write('User: %s\n' % (user))

    errorFile.write('Start Date/Time: %s\n\n' % (mgi_utils.date()))

    referenceKey = loadlib.verifyReference(reference, 0, errorFile)
    priorityKey = gxdloadlib.verifyIdxPriority(os.environ['IDXPRIORITY'], 0, errorFile)

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

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global indexKey

    results = db.sql('select maxKey = max(index_id) + 1 from GXD_Index', 'auto')
    indexKey = results[0]['maxKey']

# Purpose:  BCPs the data into the database
# Returns:  nothing
# Assumes:  nothing
# Effects:  BCPs the data into the database
# Throws:   nothing

def bcpFiles():

    if DEBUG or not bcpon:
        return

    outIndexFile.close()
    outStagesFile.close()

    bcpI = 'cat %s | bcp %s..' % (passwordFileName, db.get_sqlDatabase())
    bcpII = '-c -t\"%s' % (bcpdelim) + '" -S%s -U%s' % (db.get_sqlServer(), db.get_sqlUser())
    truncateDB = 'dump transaction %s with truncate_only' % (db.get_sqlDatabase())

    bcp1 = '%s%s in %s %s' % (bcpI, indexTable, outIndexFileName, bcpII)
    bcp2 = '%s%s in %s %s' % (bcpI, stagesTable, outStagesFileName, bcpII)

    for bcpCmd in [bcp1, bcp2]:
	diagFile.write('%s\n' % bcpCmd)
	os.system(bcpCmd)
	db.sql(truncateDB, None)

    # update statistics
    db.sql('update statistics %s' % (indexTable), None)
    db.sql('update statistics %s' % (stagesTable), None)

    return

# Purpose:  processes assay data
# Returns:  nothing
# Assumes:  nothing
# Effects:  reads in the appropriate assay data to create the output files
# Throws:   nothing

def processAssay():

    global indexKey
    indexAssay = {}

    cmds = []

    cmds.append('select distinct _Refs_key, _Marker_key ' + \
	'into #indexToAdd ' + \
	'from GXD_Assay ' + \
	'where _Refs_key = %s ' % (referenceKey))

    cmds.append('select * from #indexToAdd')

    cmds.append('select distinct i._Marker_key, a._AssayType_key, s.age, s.hybridization ' + \
	'from #indexToAdd i, GXD_Assay a, GXD_Specimen s ' + \
	'where i._Refs_key = a._Refs_key ' + \
	'and i._Marker_key = a._Marker_key ' + \
	'and a._Assay_key = s._Assay_key ' + \
	'union ' + \
	'select distinct i._Marker_key, a._AssayType_key, s.age, "NA" ' + \
	'from #indexToAdd i, GXD_Assay a, GXD_GelLane s ' + \
	'where i._Refs_key = a._Refs_key ' + \
	'and i._Marker_key = a._Marker_key ' + \
	'and a._Assay_key = s._Assay_key ' + \
	'order by i._Marker_key, a._AssayType_key, s.age')

    results = db.sql(cmds, 'auto')

    for r in results[1]:

	 outIndexFile.write(str(indexKey) + TAB + \
	     str(referenceKey) + TAB + \
	     str(r['_Marker_key']) + TAB + \
	     indexComments + TAB + \
	     createdBy + TAB + createdBy + TAB + \
	     loaddate + TAB + loaddate + CRT)

	 indexAssay[r['_Marker_key']] = indexKey
	 indexKey = indexKey + 1

    for r in results[2]:

	indexKey = indexAssay[r['_Marker_key']]

	if r['_AssayType_key'] == 1 and r['hybridization'] == 'whole mount':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RNA-WM', 0, errorFile)
	elif r['_AssayType_key'] == 1 and r['hybridization'] == 'section':
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RNA-sxn', 0, errorFile)
	elif r['_AssayType_key'] == 5:
	    idxAssayKey = gxdloadlib.verifyIdxAssay('RT-PCR', 0, errorFile)

	if string.find(r['age'], 'embryonic day') >= 0:
	    i = string.find(r['age'], 'embryonic day')
	    stage = r['age'][i + 14:]
	elif string.find(r['age'], 'postnatal') >= 0:
	   stage = 'A'

	idxStageKey = gxdloadlib.verifyIdxStage(stage, 0, errorFile)

        outStagesFile.write(str(indexKey) + TAB + \
            str(idxAssayKey) + TAB + \
            str(idxStageKey) + TAB + \
	    createdBy + TAB + createdBy + TAB + \
            loaddate + TAB + loaddate + CRT)

    return

#
# Main
#

init()
verifyMode()
setPrimaryKeys()
processAssay()
bcpFiles()
exit(0)

# $Log$
# Revision 1.4  2003/09/24 12:29:58  lec
# TR 5154
#
# Revision 1.3  2003/09/22 13:24:00  lec
# TR 5154
#
# Revision 1.2  2003/09/22 11:56:56  lec
# TR 5154
#
# Revision 1.1  2003/09/09 15:29:48  lec
# TR 5125
#
#
