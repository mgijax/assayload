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

#globals

DEBUG = 0		# if 0, not in debug mode
TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
bcpdelim = TAB		# bcp file delimiter

bcpon = 1		# can the bcp files be bcp-ed into the database?  default is yes.

datadir = os.environ['ASSAYLOADDATADIR']	# file which contains the data files

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

# primary keys

indexKey = 0		# GXD_ProbePrep._ProbePrep_key

# constants
indexComments = 'Age of embryo at noon of plug day not specified in reference.'

stageDict = {'0.5' : 0, \
    '1.0' : 1, \
    '1.5' : 2, \
    '2.0' : 3, \
    '2.5' : 4, \
    '3.0' : 5, \
    '3.5' : 6, \
    '4.0' : 7, \
    '4.5' : 8, \
    '5.0' : 9, \
    '5.5' : 10, \
    '6.0' : 11, \
    '6.5' : 12, \
    '7.0' : 13, \
    '7.5' : 14, \
    '8.0' : 15, \
    '8.5' : 16, \
    '9.0' : 17, \
    '9.5' : 18, \
    '10.0' : 19, \
    '10.5' : 20, \
    '11.0' : 21, \
    '11.5' : 22, \
    '12.0' : 23, \
    '12.5' : 24, \
    '13.0' : 25, \
    '13.5' : 26, \
    '14.0' : 27, \
    '14.5' : 28, \
    '15.0' : 29, \
    '15.5' : 30, \
    '16.0' : 31, \
    '16.5' : 32, \
    '17.0' : 33, \
    '17.5' : 34, \
    '18.0' : 35, \
    '18.5' : 36, \
    '19.0' : 37, \
    '19.5' : 38, \
    '20' : 39, \
    'E' : 40, \
    'A' : 41}

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


# Purpose:  verifies the input reference (J:)
# Returns:  nothing
# Assumes:  nothing
# Effects:  if the reference is not valid, exits.
#	else, sets global variable referenceKey.
# Throws: nothing

def verifyReference():

    global referenceKey

    referenceKey = accessionlib.get_Object_key(reference, 'Reference')
    if referenceKey is None:
        exit(1, 'Invalid Reference:  %s\n' % (reference))

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
	     cdate + TAB + cdate + CRT)

	 indexAssay[r['_Marker_key']] = indexKey
	 indexKey = indexKey + 1

    prevStage = ''
    insitu_protein_section = 0
    insitu_rna_section = 0
    insitu_protein_mount = 0
    insitu_rna_mount = 0
    northern = 0
    western = 0
    rt_pcr = 0
    clones = 0
    rnase = 0
    nuclease = 0
    primer_extension = 0

    for r in results[2]:

	if prevStage != r['age']:

	   if prevStage != '':
               outStagesFile.write(str(indexKey) + TAB + \
	           str(stageKey) + TAB + \
	           str(insitu_protein_section) + TAB + \
	           str(insitu_rna_section) + TAB + \
	           str(insitu_protein_mount) + TAB + \
	           str(insitu_rna_mount) + TAB + \
	           str(northern) + TAB + \
	           str(western) + TAB + \
	           str(rt_pcr) + TAB + \
	           str(clones) + TAB + \
	           str(rnase) + TAB + \
	           str(nuclease) + TAB + \
	           str(primer_extension) + TAB + \
	           cdate + TAB + cdate + CRT)

	   insitu_protein_section = 0
	   insitu_rna_section = 0
	   insitu_protein_mount = 0
	   insitu_rna_mount = 0
	   northern = 0
	   western = 0
	   rt_pcr = 0
	   clones = 0
	   rnase = 0
	   nuclease = 0
	   primer_extension = 0

	if r['_AssayType_key'] == 1 and r['hybridization'] == 'whole mount':
	    insitu_rna_mount = 1
	elif r['_AssayType_key'] == 1 and r['hybridization'] == 'section':
	    insitu_rna_section = 1
	elif r['_AssayType_key'] == 5:
	    rt_pcr = 1

	indexKey = indexAssay[r['_Marker_key']]

	if string.find(r['age'], 'embryonic day') >= 0:
	    i = string.find(r['age'], 'embryonic day')
	    stage = r['age'][i + 14:]
	elif string.find(r['age'], 'postnatal') >= 0:
	   stage = 'A'

	stageKey = stageDict[stage]
	prevStage = r['age']

    outStagesFile.write(str(indexKey) + TAB + \
        str(stageKey) + TAB + \
        str(insitu_protein_section) + TAB + \
        str(insitu_rna_section) + TAB + \
        str(insitu_protein_mount) + TAB + \
        str(insitu_rna_mount) + TAB + \
        str(northern) + TAB + \
        str(western) + TAB + \
        str(rt_pcr) + TAB + \
        str(clones) + TAB + \
        str(rnase) + TAB + \
        str(nuclease) + TAB + \
        str(primer_extension) + TAB + \
        cdate + TAB + cdate + CRT)

    return

#
# Main
#

init()
verifyMode()
verifyReference()
setPrimaryKeys()
processAssay()
bcpFiles()
exit(0)

# $Log$
# Revision 1.1  2003/09/09 15:29:48  lec
# TR 5125
#
#
