#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: J91257probe.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate InSitu Probe input files into input file
#	for the probeload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	J91257probe.py
#
# Envvars:
#
# Inputs:
#
#       tr6118/Probes_Table.txt, a tab-delimited file in the format:
#               field 1: Marker Symbol
#               field 2: Probe Name
#               field 3: RCN5
#               field 4: RCN3
#               field 5: Site 1
#               field 6: Site 2
#               field 7: Size
#               field 8: Author GenBank
#               field 9: Sequence ID
#               field 10: Author 5' primer
#               field 11: Primer 5'
#               field 12: Author 3' primer
#               field 13: Primer 3'
#               field 14: Adaptor Note
#               field 15: Note
#               field 16: Note
#
# Outputs:
#
#       1 tab-delimited file:
#
#	probe.txt
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
import db
import string

# insert into PRB_Source values (52416,63468,316367,1,433,450,315167,316335,null,null,null,'embryonic day 13.5',13.5,13.5,0,1001,1001,getdate(),getdate())

#globals

TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
NULL = ''

inProbeFile = ''	# file descriptor

probeFile = ''		# file descriptor

datadir = os.environ['INSITUDATADIR']

inProbeFileName = datadir + '/tr6118/Probes_Table.txt'

probeFileName = datadir + '/probe.txt'

# constants
segmentType = 'cDNA'
jnum = 'J:91257'
name = 'MTF#%s'
regionCovered = 'Nucleotides %s-%s.'
vectorType = 'Plasmid'
insertSite = '%s/%s'
insertSize = '%skb'
organism = 'mouse, laboratory'
strain = 'C57BL'
tissue = 'brain'
age = 'embryonic day 13.5'
gender = 'Not Specified'
cellLine = 'Not Specified'
relationship = 'E'
fnote = 'Fragment was generated by PCR with primers %s and %s.  '
sequenceLogicalDB = 'Sequence DB:%s|'
createdBy = os.environ['CREATEDBY']

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
    global inProbeFile
    global probeFile
 
    try:
        inProbeFile = open(inProbeFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inProbeFileName)

    try:
        probeFile = open(probeFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % probeFileName)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    results = db.sql('select a.accID, m.symbol from ACC_Accession a, MRK_Marker m ' + \
	'where a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart= "MGI:" ' + \
	'and a.preferred = 1 ' + \
	'and a._Object_key = m._Marker_key', 'auto')
    mgiMarkers = {}
    for r in results:
	mgiMarkers[r['symbol']] = r['accID']

    # For each line in the input file

    lineNum = 0
    for line in inProbeFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if lineNum == 1:
	    continue

	# else process an actual data line

        try:
	    markerSymbol = tokens[0]
	    mtf = tokens[1]
	    rcn5 = tokens[2]
	    rcn3 = tokens[3]
	    site1 = tokens[4]
	    site2 = tokens[5]
	    iSize = tokens[6]
#	    authorGenBank = tokens[7]
	    sequenceID = string.strip(tokens[8])
#	    author5 = tokens[9]
	    primer5 = tokens[10]
#	    author3 = tokens[11]
	    primer3 = tokens[12]
#	    adnote = tokens[13]
#	    note = tokens[14]
#	    note = tokens[15]

        except:
            print 'Invalid Line (%d), Length (%d): %s\n' % (lineNum, len(tokens), line)
	    continue

	try:
	    adnote = tokens[13]
        except:
	    adnote = ''

	if not mgiMarkers.has_key(markerSymbol):
	    print 'Invalid Marker: %s' % (markerSymbol)
	    continue

	# write the probe 

	probeFile.write(name % (mtf) + TAB + \
	    jnum + TAB + \
	    organism + TAB + \
	    strain + TAB + \
	    tissue + TAB + \
	    gender + TAB + \
	    cellLine + TAB + \
	    age + TAB + \
	    vectorType + TAB + \
	    segmentType + TAB)

        if len(rcn5) > 0:
	    probeFile.write(regionCovered % (rcn5, rcn3))
        probeFile.write(TAB)

        if len(site1) > 0:
	    probeFile.write(insertSite % (site1, site2))
        probeFile.write(TAB)

        if len(iSize) > 0:
	    probeFile.write(insertSize % (iSize))
        probeFile.write(TAB)

        probeFile.write(mgiMarkers[markerSymbol] + TAB + \
	    relationship + TAB)

        if len(sequenceID) > 0:
	    probeFile.write(sequenceLogicalDB % (sequenceID))
        probeFile.write(TAB)

	note = ''
        if len(primer5) > 0:
	    note = fnote % (primer5, primer3)

        if len(adnote) > 0:
	    note = note + adnote

	if (len(note) > 0):
	    probeFile.write(note)
        probeFile.write(TAB)

	probeFile.write(createdBy + CRT)

    # end of "for line in inProbeFile.readlines():"

#
# Main
#

init()
process()
exit(0)

# $Log$
# Revision 1.3  2004/11/05 16:00:31  lec
# TR 6118
#
# Revision 1.2  2004/10/14 16:58:10  lec
# TR 6118
#
# Revision 1.1  2004/09/08 12:41:15  lec
# TR 6118
#
