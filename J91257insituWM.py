#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: J91257insituWM.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate InSitu Whole Mount input files into input files
#	for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	J91257insituWM.py
#
# Envvars:
#
# Inputs:
#
#       E10.5_WM_Coding_Table.txt, a tab-delimited file in the format:
#               field 1: Gene Symbol
#		field 2: MTF# (ignore)
#		field 3: Image file name
#		field 4-16 : expression
#		field 17 : Image label
#
# Outputs:
#
#       4 tab-delimited files:
#
#	In_Situ_probeprep.txt
#	In_Situ_assay.txt
#	In_Situ_specimen.txt
#	In_Situ_results.txt
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

inInSituFile = ''	# file descriptor
inTissueFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['WMDATADIR']

inInSituFileName = datadir + '/tr6118/E10.5_WM_Coding_Table.txt'
inTissueFileName = datadir + '/tr6118/WM_Tissue_Translation.txt'

prepFileName = datadir + '/In_Situ_probeprep.txt'
assayFileName = datadir + '/In_Situ_assay.txt'
specimenFileName = datadir + '/In_Situ_specimen.txt'
resultsFileName = datadir + '/In_Situ_results.txt'

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'Digoxigenin'
labelCoverage = 'continuously'
visualizedWith = 'Alkaline phosphatase'

# constants for Absent
ABSENT = 'Absent'
NA = 'Not Applicable'

# constants for assay
reference = 'J:91257'
assayType = 'RNA in situ'
createdBy = os.environ['CREATEDBY']

# constants for specimen
specimenLabel = '1'
genotype = 'MGI:2166324'	# C57BL
age = 'embryonic day 10.5'
ageNote = NULL
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Not Applicable'
specimenHybridization = 'whole mount'
specimenNote = NULL

# translation of input file strengths and MGI strengths
strengthTrans = {'x':'Present', 
		 '':'Absent'
		 }

mgiProbe = {}	 # dictionary of Probe Name/MGI Acc ID
mgiMarker = {}	 # dictionary of Marker Symbol/MGI Acc ID
tissueTrans = {} # maps input tissue to MGI tissue and Theiler Stage
tissueNote = {}	 # maps input tissue to MGI Tissue Note
assays = {}	 # maps Marker MGI Acc ID to assayKey

# next available assay key
assayKey = 1

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
    global inInSituFile, inTissueFile
    global prepFile, assayFile, specimenFile, resultsFile
    global mgiProbe, mgiMarker
    global tissueTrans, tissueNote
 
    try:
        inInSituFile = open(inInSituFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName)

    try:
        inTissueFile = open(inTissueFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inTissueFileName)

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

    results = db.sql('select a.accID, p.name from ACC_Accession a, PRB_Probe p ' + \
	'where p.name like "MTF#%" ' + \
	'and p._Probe_key = a._Object_key ' + \
	'and a._MGIType_key = 3 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
    for r in results:
	mgiProbe[r['name']] = r['accID']

    results = db.sql('select a.accID, m.symbol from ACC_Accession a, MRK_Marker m ' + \
	'where m._Organism_key = 1 ' + \
	'and m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key = 1 ' + \
	'and a.prefixPart = "MGI:" ' + \
	'and a.preferred = 1', 'auto')
    for r in results:
	mgiMarker[r['symbol']] = r['accID']

    lineNum = 0
    for line in inTissueFile.readlines():

	lineNum = lineNum + 1

	if lineNum == 1:
	    continue

	tokens = string.split(line[:-1], TAB)
	
	try:
	    inTissue = tokens[0]
	    inAge = tokens[1]
	    theilerStage = tokens[2]
	    gxdTissue1 = tokens[3]
	    gxdTissue2 = tokens[4]
	    inNote = tokens[5]
        
	    key = inTissue

	    value = gxdTissue1 + '|' + theilerStage
	    if not tissueTrans.has_key(key):
	        tissueTrans[key] = []
	    tissueTrans[key].append(value)
    
	    if len(gxdTissue2) > 0:
	        value = gxdTissue2 + '|' + theilerStage
	        if not tissueTrans.has_key(key):
	            tissueTrans[key] = []
	        tissueTrans[key].append(value)

	    if len(inNote) > 0:
	        value = inNote
	        tissueNote[key] = value
        except:
	    pass

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    global assayKey, assays

    # For each line in the input file

    specimenKey = 1
    lineNum = 0
    for line in inInSituFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if lineNum == 1:
	    tissueLabels = tokens[3:16]
	    continue

	# else process an actual data line

        try:
	    mouseGene = tokens[0]
	    mtf = tokens[1]
	    eResults = tokens[3:16]

        except:
            print 'Invalid Line (%d): %s\n' % (lineNum, line)

	if not mgiMarker.has_key(mouseGene):
	    continue

	# create one assay per record

	probeName = 'MTF#' + mtf

	# write the probe prep information

	prepFile.write(str(assayKey) + TAB + \
	    mgiProbe[probeName] + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assayKey) + TAB + \
	    mgiMarker[mouseGene] + TAB + \
	    reference + TAB + \
	    assayType + TAB + \
	    TAB + \
	    TAB + \
	    createdBy + CRT)

        # write the specimen (one for each Assay)

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenLabel + TAB + \
	    genotype + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# one result for each Tissue 

	resultKey = 1

	for i in range(len(tissueLabels)):
    
	    for t in tissueTrans[tissueLabels[i]]:

		#
		# -v = ventral
		# -d = dorsal
		#
		# we only want to store 1 result for each -v, -d pair.
		# so, we ignore the -v for Telencephalon, Diencephalon, Spinal cord
		#
		# Mesencephalon = midbrain
		# Rhombencephalon = hindbrain
		# Mid-hindbrain boundary incorporates both tissues, midbrain and hindbrain
		#
		# we only want to store 1 result for midbrain and 1 result for hindbrain
		# but the data to determine this is encompassed by 5 columns of data.
		# so, we ignore the Mesencephalon and Rhombencephalon columns and only
		# "process" the Mid-hindbrain boundary column.  the logic when we hit
		# Mid-hindbrain boundary takes into account the Mesencephalon and
		# Rhombencephalon columns.
		#

	        if tissueLabels[i] in ['Telencephalon - v', 'Diencephalon - v', 'Spinal cord - v',
		    'Mesencephalon - v', 'Mesencephalon - d', 'Rhombencephalon - v', 'Rhombencephalon - d']:
                    continue

	        # Translate the Tissue into a Tissue and Stage
	        [tissue, theilerStage] = string.split(t, '|')

		## Strength

		strength = strengthTrans[eResults[i]]

	        if tissueLabels[i] in ['Telencephalon - d', 'Diencephalon - d', 'Spinal cord - d']:
		    if eResults[i] == 'x' or eResults[i - 1] == 'x':
			strength = 'Present'
		    else:
			strength = 'Absent'

	        if tissueLabels[i] == 'Mid-hindbrain boundary' and tissue == 'midbrain':
		    if eResults[i] == 'x' or eResults[i - 1] == 'x' or eResults[i - 2] == 'x':
			strength = 'Present'
		    else:
			strength = 'Absent'

	        if tissueLabels[i] == 'Mid-hindbrain boundary' and tissue == 'hindbrain':
		    if eResults[i] == 'x' or eResults[i + 1] == 'x' or eResults[i + 2] == 'x':
			strength = 'Present'
		    else:
			strength = 'Absent'

		## Pattern

		if strength == 'Absent':
		    pattern = 'Not Applicable'
		elif strength == 'Present':
		    pattern = 'Not Specified'

	        if tissueLabels[i] in ['Telencephalon - d', 'Diencephalon - d', 'Spinal cord - d']:
		    if eResults[i] == 'x' and eResults[i - 1] != 'x':
			pattern = 'Regionally restricted'
		    elif eResults[i - 1] == 'x' and eResults[i] != 'x':
			pattern = 'Regionally restricted'

	        if tissueLabels[i] == 'Mid-hindbrain boundary':
		    if eResults[i] == 'x':
		        pattern = 'Regionally restricted'
		    elif tissue == 'midbrain' and eResults[i - 2] == 'x' and eResults[i - 1] != 'x':
		        pattern = 'Regionally restricted'
		    elif tissue == 'midbrain' and eResults[i - 2] != 'x' and eResults[i - 1] == 'x':
		        pattern = 'Regionally restricted'
		    elif tissue == 'hindbrain' and eResults[i + 1] == 'x' and eResults[i + 2] != 'x':
		        pattern = 'Regionally restricted'
		    elif tissue == 'hindbrain' and eResults[i + 1] != 'x' and eResults[i + 2] == 'x':
		        pattern = 'Regionally restricted'

		## Note

	        resultNote = ''

	        if tissueLabels[i] in ['Optic vesicle', 'Lens']:
		    resultNote = tissueNote[tissueLabels[i]]

	        if tissueLabels[i] == 'Telencephalon - d':
		    if eResults[i - 1] == 'x' and eResults[i] != 'x':
			resultNote = 'Expression was detected in ventral telencephalon.'
		    elif eResults[i] == 'x' and eResults[i - 1] != 'x':
			resultNote = 'Expression was detected in dorsal telencephalon.'
		    elif eResults[i] == 'x' and eResults[i - 1] == 'x':
			resultNote = 'Expression was detected in dorsal and ventral telencephalon.'

	        if tissueLabels[i] == 'Diencephalon - d':
		    if eResults[i - 1] == 'x' and eResults[i] != 'x':
			resultNote = 'Expression was detected in ventral diencephalon.'
		    elif eResults[i] == 'x' and eResults[i - 1] != 'x':
			resultNote = 'Expression was detected in dorsal diencephalon.'
		    elif eResults[i] == 'x' and eResults[i - 1] == 'x':
			resultNote = 'Expression was detected in dorsal and ventral diencephalon.'

		## Note

	        if tissueLabels[i] == 'Mid-hindbrain boundary':

		    note1 = ''
		    note2 = ''

		    if eResults[i] == 'x':
			note1 = '  Expression was detected in the mid-hindbrain boundary.'

		    if tissue == 'midbrain':
		        if eResults[i - 2] == 'x' and eResults[i - 1] != 'x':
			    note2 = 'Expression was detected in ventral mesencephalon.'
		        elif eResults[i - 2] != 'x' and eResults[i - 1] == 'x':
			    note2 = 'Expression was detected in dorsal mesencephalon.'
		        elif eResults[i - 2] == 'x' and eResults[i - 1] == 'x':
			    note2 = 'Expression was detected in dorsal and ventral mesencephalon.'
		    elif tissue == 'hindbrain':
		        if eResults[i + 1] == 'x' and eResults[i + 2] != 'x':
			    note2 = 'Expression was detected in ventral rhombencephalon.'
		        elif eResults[i + 1] != 'x' and eResults[i + 2] == 'x':
			    note2 = 'Expression was detected in dorsal rhombencephalon.'
		        elif eResults[i + 1] == 'x' and eResults[i + 2] == 'x':
			    note2 = 'Expression was detected in dorsal and ventral rhombencephalon.'

		    resultNote = note2 + note1

	        if tissueLabels[i] == 'Spinal cord - d':
		    if eResults[i - 1] == 'x' and eResults[i] != 'x':
			resultNote = 'Expression was detected in ventral spinal cord.'
		    elif eResults[i] == 'x' and eResults[i - 1] != 'x':
			resultNote = 'Expression was detected in dorsal spinal cord.'
		    elif eResults[i] == 'x' and eResults[i - 1] == 'x':
			resultNote = 'Expression was detected in dorsal and ventral spinal cord.'

	        resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)

	        resultKey = resultKey + 1

	assays[probeName] = assayKey
	assayKey = assayKey + 1

    # end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)

# $Log$
# Revision 1.2  2004/11/12 20:20:54  lec
# TR 6118
#
# Revision 1.1  2004/11/12 17:39:24  lec
# TR 6118
#
