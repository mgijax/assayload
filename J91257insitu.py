#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: J91257insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate InSitu input files into input files
#	for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	J91257insitu.py
#
# Envvars:
#
# Inputs:
#
#       E13.5_In_Situ_Coding_Table.txt, a tab-delimited file in the format:
#               field 1: Gene Symbol
#		field 2: MTF# (ignore)
#		field 3: Image file name
#		field 4: Informativity (ignore)
#		field 5: E-expression
#		field 6: E-specificity
#		field 7: E-CNS
#	        field 8-19: remaining expression
#
#       P0_In_Situ_Coding_Table.txt, a tab-delimited file in the format:
#               field 1: Gene Symbol
#		field 2: MTF# (ignore)
#		field 3: Image file name
#		field 4: Informativity (ignore)
#		field 5: P-expression
#		field 6: P-specificity
#		field 7: E-CNS
#	        field 8-21: remaining expression
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

inInSituFile1 = ''	# file descriptor
inInSituFile2 = ''	# file descriptor
inTissueFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['INSITUDATADIR']

inInSituFileName1 = datadir + '/tr6118/E13.5_In_Situ_Coding_Table.txt'
inInSituFileName2 = datadir + '/tr6118/P0_In_Situ_Coding_Table.txt'
inTissueFileName = datadir + '/tr6118/Tissue_Translation.txt'

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
specimenLabel1 = '1'
specimenLabel2 = '2'
genotype = 'MGI:2166324'	# C57BL
age1 = 'embryonic day 13.5'
age2 = 'postnatal newborn'
ageNote = NULL
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = NULL

# E-expression, P-expression

# translation of input file strengths and MGI strengths
strengthTrans1 = {'1':'Present', 
		 '2':'Ambiguous', 
		 '3':'Absent',
		 ' ':'Absent',
		 NULL:'Absent'
		 }

# translation of input file patterns and MGI patterns
patternTrans1 = {'1':'Regionally restricted', 
                '2':'Regionally restricted', 
		'3':'Not Specified', 
		'4':'Not Specified',
		' ':'Not Specified',
		 NULL:'Not Specified',
		 'Absent':'Not Applicable'
		} 

resultNoteTrans1 = {'1':'Expression was neural specific.', 
		   '2':'Expression was non-neural specific.', 
		   '3':'Expression was detected in neural and non-neural tissue.',
		   '4':NULL,
		   'Not Specified':NULL,
		   'Not Applicable':NULL,
		   NULL:NULL,
		   }

# E-CNS, P-CNS

strengthTrans2 = {'1':'Present', 
		 '2':'Present', 
		 NULL:'Absent'
		 }

patternTrans2 = {'1':'Single cells', 
		 '2':'Homogeneous', 
		 NULL:'Not Applicable'
		}

resultNoteTrans2 = {'Single cells':'Expression was detected in some cells.', 
		   'Homogeneous':'Expression was ubiquitous.',
		   'Not Applicable':NULL,
		   NULL:NULL
		   }

# other tissues

strengthTrans3 = {'1':'Present', 
		 '2':'Present', 
		 NULL:'Absent',
		 ' ':'Absent',
		 '  ':'Absent'
		 }

patternTrans3 = {'1':'Single cells', 
		 '2':'Homogeneous', 
		 NULL:'Not Applicable',
		 ' ':'Not Applicable',
		 '  ':'Not Applicable'
		}

resultNoteTrans3 = {'Single cells':'Expression was detected in some cells.', 
		   'Homogeneous':'Expression was detected in all cells.', 
		   'Not Applicable':NULL,
		   NULL:NULL
		   }

mgiProbe = {}	 # dictionary of Probe Name/MGI Acc ID
mgiMarker = {}	 # dictionary of Marker Symbol/MGI Acc ID
tissueTrans1 = {} # maps input tissue to MGI tissue and Theiler Stage
tissueTrans2 = {} # maps input tissue to MGI tissue and Theiler Stage
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
    global inInSituFile1, inInSituFile2, inTissueFile
    global prepFile, assayFile, specimenFile, resultsFile
    global mgiProbe, mgiMarker
    global tissueTrans1, tissueTrans2, tissueNote
 
    try:
        inInSituFile1 = open(inInSituFileName1, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName1)

    try:
        inInSituFile2 = open(inInSituFileName2, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName2)

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

	tokens = string.split(line[:-1], TAB)
	inTissue = tokens[0]
	inAge = tokens[1]
	theilerStage = tokens[2]
	superStructure = tokens[3]
	gxdTissue1 = tokens[4]
	gxdTissue2 = tokens[5]
	inNote = tokens[6]

	if lineNum == 1:
	    continue

	if inAge == 'E13.5':
	    key = inTissue

	    value = gxdTissue1 + '|' + theilerStage
	    if not tissueTrans1.has_key(key):
	        tissueTrans1[key] = []
	    tissueTrans1[key].append(value)

	    if len(gxdTissue2) > 0:
	        value = gxdTissue2 + '|' + theilerStage
	        if not tissueTrans1.has_key(key):
	            tissueTrans1[key] = []
	        tissueTrans1[key].append(value)
	else:
	    key = inTissue

	    value = gxdTissue1 + '|' + theilerStage
	    if not tissueTrans2.has_key(key):
	        tissueTrans2[key] = []
	    tissueTrans2[key].append(value)

	    if len(gxdTissue2) > 0:
	        value = gxdTissue2 + '|' + theilerStage
	        if not tissueTrans2.has_key(key):
	            tissueTrans2[key] = []
	        tissueTrans2[key].append(value)

	if len(inNote) > 0:
	    value = inNote
	    tissueNote[key] = value

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process1():

    global assayKey, assays

    # For each line in the input file

    specimenKey = 1
    lineNum = 0
    for line in inInSituFile1.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if lineNum == 1:
	    tissueLabels = tokens[7:19]
	    continue

	# else process an actual data line

	if len(tokens) > 5:
	    mouseGene = tokens[0]
	    mtf = tokens[1]
	    eExpression = tokens[4]
	    eSpecificity = tokens[5]
	    eCNS = tokens[6]
	    eResults = tokens[7:19]
        else:
	    mouseGene = tokens[0]
	    mtf = tokens[1]
	    eExpression = tokens[4]
	    eSpecificity = tokens[5]
	    eCNS = ''
	    eResults = ''

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
	    specimenLabel1 + TAB + \
	    genotype + TAB + \
	    age1 + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# one result for each Tissue 

	resultKey = 1

	# E-expression
	t = tissueTrans1['E-expression'][0]
	[tissue, theilerStage] = string.split(t, '|')
	strength = strengthTrans1[eExpression]

	if strength == 'Absent':
	    pattern = patternTrans1[strength]
	else:
	    pattern = patternTrans1[eSpecificity]

	if strength == 'Ambiguous':
	    resultNote = 'Expression was low or absent.  '
        else:
	    resultNote = ''

	resultNote = resultNote + resultNoteTrans1[eSpecificity]

	resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)
	resultKey = resultKey + 1

	# E-CNS
	t = tissueTrans1['E-CNS overall'][0]
	[tissue, theilerStage] = string.split(t, '|')

	strength = strengthTrans2[eCNS]
	pattern = patternTrans2[eCNS]
	resultNote = resultNoteTrans2[pattern]

	resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)
	resultKey = resultKey + 1

	for i in range(len(tissueLabels)):
    
	    for t in tissueTrans1[tissueLabels[i]]:

	        # Translate the Tissue into a Tissue and Stage
	        [tissue, theilerStage] = string.split(t, '|')

		if len(eResults) > 0:
		    strength = strengthTrans3[eResults[i]]
	            pattern = patternTrans3[eResults[i]]
	            resultNote = resultNoteTrans3[pattern]
                else:
		    strength = strengthTrans3[eResults]
	            pattern = patternTrans3[eResults]
	            resultNote = resultNoteTrans3[pattern]

                if tissueLabels[i] in ['E-ear', 'E-OLF']:
	            tNote = tissueNote[tissueLabels[i]]
		    if len(resultNote) > 0:
	                resultNote = tNote + '  ' + resultNote
		    else:
			resultNote = tNote

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

    # end of "for line in inInSituFile1.readlines():"

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process2():

    # For each line in the input file

    specimenKey = 2
    lineNum = 0
    for line in inInSituFile2.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if lineNum == 1:
	    tissueLabels = tokens[7:21]
	    continue

	# else process an actual data line

        if len(tokens) > 5:
	    mouseGene = tokens[0]
	    mtf = tokens[1]
	    pExpression = tokens[4]
	    pSpecificity = tokens[5]
	    pCNS = tokens[6]
	    pResults = tokens[7:21]
        else:
	    mouseGene = tokens[0]
	    mtf = tokens[1]
	    pExpression = tokens[4]
	    pSpecificity = ''
	    pCNS = ''
	    pResults = ''

	if not mgiMarker.has_key(mouseGene):
	    continue

	probeName = 'MTF#' + mtf
	assayKey = assays[probeName]

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenLabel2 + TAB + \
	    genotype + TAB + \
	    age2 + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# one result for each Tissue 

	resultKey = 1

	# P-expression
	t = tissueTrans2['P-expression'][0]
	[tissue, theilerStage] = string.split(t, '|')
	strength = strengthTrans1[pExpression]

	if strength == 'Absent':
	    pattern = patternTrans1[strength]
	else:
	    pattern = patternTrans1[pSpecificity]

	if strength == 'Ambiguous':
	    resultNote = 'Expression was low or absent.  '
        else:
	    resultNote = ''

	resultNote = resultNote + resultNoteTrans1[pSpecificity]

	resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)
	resultKey = resultKey + 1

	# P-CNS
	t = tissueTrans2['P-CNS overall'][0]
	[tissue, theilerStage] = string.split(t, '|')
	strength = strengthTrans2[pCNS]
	pattern = patternTrans2[pCNS]
	resultNote = resultNoteTrans2[pattern]

	resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)
	resultKey = resultKey + 1

	for i in range(len(tissueLabels)):
    
	    for t in tissueTrans2[tissueLabels[i]]:

	        # Translate the Tissue into a Tissue and Stage
	        [tissue, theilerStage] = string.split(t, '|')

		if len(pResults) > 0:
		    if i > len(pResults) - 1:
		        strength = strengthTrans3['']
		        pattern = patternTrans3['']
		        resultNote = resultNoteTrans3[pattern]
		    else:
		        strength = strengthTrans3[pResults[i]]
		        pattern = patternTrans3[pResults[i]]
		        resultNote = resultNoteTrans3[pattern]
		else:
		    strength = strengthTrans3[pResults]
		    pattern = patternTrans3[pResults]
		    resultNote = resultNoteTrans3[pattern]

                if tissueLabels[i] in ['P-ear', 'P-OLF']:
	            tNote = tissueNote[tissueLabels[i]]
		    if len(resultNote) > 0:
	                resultNote = tNote + ' ' + resultNote
		    else:
			resultNote = tNote

	        resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)

	    resultKey = resultKey + 1

    # end of "for line in inInSituFile2.readlines():"

#
# Main
#

init()
process1()
process2()
exit(0)

# $Log$
# Revision 1.10  2004/11/19 18:28:00  lec
# TR 6118
#
# Revision 1.9  2004/11/19 16:50:33  lec
# TR 6118
#
# Revision 1.8  2004/11/19 16:26:32  lec
# TR 6118
#
# Revision 1.7  2004/11/15 13:31:25  lec
# assayload-2-0-0
#
# Revision 1.6  2004/11/15 13:21:13  lec
# assayload-2-0-0
#
# Revision 1.5  2004/11/11 20:37:45  lec
# TR 6118
#
# Revision 1.4  2004/10/14 16:58:10  lec
# TR 6118
#
# Revision 1.3  2004/09/16 13:25:03  lec
# TR 6118
#
# Revision 1.2  2004/09/08 17:09:30  lec
# TR 6118
#
# Revision 1.1  2004/09/08 12:41:15  lec
# TR 6118
#
