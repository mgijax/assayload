#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: J80501rnainsitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RNA InSitu 10 and 14 input files into input files
#	for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	J80501rnainsitu.py
#
# Envvars:
#
# Inputs:
#
#       9.5_dpc_in_situ_results.txt, a tab-delimited file in the format:
#               field 1: Paper ID
#		field 2: Clone ID
#		field 3: MGI clone ID
#		field 4: MGI gene symbol
#		field 5: MGI gene ID
#		field 6: overall expression
#		field 7-24: expression
#		field 25-27: image file names
#
#	9.5dpc_in_situ_tissue.txt, a tab-delimited file in the format:
#		field 1: Reported Tissue
#		field 2: MGI Tissue
#		field 3: Theiler Stage
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

datadir = os.environ['INSITUDATADIR']

inInSituFileName = datadir + '/tr5154/9.5_dpc_in_situ_results.txt'
inTissueFileName = datadir + '/tr5154/9.5dpc_in_situ_tissue.txt'

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
reference = 'J:80501'
assayType = 'RNA In Situ'
createdBy = os.environ['CREATEDBY']

# constants for specimen
specimenLabel = '%s 9.5dpc'
genotype = 'MGI:2166653'	# NMRI
age = 'embryonic day 9.5'
ageNote = 'Staging was performed using morphological criteria.'
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Not Applicable'
specimenHybridization = 'whole mount'
specimenNote = NULL

# translation of input file strengths and MGI strengths
strengthTrans = {'(+)':'Present', '(+)rr':'Present', NULL:'Absent', 'ND':'Absent'}
# translation of input file patterns and MGI patterns
patternTrans = {'(+)':'Not Specified', '(+)rr':'Regionally restricted', NULL:'Not Specified', 'ND':'Not Specified'}
resultNoteTrans = {'Wh_emb':'The expression was widespread.'}

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

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    global assayKey

    tissueTrans = {}	# maps input tissue to MGI tissue and Theiler Stage

    lineNum = 0
    for line in inTissueFile.readlines():
	lineNum = lineNum + 1

	tokens = string.split(line[:-1], TAB)
	badTissue = tokens[0]
	goodTissue = tokens[1]
	theilerStage = tokens[2]

	if lineNum == 1:
	    continue

	key = badTissue
	value = goodTissue + '|' + theilerStage
	if not tissueTrans.has_key(key):
	    tissueTrans[key] = []
	tissueTrans[key].append(value)

    # For each line in the input file

    lineNum = 0
    for line in inInSituFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if lineNum == 1:
	    tissueLabels = tokens[5:24]
	    continue

	# else process an actual data line

        try:
	    paperID = tokens[0]
	    cloneID = tokens[1]
	    mgiCloneID = tokens[2]
	    mouseGene = tokens[3]
	    accID = tokens[4]
	    results = tokens[5:24]
	    imageFileName = tokens[25:27]

        except:
            print 'Invalid Line (%d): %s\n' % (lineNum, line)

	if len(mouseGene) == 0:
	    continue

	# create one assay per record

	# write the probe prep information

	prepFile.write(str(assayKey) + TAB + \
	    mgiCloneID + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assayKey) + TAB + \
	    accID + TAB + \
	    reference + TAB + \
	    assayType + TAB + \
	    TAB + \
	    createdBy + CRT)

        # write the specimen (one for each Assay)

	specimenKey = 1

	specimenFile.write(str(assayKey) + TAB + \
	    str(specimenKey) + TAB + \
	    specimenLabel % (mouseGene) + TAB + \
	    genotype + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    sex + TAB + \
	    fixation + TAB + \
	    embedding + TAB + \
	    specimenHybridization + TAB + \
	    specimenNote + CRT)

	# one result for each Tissue 

	ignoreRemainingTissues = 0
	resultKey = 1
	for i in range(len(tissueLabels)):
    
	    for t in tissueTrans[tissueLabels[i]]:

	        # Translate the Tissue into a Tissue and Age
	        [tissue, theilerStage] = string.split(t, '|')
	        resultNote = NULL

	        # if Wh_emb = (+) and other tissues have no results, then ignore other tissues
	        # if Wh_emb = (+) and other tissues have results, then process other tissues
	        # if Wh_emb = ND, then ignore other tissues
	        # if Wh_emb = null, then process other tissues

	        if tissueLabels[i] == 'Wh_emb':

		    if results[i] == '':
		        continue	# skip processing of Wh_emb

		    # deterimne if remaining tissues have results

		    ignoreRemainingTissues = 1

	            if results[i] == '(+)':
	                resultNote = resultNoteTrans[tissueLabels[i]]
		        for j in range(2, len(tissueLabels)):
	                    if results[j] != '':
		                ignoreRemainingTissues = 0

	        # process tissue

	        strength = strengthTrans[results[i]]
	        pattern = patternTrans[results[i]]

	        resultsFile.write(str(assayKey) + TAB + \
	            str(specimenKey) + TAB + \
	            str(resultKey) + TAB + \
		    strength + TAB + \
		    pattern + TAB + \
		    tissue + TAB + \
		    theilerStage + TAB + \
		    resultNote + CRT)

	    resultKey = resultKey + 1

	    if ignoreRemainingTissues:
		break

	specimenKey = specimenKey + 1
	assayKey = assayKey + 1

    # end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)

# $Log$
# Revision 1.4  2003/09/22 13:00:30  lec
# TR 5154
#
# Revision 1.3  2003/09/19 19:32:58  lec
# TR 5154
#
# Revision 1.2  2003/09/19 18:42:55  lec
# TR 5154
#
# Revision 1.1  2003/09/19 18:09:01  lec
# TR 5154
#
#
