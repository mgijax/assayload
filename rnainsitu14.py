#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: rnainsitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RNA InSitu input files into input files
#	for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	rnainsitu.py
#
# Envvars:
#
# Inputs:
#
#       In_Situ.txt, a tab-delimited file in the format:
#		field 1: Vial #
#		field 2: Mouse Gene Symbol
#               field 3: MGI Marker Accession ID
#		field 4: Human Gene
#               field 5-16: Tissues
#
#	In_Situ_Tissues.txt, a tab-delimited file in the format:
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

inInSituFile = ''	# file descriptor
inTissueFile = ''	# file descriptor
prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['INSITU14DATADIR']

inInSituFileName = datadir + '/tr4800/14.5_In_situ.txt'
inTissueFileName = datadir + '/tr4800/14.5_In_Situ_tissues.txt'
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

# constants for assay
reference = 'J:80502'
assayType = 'RNA In Situ'

# constants for specimen
specimenLabel = '%s 14.5dpc'
genotype = 'MGI:2166653'
age = 'embryonic day 14.5'
ageNote = 'Age of embryo at noon of plug day not specified in reference.'
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = ''
resultNote = ''

# translation of patterns

patternTrans = {'U':'Homogeneous', 'R':'Regionally restricted', 'I':'Regionally restricted', \
   'U,R':'Not Specified', 'R,U':'Not Specified', 'i':'Regionally restricted'}

patternNA = 'Not Applicable'

# translation of input file strengths and MGI strengths

strengthTrans = {'*':'Weak', '**':'Moderate', '***':'Strong', '':'Absent', \
    'U,R':'Present', 'R,U':'Present', '-':'Not Applicable'}
presentStrength = ['U,R', 'R,U']

mgiTypeKey = 8		# Assay
mgiPrefix = "MGI:"

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
    global inInSituFile, inTissueFile, prepFile, assayFile, specimenFile, resultsFile
 
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

    tissueTrans = {}	# maps input tissue to MGI tissue and Theiler Stage

    for line in inTissueFile.readlines():
	tokens = string.split(line[:-1], TAB)
	badTissue = tokens[0]
	goodTissue = tokens[1]
	theilerStage = tokens[2]

	key = badTissue
	value = goodTissue + '|' + theilerStage
	tissueTrans[key] = value

    assay = 0	# unique Assay ID

    # For each line in the input file

    for line in inInSituFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	if assay == 0:
	    tissueLabels = tokens[7:-1]
	    assay = assay + 1
	    continue

	# else process an actual data line

        try:
	    humanGene = tokens[0]
	    mouseGene = tokens[1]
	    accID = tokens[2]
	    ishNumber = tokens[3]
	    specimen = tokens[4]
	    tissueQuality = tokens[5]
	    overallExpression = tokens[6]
	    results = tokens[7:50]
	    probeID = 'MGI:35046'
#	    probeID = tokens[]

        except:
            print 'Invalid Line (%d): %s\n' % (assay, line)

	if len(mouseGene) == 0:
	    continue

	# write the probe prep information

	prepFile.write(str(assay) + TAB + \
	    probeID + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assay) + TAB + \
	    accID + TAB + \
	    reference + TAB + \
	    assayType + CRT)

	# write the specimen (one for each Assay)

	specimen = 1

	specimenFile.write(str(assay) + TAB + \
		str(specimen) + TAB + \
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

	result = 1
	for i in range(len(tissueLabels)):

	    # Translate the Tissue into a Tissue and Age
	    [tissue, theilerStage] = string.split(tissueTrans[tissueLabels[i]], '|')

	    if len(results) < i + 1:
                # not every tissue was assayed
#		print 'No Results for %s:%s' % (mouseGene, tissueLabels[i])
		continue

	    if len(results[i]) > 1:
	        [inStrength, inPattern] = string.split(results[i], ' ')
		if inPattern in presentStrength:
		  strength = strengthTrans[inPattern]
		  pattern = patternTrans[inPattern]
		else:
	          strength = strengthTrans[inStrength]
	          pattern = patternTrans[inPattern]

	    # no strength or pattern given

	    else:
	        strength = strengthTrans[results[i]]
	        pattern = patternNA

	    resultsFile.write(str(assay) + TAB + \
	        str(specimen) + TAB + \
	        str(result) + TAB + \
		strength + TAB + \
		pattern + TAB + \
		tissue + TAB + \
		theilerStage + TAB + \
		resultNote + CRT)

	    result = result + 1

	specimen = specimen + 1
	assay = assay + 1

    #	end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)

# $Log$
# Revision 1.2  2003/06/18 15:38:25  lec
# new
#
# Revision 1.1  2003/06/18 13:34:23  lec
# new
#
# Revision 1.2  2003/06/18 13:19:33  lec
# TR 4800
#
# Revision 1.1  2003/06/17 19:24:10  lec
# new
#
#
