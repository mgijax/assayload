#!/usr/local/bin/python

#
# Program: tr7819insitu.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RNA InSitu input file into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr7819insitu.py
#
# Envvars:
#
# Inputs:
#
#	DgenTable.txt , a tab-delimited file in the format:
#		field 1: Image (w/out .jpg)
#		field 2: Hybridization
#		field 3: Embedding
#		field 4: Marker
#		field 5: Genotype
#		field 6: Age_range
#		field 7: Sex
#		field 8: Organ System
#		field 9: Organ
#		field 10: Organ Structure
#		field 11: Cell Type
#		field 12: Intensity
#		field 13: Distribution
#		field 14: Mouse
#
#       DgenResults.txt, a tab-delimited file in the format:
#               field 1: Organ System
#		field 2: Organ
#		field 3: Organ Structure
#		field 4: Cell Type
#		field 5: Intensity
#		field 6: Distribution
#		field 7: TS28_structure1
#		field 8: TS28_structure2
#		field 9: Pattern
#		field 10: Strength
#		field 11: Results Notes
#
# Outputs:
#
#       3 tab-delimited files:
#
#	In_Situ_assay.txt
#	In_Situ_specimen.txt
#	In_Situ_results.txt
#	In_Situ_probeprep.txt
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

inResultsFile = ''	# file descriptor
inTableFile = ''	# file descriptor
inPixFile = ''		# file descriptor

assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['DATADIR']

inResultsFileName = datadir + '/DgenResults.txt'
inTableFileName = datadir + '/DgenTable.txt'
inPixFileName = datadir + '/pix7819.txt'

assayFileName = datadir + '/In_Situ_assay.txt'
specimenFileName = datadir + '/In_Situ_specimen.txt'
resultsFileName = datadir + '/In_Situ_results.txt'
prepFileName = datadir + '/In_Situ_probeprep.txt'

# constants for Absent
ABSENT = 'Absent'
NA = 'Not Applicable'

# constants for assay
reference = 'J:101679'
assayType = 'In situ reporter (knock in)'
createdBy = os.environ['CREATEDBY']

# constants for knock-in
reporterGene = 'lacZ'
# since this is 'direct detection' there is no antibody prep or probe prep

# constants for specimen
age = 'postnatal day %s'
ageNote = NULL
fixation = 'Not Specified'
specimenNote = NULL
theilerStage = 28

# dictionary to lookup images
pixelDict = {}

# diciontary to lookup results
resultsDict = {}

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
 
    inTableFile.close()
    assayFile.close()
    specimenFile.close()
    resultsFile.close()
    prepFile.close()

    sys.exit(status)
 
# Purpose: initialize
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global inResultsFile, inTableFile, inPixFile
    global assayFile, specimenFile, resultsFile, prepFile
    global pixelDict, resultsDict
 
    try:
        inTableFile = open(inTableFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inTableFileName)

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

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    #
    # lookup file for images
    #

    try:
        inPixFile = open(inPixFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inPixFileName)

    for line in inPixFile.readlines():
        tokens = string.split(line[:-1], TAB)
        pixFileName = tokens[0]
        pixID = tokens[1]
        pixelDict[pixFileName] = pixID
    inPixFile.close()

    #
    # cache results
    #

    try:
        inResultsFile = open(inResultsFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inResultsFileName)

    for line in inResultsFile.readlines():

	tokens = string.split(line[:-1], TAB)

        organSystem = tokens[0]
	organ = tokens[1]
	organStructure = tokens[2]
	cellType = tokens[3]
	intensity = tokens[4]
	distribution = tokens[5]
	structure1 = tokens[6]
	structure2 = tokens[7]
	pattern = tokens[8]
	strength = tokens[9]
	resultsNotes = tokens[10]

	key = organSystem + organ + organStructure + cellType + intensity + distribution
	value = {}
	value['structure1'] = structure1
	value['structure2'] = structure2
	value['pattern'] = pattern
	value['strength'] = strength
	value['resultsNotes'] = resultsNotes

	resultsDict[key] = value

    inResultsFile.close()

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0
    specimenKey = 0
    assayList = []
    specimenList = []
    specimenKeyList = {}

    # For each line in the input file

    lineNum = 0
    for line in inTableFile.readlines():

	lineNum = lineNum + 1

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)

	if lineNum == 1:
	    continue

	# else process an actual data line

        try:
	    image = tokens[0]
	    imageJPG = image + '.jpg'
	    hybridization = tokens[1]
	    embedding = tokens[2]
	    marker = tokens[3]
	    genotype = tokens[4]
	    ageRange = tokens[5]
	    sex = tokens[6]
	    organSystem = tokens[7]
	    organ = tokens[8]
	    organStructure = tokens[9]
	    cellType = tokens[10]
	    intensity = tokens[11]
	    distribution = tokens[12]
	    mouse = tokens[13]

        except:
            print 'Invalid Line (%d): %s\n' % (lineNum, line)

	# write out distinct assays, based on marker

	if marker not in assayList:

	    assayKey = assayKey + 1
	    specimenKey = 0

	    assayFile.write(str(assayKey) + TAB + \
	        marker + TAB + \
	        reference + TAB + \
	        assayType + TAB + \
	        reporterGene + TAB + \
	        TAB + \
	        createdBy + CRT)

	    assayList.append(marker)

        # write distinct specimens, based on specimenLabel

	if pixelDict.has_key(imageJPG):
	    specimenLabel = 'Mouse ' + mouse + '-' + image
	else:
	    specimenLabel = 'Mouse ' + mouse

	if specimenLabel not in specimenList:

	    specimenKey = specimenKey + 1

	    specimenFile.write(str(assayKey) + TAB + \
	        str(specimenKey) + TAB + \
	        specimenLabel + TAB + \
	        genotype + TAB + \
	        age % (ageRange) + TAB + \
	        ageNote + TAB + \
	        sex + TAB + \
	        fixation + TAB + \
	        embedding + TAB + \
	        hybridization + TAB + \
	        specimenNote + CRT)

	    specimenList.append(specimenLabel)

	resultKey = organSystem + organ + organStructure + cellType + intensity + distribution
	key = '%s:%s:%s' % (assayKey, specimenKey, imageJPG)

	if not specimenKeyList.has_key(key):
	    specimenKeyList[key] = []
	specimenKeyList[key].append(resultKey)

    # end of "for line in inTable.readlines():"

    # results

    specimenKeyList.keys().sort

    for k in specimenKeyList.keys():

        resultKey = 0

	for rKey in specimenKeyList[k]:

	    (aKey, sKey, imageJPG) = string.split(k, ':')

	    r = resultsDict[rKey]

            resultKey = resultKey + 1
    
            resultsFile.write(str(aKey) + TAB + \
	            str(sKey) + TAB + \
	            str(resultKey) + TAB + \
		    r['strength'] + TAB + \
		    r['pattern'] + TAB + \
		    r['structure1'] + TAB + \
		    str(theilerStage) + TAB + \
		    r['resultsNotes'] + TAB + \
		    imageJPG + CRT)

	    if len(r['structure2']) > 0:

	        resultKey = resultKey + 1

	        resultsFile.write(str(aKey) + TAB + \
	                str(sKey) + TAB + \
	                str(resultKey) + TAB + \
		        r['strength'] + TAB + \
		        r['pattern'] + TAB + \
		        r['structure2'] + TAB + \
		        str(theilerStage) + TAB + \
		        r['resultsNotes'] + TAB + \
		        imageJPG + CRT)

#
# Main
#

init()
process()
exit(0)

