#!/usr/local/bin/python

# $Header$
# $Name$

#
# Program: J80501rtpcr.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RT-PCR input files into input files
#	for the assayload/gelload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	J80501rtpcr.py
#
# Envvars:
#
# Inputs:
#
#       RT-PCR.txt, a tab-delimited file in the format:
#		field 1: Human Gene Symbol
#		field 2: Mouse Gene Symbol
#               field 3: MGI Marker Accession ID
#               field 4: Sequence 1
#               field 5: Sequence 2
#               field 6-7: Tissue expression
#
#	RT_PCR_tissue.txt, a tab-delimited file in the format:
#		field 1: Reported Tissue
#		field 2: MGI Tissue
#		field 3: Theiler Stage
#
# Outputs:
#
#       5 tab-delimited files:
#
#	RT_PCR_primer.txt
#	RT_PCR_probeprep.txt
#	RT_PCR_assay.txt
#	RT_PCR_gellane.txt
#	RT_PCR_gelband.txt
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

inRTPCRFile = ''	# file descriptor
inTissueFile = ''	# file descriptor
primerFile = ''         # file descriptor
prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
gellaneFile = ''        # file descriptor
gelbandFile = ''        # file descriptor

datadir = os.environ['RTPCRDATADIR']

inRTPCRFileName = datadir + '/tr5154/RT-PCR.txt'
inTissueFileName = datadir + '/tr5154/RT_PCR_tissue.txt'
primerFileName = datadir + '/RT_PCR_primer.txt'
prepFileName = datadir + '/RT_PCR_probeprep.txt'
assayFileName = datadir + '/RT_PCR_assay.txt'
gellaneFileName = datadir + '/RT_PCR_gellane.txt'
gelbandFileName = datadir + '/RT_PCR_gelband.txt'

# constants for primers
regionCovered = NULL
repeatUnit = NULL
moreProduct = '0'
productSize = NULL

# constants for probe prep
prepType = 'Not Specified'
hybridization = 'Not Applicable'
labelledWith = 'Not Applicable'
labelCoverage = 'Not Applicable'
visualizedWith = 'Not Specified'

# constants for assay
reference = 'J:80501'
assayType = 'RT-PCR'
createdBy = os.environ['CREATEDBY']

# constants for gel lanes
control = 'No'
genotype = 'MGI:2166310'	# Not Specified
sampleAmount = '2'
rnaType = 'total'
embryonicAgeNote = 'Age of embryo at noon of plug day not specified in reference.'
sex = 'Not Specified'
laneNote = NULL
rowNote = NULL
bandNote = NULL

# constants for gel rows
gelsize = NULL
gelunits = 'Not Specified'

mgiTypeKey = 8		# Assay
mgiPrefix = "MGI:"

# translation of ages to theiler stages

ageTrans = {'28':'postnatal day 2', \
    '23':'embryonic day 15.5'}

# translation of input file band strengths and MGI band strengths

bandTrans = {'+':'Present', '-':'Absent'}

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
    global inRTPCRFile, inTissueFile, primerFile, prepFile, assayFile, gellaneFile, gelbandFile
 
    try:
        inRTPCRFile = open(inRTPCRFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inRTPCRFileName)

    try:
        inTissueFile = open(inTissueFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inTissueFileName)

    try:
        primerFile = open(primerFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % primerFileName)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        gellaneFile = open(gellaneFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % gellaneFileName)

    try:
        gelbandFile = open(gelbandFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % gelbandFileName)

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

#    assay = 0	# unique Assay ID, if heading exists
    assay = 1	# unique Assay ID

    # For each line in the input file

    for line in inRTPCRFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	# processing first line (header)
	# grab the Tissue headings

	# hard-code headings for this file
	tissueLabels = ['E15.5', 'P2']

	if assay == 0:
	    tissueLabels = tokens[5:-1]
	    assay = assay + 1
	    continue

	# else process an actual data line

        try:
	    humanGene = tokens[0]
	    mouseGene = tokens[1]
	    accID = tokens[2]
	    seq1 = tokens[3]
	    seq2 = tokens[4]
	    bands = tokens[5:]

        except:
            print 'Invalid Line (%d): %s\n' % (assay, line)

	if len(mouseGene) == 0:
	    continue

	# write the primer information

	primerFile.write(str(assay) + TAB + \
	    accID + TAB + \
	    mouseGene + '-pA,' + mouseGene + '-pB' + TAB + \
	    reference + TAB + \
	    regionCovered + TAB + \
	    seq1 + TAB + \
	    seq2 + TAB + \
	    repeatUnit + TAB + \
	    moreProduct + TAB + \
	    productSize + CRT)

	# write the probe prep information

	prepFile.write(str(assay) + TAB + \
	    prepType + TAB + \
	    hybridization + TAB + \
	    labelledWith + TAB + \
	    labelCoverage + TAB + \
	    visualizedWith + CRT)

	# write the assay information

	assayFile.write(str(assay) + TAB + \
	    accID + TAB + \
	    reference + TAB + \
	    assayType + TAB + \
	    createdBy + CRT)

	# write the gel lane and gel band information

	lane = 1
	row = 1

	# for each Tissue from the heading (Brain, Heart, etc.)

	for i in range(len(tissueLabels)):

	    # Translate the Tissue into a Tissue and Age
	    [tissue, theilerStage] = string.split(tissueTrans[tissueLabels[i]], '|')

	    if theilerStage == '28':
		ageNote = NULL
            else:
		ageNote = embryonicAgeNote

	    gellaneFile.write(str(assay) + TAB + \
		str(lane) + TAB + \
		tissueLabels[i] + TAB + \
		genotype + TAB + \
		rnaType + TAB + \
		control + TAB + \
		sampleAmount + TAB + \
		sex + TAB + \
		ageTrans[theilerStage] + TAB + \
		ageNote + TAB + \
		laneNote + TAB + \
		tissue + TAB + \
		theilerStage + CRT)

	    gelbandFile.write(str(assay) + TAB + \
	        str(lane) + TAB + \
	        str(row) + TAB + \
	        gelsize + TAB + \
	        gelunits + TAB + \
	        bandTrans[bands[i]] + TAB + \
	        rowNote + TAB + \
	        bandNote + CRT)

	    lane = lane + 1

	assay = assay + 1

    #	end of "for line in inRTPCRFile.readlines():"

#
# Main
#

init()
process()
exit(0)

# $Log$
# Revision 1.1  2003/09/19 18:09:02  lec
# TR 5154
#
#
