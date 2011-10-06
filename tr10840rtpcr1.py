#!/usr/local/bin/python

#
# Program: tr10840rtpcr1.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10840/SuraniTableS7.txt
#	into input files for the gelload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10840rtpcr1.py
#
# Envvars:
#
#	LOADFILE1
#	PROBEPREP_FILE
#	ASSAY_FILE
#	LANE_FILE
#	BAND_FILE
#	REFERENCE
#	CREATEDBY
#
# Inputs:
#
#       SuraniTableS7.txt (LOADFILE1), a tab-delimited file in the format:
#
#               field 0: Marker ID
#               field 1: Probe Name
#               field 2: Lane Label
#               field 3: Genotype
#	        field 4: RNAType
#	        field 5: Gel Control
#		field 6: Age
#		field 7: Lane Note
#		field 8: Structure Name
#		field 9: Theiler Stage
#		field 10 Band Strength
#		field 11: Band Note
#
# Outputs:
#
#       4 tab-delimited files:
#
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
import db

#globals

TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline
NULL = ''

inAssayFile = ''	# file descriptor

prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
laneFile = ''       # file descriptor
bandFile = ''        # file descriptor

inAssay = os.environ['LOADFILE1']

prepFileName = os.environ['PROBEPREP_FILE']
assayFileName = os.environ['ASSAY_FILE']
laneFileName = os.environ['LANE_FILE']
bandFileName = os.environ['BAND_FILE']

assayNote = '''These results were reported in Table S7.  They were obtained using TaqMan real time PCR using individual blastomeres from either 2-, 4-, or 8-cell embryos.  Two different primer sets were used with each sample in order to examine allele specific expression. In the Supplementary Materials and Methods, the authors reported that a These results were reported in Table S7.  They were obtained using TaqMan real time PCR using individual blastomeres from either 2-, 4-, or 8-cell embryos.  Two different primer sets were used with each sample in order to examine allele specific expression. In the Supplementary Materials and Methods, the authors reported that a These results were reported in Table S7.  They were obtained using TaqMan real time PCR using individual blastomeres from either 2-, 4-, or 8-cell embryos.  Two different primer sets were used with each sample in order to examine allele specific expression. In the Supplementary Materials and Methods, the authors reported that a gene was detected when the real-time PCR reaction had a Ct value of <32.  Thus, we have interpreted Ct values of <=31.5 as present and >=32.5 as absent.'''

reference = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']

# constants for primer prep
assayType = 'RT-PCR'
prepType = 'Not Specified'
hybridization = 'Not Applicable'
labelledWith = 'Not Applicable'
visualizedWith = 'Not Specified'

# gel lane
sampleAmount = ''
sex = 'Not Specified'
ageNote = ''

# gel row/band
bandSize = ''
bandUnits = 'Not Specified'
rowNote = ''

# primer lookup
primerLookup = {}

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
    global inAssayFile
    global prepFile, assayFile, laneFile, bandFile
    global primerLookup
 
    try:
        inAssayFile = open(inAssay, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inAssay)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        laneFile = open(laneFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % laneFileName)

    try:
        bandFile = open(bandFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % bandFileName)

    results = db.sql('''
	select a.accID, p.name 
	from ACC_Accession a, PRB_Probe p, PRB_Reference r
	where a._MGITYpe_key = 3
	and a._Object_key = p._Probe_key
	and p._Probe_key = r._Probe_key
	and r._Refs_key = 175863
	''', 'auto')
    for r in results:
	primerID = r['accID']
	primerName = r['name']

	primerName = string.replace(primerName, 'Kif20a,Cdc23-assay_14/19-FP_first allele', \
						'Kif20a,Cdc23-assay_14/19-FP_first allele, RP')
	primerName = string.replace(primerName, 'Kif20a,Cdc23-assay_14/19-FP_second allel', \
						'Kif20a,Cdc23-assay_14/19-FP_second allele, RP')
	primerName = string.replace(primerName, 'Klf17-assay_12/16/21-FP_second allele, R', \
						'Klf17-assay_12/16/21-FP_second allele, RP')

	primerLookup[primerName] = []
	primerLookup[primerName].append(primerID)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    assayKey = 0
    laneKey = 0
    bandKey = 1
    prevMarker = ''
    prevPrimer = ''

    # For each line in the input file

    for line in inAssayFile.readlines():

        # Split the line into tokens
        tokens = string.split(line[:-1], TAB)

	markerID = tokens[0]
	primerName = tokens[1]
	laneLabel = tokens[2]
	genotype = tokens[3]
	rnaType = tokens[4]
	gelControl = tokens[5]
	age = tokens[6]
	laneNote = tokens[7]
	structureName = tokens[8]
	structureTheilerStage = tokens[9]
	bandStrength = tokens[10]
	bandNote = tokens[11]

	if markerID != prevMarker:
	    prevMarker = markerID
	    prevPrimer = ''

	primerID = primerLookup[primerName][0]
	if primerID != prevPrimer:

	    assayKey = assayKey + 1
	    prevPrimer = primerID
            laneKey = 0

	    prepFile.write(str(assayKey) + TAB + \
	        primerID + TAB + \
	        prepType + TAB + \
	        hybridization + TAB + \
	        labelledWith + TAB + \
	        visualizedWith + CRT)

	    # write the assay information

            assayFile.write(str(assayKey) + TAB + \
                markerID + TAB + \
                reference + TAB + \
                assayType + TAB + \
                TAB + \
                assayNote + TAB + \
                createdBy + CRT)

        laneKey = laneKey + 1
	laneFile.write(str(assayKey) + TAB + \
	    str(laneKey) + TAB + \
	    laneLabel + TAB + \
	    genotype + TAB + \
	    rnaType + TAB + \
	    gelControl + TAB + \
	    sampleAmount + TAB + \
	    sex + TAB + \
	    age + TAB + \
	    ageNote + TAB + \
	    laneNote + TAB + \
	    structureName + TAB + \
	    structureTheilerStage + CRT)

	bandKey = 1
	bandFile.write(str(assayKey) + TAB + \
	    str(laneKey) + TAB + \
	    str(bandKey) + TAB + \
	    bandSize + TAB + \
	    bandUnits + TAB + \
	    bandStrength + TAB + \
	    rowNote + TAB +
	    bandNote + CRT)

    inAssayFile.close()
    assayFile.close()
    laneFile.close()
    bandFile.close()

    # end of "for line in inAssayFile.readlines():"

# Main
#

init()
process()
exit(0)

