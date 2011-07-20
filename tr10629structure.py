#!/usr/local/bin/python

#
# Program: tr10629structure.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate 
#		10629/LoaderDocuments/LoadFile3.txt 
#		10629/LoaderDocuments/StructureLookup.txt 
#	into input files for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	tr10629structure.py
#
# Envvars:
#
#	LOADFILE3
#	STRUCTURE_LOOKUP
#
# Inputs:
#
#
#       LoadFile3.txt, a tab-delimited file in the format:
#
#               field 0: Specimen Label
#		field 1: EMAP_ID Structure_Name
#		field 2: EMAP name
#		field 3: Strength
#		field 4: Pattern
#		field 5: Result Note
#		field 6: Image Names
#
#	StructureLookup.txt, a tab-delimited file in the format:
#
#		field 0: EMAP id
#		field 1: GUDMAP name
#		field 2: MGI ID of the structure
#		field 3: Theiler Stage
#		field 4: GXD Structure Print Name
#
# Outputs:
#
#       LoadFile4.txt with (LoadFile3.txt plus 2 new columns):
#
#               field 0: Specimen Label
#		field 1: EMAP_ID Structure_Name
#		field 2: EMAP name
#		field 3: Strength
#		field 4: Pattern
#		field 5: Result Note
#		field 6: Image Names
#		field 7: MGD Structure Print Name
#		field 8: MGD Structure Theiler Stage
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

inLoadFile3 = ''
inStructureFile = ''
outLoadFile4 = ''
inLoadFile3Name = os.environ['LOADFILE3']
inStructureFileName = os.environ['STRUCTURE_LOOKUP']
outLoadFile4Name = os.environ['LOADFILE4']

# structure lookup
structureLookup = {}

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
    global inLoadFile3, inStructureFile, outLoadFile4
    global structureLookup
 
    try:
        inLoadFile3 = open(inLoadFile3Name, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inLoadFile3Name)

    try:
        inStructureFile = open(inStructureFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inStructureFileName)

    try:
        outLoadFile4 = open(outLoadFile4Name, 'w')
    except:
        exit(1, 'Could not open file %s\n' % outLoadFile4Name)

    # structure lookup
    for line in inStructureFile.readlines():
	tokens = string.split(line[:-1], TAB)
	structureID = int(string.replace(tokens[0], 'EMAP:', ''))
	structureLookup[structureID] = []
	structureLookup[structureID].append(tokens)
    inStructureFile.close()

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    resultLookup = {}

    for rline in inLoadFile3.readlines():

        # Split the line into tokens
        tokens = string.split(rline[:-1], TAB)

	specimenID = tokens[0]
	emapID = int(string.replace(tokens[1], 'EMAP:', ''))
	emapName = tokens[2]
	strength = tokens[3]
	pattern = tokens[4]
	resultNote = tokens[5]
	imageName = tokens[6]

	r = structureLookup[emapID]
	structureName = r[0][4]
	structureTheilerStage = r[0][3]

	outLoadFile4.write(specimenID + TAB + \
	    'EMAP:' + str(emapID) + TAB + \
	    emapName + TAB + \
	    strength + TAB + \
	    pattern + TAB + \
	    resultNote + TAB + \
	    imageName + TAB + \
	    structureName + TAB + \
	    structureTheilerStage + CRT)

    inLoadFile3.close()
    outLoadFile4.close()

    # end of "for line in inLoadFile3.readlines():"
#
# Main
#

init()
process()
exit(0)

