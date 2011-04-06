#!/usr/local/bin/python

#
#  tr10537.py
###########################################################################
#
#  Purpose:
#
#      This script will create the input files for the GXD in situ load.
#
#  Usage:
#
#      tr10537.py
#
#  Env Vars:
#
#      PROBEPREP_DATA
#      ASSAY_DATA
#      SPECIMEN_DATA
#      RESULTS_DATA
#      STRENGTH_LOOKUP
#      PATTERN_LOOKUP
#      STRUCTURE_LOOKUP
#      IMAGE_LOOKUP
#      PROBEPREP_FILE
#      ASSAY_FILE
#      SPECIMEN_FILE
#      RESULTS_FILE
#      REFERENCE
#      CREATEDBY
#
#  Inputs:
#
#      ProbePrep.txt - Tab-delimited fields:
#
#          1) Assay ID ("g#####")
#          2) Probe MGI ID
#          3) Probe Prep Type (all are "RNA")
#          4) Hybridization (all are "Antisense")
#          5) Labelled With (all are "P33")
#          6) Visualized With (all are "Autoradiography")
#
#      Assays.txt - Tab-delimited fields:
#
#          1) Assay ID ("g#####")
#          2) MGI ID for the assayed gene
#          3) Reference (all are "J:162220")
#          4) Assay Type (all are "RNA in situ")
#          5) Reporter Gene (all are null)
#          6) Assay Note (all are null)
#          7) Created By (all are "jfinger")
#
#      Specimen.txt - Tab-delimited fields:
#
#          1) Assay ID ("g#####")
#          2) Specimen ID ("g##### " + "E11.5", "E15.5", "P7" or "Adult")
#          3) Specimen Label (same as specimen number)
#          4) Genotype (all are "MGI:2166522")
#          5) Age ("embryonic day 11.5", "embryonic day 15.5",
#                  "postnatal day 7" or "postnatal day 42")
#          6) Age Note (all are null)
#          7) Sex (all are "Not Specified")
#          8) Fixation (all are "4% Paraformaldehyde")
#          9) Embedding Method (all are "Cryosection")
#          10) Hybridization (all are "section")
#          11) Specimen Note (all are null)
#
#      SpecimenResults.txt - Tab-delimited fields:
#
#          1) Assay ID ("g#####")
#          2) Specimen ID ("g##### " + "E11.5", "E15.5", "P7" or "Adult")
#          3) Result ID (all are null)
#          4) BGEM Strength (code from StrengthLookupTable.txt)
#          5) BGEM Pattern (code from PatternLookupTable.txt)
#          6) BGEM Structure (code from StructureLookupTable.txt)
#          7) MGI Structure Theiler Stage ("19", "23" or "28")
#          8) Result Note (all are null)
#
#      StrengthLookupTable.txt - Tab-delimited fields:
#
#          1) BGEM Strength Code
#          2) GXD Strength Name
#
#      PatternLookupTable.txt - Tab-delimited fields:
#
#          1) BGEM Pattern Code
#          2) GXD Pattern Name
#
#      StructureLookupTable.txt - Tab-delimited fields:
#
#          1) BGEM Structure
#          2) GXD Print Name
#
#      ImageLookupTable.txt - Tab-delimited fields:
#
#          1) Specimen ID ("g##### " + "Adult","P7","E15.5" or "E11.5")
#          2) Figure Labels for images (comma-separated)
#
#  Outputs:
#
#      In_Situ_probeprep.txt - Tab-delimited fields:
#
#          1) Assay Number
#          2) Probe MGI ID
#          3) Probe Prep Type
#          4) Hybridization
#          5) Labelled With
#          6) Visualized With
#
#      In_Situ_assay.txt - Tab-delimited fields:
#
#          1) Assay Number
#          2) Marker MGI ID
#          3) Reference (J:#####)
#          4) Assay Type
#          5) Reporter Gene
#          6) Assay Note
#          7) Created By
#
#      In_Situ_specimen.txt - Tab-delimited fields:
#
#          1) Assay Number
#          2) Specimen Number
#          3) Specimen Label
#          4) Genotype ID
#          5) Age
#          6) Age Note
#          7) Sex
#          8) Fixation
#          9) Embedding Method
#          10) Hybridization
#          11) Specimen Note
#
#      In_Situ_results.txt - Tab-delimited fields:
#
#          1) Assay Number
#          2) Specimen Number
#          3) Result Number
#          4) Strength
#          5) Pattern
#          6) MGI Structure Name
#          7) MGI Structure Theiler Stage
#          8) Result Note
#          9) Image Figure Labels (comma-separated)
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes:  None
#
###########################################################################

import sys
import os
import string
import re
import db

#
#  GLOBALS
#
probePrepData = os.environ['PROBEPREP_DATA']
assayData = os.environ['ASSAY_DATA']
specimenData = os.environ['SPECIMEN_DATA']
resultsData = os.environ['RESULTS_DATA']

strengthLookupFile = os.environ['STRENGTH_LOOKUP']
patternLookupFile = os.environ['PATTERN_LOOKUP']
structureLookupFile = os.environ['STRUCTURE_LOOKUP']
imageLookupFile = os.environ['IMAGE_LOOKUP']

probePrepFile = os.environ['PROBEPREP_FILE']
assayFile = os.environ['ASSAY_FILE']
specimenFile = os.environ['SPECIMEN_FILE']
resultsFile = os.environ['RESULTS_FILE']

jNumber = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']


#
# Purpose: Build a dictionary lookup for each of the lookup files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables
# Throws: Nothing
#
def buildLookups ():
    global strengthLookup, patternLookup, structureLookup, imageLookup

    strengthLookup = {}
    patternLookup = {}
    structureLookup = {}
    imageLookup = {}

    #
    # Open the lookup files.
    #
    try:
        fpStrengthLookup = open(strengthLookupFile, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + strengthLookupFile + '\n')
        sys.exit(1)

    try:
        fpPatternLookup = open(patternLookupFile, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + patternLookupFile + '\n')
        sys.exit(1)

    try:
        fpStructureLookup = open(structureLookupFile, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + structureLookupFile + '\n')
        sys.exit(1)

    try:
        fpImageLookup = open(imageLookupFile, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + imageLookupFile + '\n')
        sys.exit(1)

    #
    # Build the lookups.
    #
    line = fpStrengthLookup.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        strengthLookup[tokens[0]] = tokens[1]
        line = fpStrengthLookup.readline()

    line = fpPatternLookup.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        patternLookup[tokens[0]] = tokens[1]
        line = fpPatternLookup.readline()

    line = fpStructureLookup.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        structureLookup[tokens[0]] = tokens[1]
        line = fpStructureLookup.readline()

    line = fpImageLookup.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        imageLookup[tokens[0]] = tokens[1]
        line = fpImageLookup.readline()

    #
    # Close the lookup files.
    #
    fpStrengthLookup.close()
    fpPatternLookup.close()
    fpStructureLookup.close()
    fpImageLookup.close()

    return


#
# Purpose: Open the files.
# Returns: Nothing
# Assumes: The names of the files are set in the environment.
# Effects: Sets global variables
# Throws: Nothing
#
def openFiles ():
    global fpProbePrepData, fpAssayData, fpSpecimenData, fpResultsData
    global fpProbePrepFile, fpAssayFile, fpSpecimenFile, fpResultsFile

    #
    # Open the input files.
    #
    try:
        fpProbePrepData = open(probePrepData, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + probePrepData + '\n')
        sys.exit(1)

    try:
        fpAssayData = open(assayData, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + assayData + '\n')
        sys.exit(1)

    try:
        fpSpecimenData = open(specimenData, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + specimenData + '\n')
        sys.exit(1)

    try:
        fpResultsData = open(resultsData, 'r')
    except:
        sys.stderr.write('Cannot open file: ' + resultsData + '\n')
        sys.exit(1)

    #
    # Open the output files.
    #
    try:
        fpProbePrepFile = open(probePrepFile, 'w')
    except:
        sys.stderr.write('Cannot open file: ' + probePrepFile + '\n')
        sys.exit(1)

    try:
        fpAssayFile = open(assayFile, 'w')
    except:
        sys.stderr.write('Cannot open file: ' + assayFile + '\n')
        sys.exit(1)

    try:
        fpSpecimenFile = open(specimenFile, 'w')
    except:
        sys.stderr.write('Cannot open file: ' + specimenFile + '\n')
        sys.exit(1)

    try:
        fpResultsFile = open(resultsFile, 'w')
    except:
        sys.stderr.write('Cannot open file: ' + resultsFile + '\n')
        sys.exit(1)

    return


#
# Purpose: Close the files.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def closeFiles ():
    fpProbePrepData.close()
    fpAssayData.close()
    fpSpecimenData.close()
    fpResultsData.close()
    fpProbePrepFile.close()
    fpAssayFile.close()
    fpSpecimenFile.close()
    fpResultsFile.close()

    return


#
# Purpose: Use the records in the input files to create the probe prep,
#          assay, specimen and results output files.
# Returns: 0 if successful, 1 for an error
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def process ():

    #
    # Initialize assay number and lookup.
    #
    # There is a 1-to-1 relationship between records in the probe prep
    # input file and records in the assay input file. The records are
    # linked by the assay ID in column 1. We'll generate the sequential
    # assay numbers while processing the probe prep records and create
    # an assay number lookup that can be used to find the proper assay
    # number when processing the other files.
    #
    assayNumber = 0
    assayLookup = {}

    #
    # Process the input records from the probe prep input file.
    #
    line = fpProbePrepData.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        assayID = tokens[0]
        probeID = tokens[1]
        probePrepType = tokens[2]
        hybridization = tokens[3]
        labelledWith = tokens[4]
        visualizedWith = tokens[5]

        #
        # Use the next unique assay number for each output record.
        #
        assayNumber += 1

        #
        # Write a record to the probe prep file.
        #
        fpProbePrepFile.write(str(assayNumber) + '\t' +
                          probeID + '\t' +
                          probePrepType + '\t' +
                          hybridization + '\t' +
                          labelledWith + '\t' +
                          visualizedWith + '\n')

        #
        # Save the assay number that was used for this assay ID in the
        # lookup dictionary.
        #
        assayLookup[assayID] = assayNumber

        line = fpProbePrepData.readline()

    #
    # Process the input records from the assay input file.
    #
    line = fpAssayData.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        assayID = tokens[0]
        markerID = tokens[1]
        jNumber = tokens[2]
        assayType = tokens[3]
        reporterGene = tokens[4]
        assayNote = tokens[5]
        createdBy = tokens[6]

        #
        # Look up the assay number to be used for this assay ID.
        #
        assayNumber = assayLookup[assayID]

        #
        # Write a record to the assay file.
        #
        fpAssayFile.write(str(assayNumber) + '\t' +
                      markerID + '\t' +
                      jNumber + '\t' +
                      assayType + '\t' +
                      reporterGene + '\t' +
                      assayNote + '\t' +
                      createdBy + '\n')

        line = fpAssayData.readline()

    #
    # Initialize specimen number and lookup.
    #
    assayID_old = ''
    specimenNumber = 0
    specimenLookup = {}

    #
    # Process the input records from the specimen input file.
    #
    line = fpSpecimenData.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        assayID = tokens[0]
        specimenID = tokens[1]
        specimenLabel = tokens[2]
        genotypeID = tokens[3]
        age = tokens[4]
        ageNote = tokens[5]
        sex = tokens[6]
        fixation = tokens[7]
        embeddingMethod = tokens[8]
        hybridization = tokens[9]
        specimenNote = tokens[10]

        #
        # Look up the assay number to be used for this assay ID.
        #
        assayNumber = assayLookup[assayID]

        #
        # If this record does not have a new assay, use the next available
        # sequential number for the specimen. Otherwise, reset the specimen
        # number for a new assay.
        #
        if assayID == assayID_old:
            specimenNumber += 1
        else:
            specimenNumber = 1
            assayID_old = assayID

        #
        # Write a record to the specimen file.
        #
        fpSpecimenFile.write(str(assayNumber) + '\t' +
                         str(specimenNumber) + '\t' +
                         specimenLabel + '\t' +
                         genotypeID + '\t' +
                         age + '\t' +
                         ageNote + '\t' +
                         sex + '\t' +
                         fixation + '\t' +
                         embeddingMethod + '\t' +
                         hybridization  + '\t' +
                         specimenNote + '\n')

        #
        # Save the specimen number that was used for this specimen ID in the
        # lookup dictionary.
        #
        specimenLookup[specimenID] = specimenNumber

        line = fpSpecimenData.readline()

    #
    # Initialize result number.
    #
    specimenID_old = ''
    resultNumber = 0

    #
    # Process the input records from the results input file.
    #
    line = fpResultsData.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        assayID = tokens[0]
        specimenID = tokens[1]
        resultID = tokens[2]
        strengthCode = tokens[3]
        patternCode = tokens[4]
        structureCode = tokens[5]
        theilerStage = tokens[6]
        resultNote = tokens[7]

        #
        # Look up the assay number to be used for this assay ID.
        #
        assayNumber = assayLookup[assayID]

        #
        # Look up the specimen number to be used for this specimen ID.
        #
        specimenNumber = specimenLookup[specimenID]

        #
        # If this record does not have a new specimen, use the next available
        # sequential number for the result. Otherwise, reset the result
        # number for a new specimen.
        #
        if specimenID == specimenID_old:
            resultNumber += 1
        else:
            resultNumber = 1
            specimenID_old = specimenID

        #
        # Look up the strength, pattern, structure and image list.
        #
        strength = strengthLookup[strengthCode]
        pattern = patternLookup[patternCode]
        structure = structureLookup[structureCode]
        if imageLookup.has_key(specimenID):
            imageList = imageLookup[specimenID]
        else:
            imageList = ''

        #
        # Write a record to the results file.
        #
        fpResultsFile.write(str(assayNumber) + '\t' +
                        str(specimenNumber) + '\t' +
                        str(resultNumber) + '\t' +
                        strength + '\t' +
                        pattern + '\t' +
                        structure + '\t' +
                        theilerStage + '\t' +
                        resultNote + '\t' +
                        imageList + '\n')

        line = fpResultsData.readline()

    return 0


#
# Main
#
buildLookups()
openFiles()
process()
closeFiles()

sys.exit(0)
