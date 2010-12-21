#!/usr/local/bin/python

#
#  tr10407.py
###########################################################################
#
#  Purpose:
#
#      This script will create the input files for the GXD in situ load.
#
#  Usage:
#
#      tr10407.py
#
#  Env Vars:
#
#      ASSAY_DATA
#      SPEC_RESULTS_DATA
#      PROBEPREP_FILE
#      ASSAY_FILE
#      SPECIMEN_FILE
#      RESULTS_FILE
#      REFERENCE
#      CREATEDBY
#
#  Inputs:
#
#      Assays.txt - Tab-delimited fields:
#
#          1) MGI ID for the assayed gene
#          2) Probe Prep Type (all are "RNA")
#          3) Hybridization (all are "Antisense")
#          4) Labelled With (all are "digoigenin")
#          5) Visualized With (all are "Not Specified")
#          6) Reference (all are "J:156017")
#          7) Assay Type (all are "RNA in situ")
#          8) Reporter Gene (all are null)
#          9) Assay Note (all are null)
#          10) Created By (all are "cms")
#
#      SpecimenResults.txt - Tab-delimited fields:
#
#          1) MGI ID for the assayed gene
#          2) Specimen Label
#          3) Genotype (all are "ICR")
#          4) Age (either embryonic day 9.5, 10.5 or 11.5)
#          5) Age Note (all are null)
#          6) Sex (all are "Not Specified")
#          7) Fixation (all are "4% Paraformaldehyde")
#          8) Embedding Method (all are "Not Applicable")
#          9) Hybridization (all are "whole mount")
#          10) Specimen Note (all are null)
#          11) Strength (for result 1)
#          12) Pattern (for result 1)
#          13) MGI Structure (for result 1)
#          14) Theiler Stage (for result 1)
#          15) Note (for result 1)
#          .
#          .  There are 20 sets of result data (repeat columns 11-15).
#          .
#          106) Strength (for result 20)
#          107) Pattern (for result 20)
#          108) MGI Structure (for result 20)
#          109) Theiler Stage (for result 20)
#          110) Note (for result 20)
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
#  CONSTANTS
#

# Column indexes in the Specimen/Result file
FIRST_RESULT_IDX = 10     # Index where the first result starts
LAST_RESULT_IDX = 105     # Index where the last result starts

# Date the probes were created
PROBE_CREATION_DATE = "12/14/2010"

#
#  GLOBALS
#
inputFile1 = os.environ['PROBEPREP_ASSAY_DATA']
inputFile2 = os.environ['SPEC_RESULTS_DATA']

probePrepFile = os.environ['PROBEPREP_FILE']
assayFile = os.environ['ASSAY_FILE']
specimenFile = os.environ['SPECIMEN_FILE']
resultsFile = os.environ['RESULTS_FILE']

jNumber = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']


#
# Purpose: Perform initialization tasks.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables
# Throws: Nothing
#
def initialize ():
    global refKey

    #
    # Get the reference key for the J-Number.
    #
    results = db.sql('select _Object_key ' + \
                     'from ACC_Accession ' + \
                     'where accID = "' + jNumber + '" and ' + \
                           '_LogicalDB_key = 1 and ' + \
                           '_MGIType_key = 1', 'auto')
    refKey = results[0]['_Object_key']


#
# Purpose: Build a lookup to find the MGI ID of a probe that was created
#          for a given marker MGI ID.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables
# Throws: Nothing
#
def buildProbeLookup ():
    global probeLookup

    probeLookup = {}

    #
    # Get the MGI IDs of the new probes and the corresponding markers.
    #
    cmd = 'select a1.accID "probeID", ' + \
                 'a2.accID "markerID" ' + \
          'from PRB_Probe p, ' + \
               'PRB_Reference pr, ' + \
               'PRB_Marker pm, ' + \
               'ACC_Accession a1, ' + \
               'ACC_Accession a2 ' + \
          'where p.creation_date = "%s" and ' + \
                'p._Probe_key = pr._Probe_key and ' + \
                'pr._Refs_key = %d and ' + \
                'p._Probe_key = pm._Probe_key and ' + \
                'pm._Probe_key = a1._Object_key and ' + \
                'a1._MGIType_key = 3 and ' + \
                'a1._LogicalDB_key = 1 and ' + \
                'a1.prefixPart = "MGI:" and ' + \
                'a1.preferred = 1 and ' + \
                'pm._Marker_key = a2._Object_key and ' + \
                'a2._MGIType_key = 2 and ' + \
                'a2._LogicalDB_key = 1 and ' + \
                'a2.prefixPart = "MGI:" and ' + \
                'a2.preferred = 1'
    results = db.sql(cmd % (PROBE_CREATION_DATE,refKey),'auto')

    #
    # Add each marker/probe to the dictionary.
    #
    for r in results:
        probeLookup[r['markerID']] = r['probeID']

    return


#
# Purpose: Open the files.
# Returns: Nothing
# Assumes: The names of the files are set in the environment.
# Effects: Sets global variables
# Throws: Nothing
#
def openFiles ():
    global fp1, fp2
    global fpProbePrepFile, fpAssayFile, fpSpecimenFile, fpResultsFile

    #
    # Open the input files.
    #
    try:
        fp1 = open(inputFile1, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + inputFile1 + '\n')
        sys.exit(1)

    try:
        fp2 = open(inputFile2, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + inputFile2 + '\n')
        sys.exit(1)

    #
    # Open the output files.
    #
    try:
        fpProbePrepFile = open(probePrepFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + probePrepFile + '\n')
        sys.exit(1)

    try:
        fpAssayFile = open(assayFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + assayFile + '\n')
        sys.exit(1)

    try:
        fpSpecimenFile = open(specimenFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + specimenFile + '\n')
        sys.exit(1)

    try:
        fpResultsFile = open(resultsFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + resultsFile + '\n')
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
    fp1.close()
    fp2.close()
    fpProbePrepFile.close()
    fpAssayFile.close()
    fpSpecimenFile.close()
    fpResultsFile.close()

    return


#
# Purpose: Use the records in input file 1 to create the probe prep and
#          assay output files. Use the records in input file 2 to create
#          the specimen and results output files.
# Returns: 0 if successful, 1 for an error
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def process ():

    #
    # Process the records in input file 1 to create the probe prep and assay
    # output files.
    #
    assayNumber = 0
    assayLookup = {}

    line = fp1.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        markerID = tokens[0]
        probPrepType = tokens[1]
        hybridization = tokens[2]
        labelledWith = tokens[3]
        visualizedWith = tokens[4]
        # skip reference field (tokens[5]) and use the config setting.
        assayType = tokens[6]
        reportedGene = tokens[7]
        assayNote = tokens[8]
        # skip reference field (tokens[9]) and use the config setting.

        probeID = probeLookup[markerID]

        #
        # Use the next unique assay number for this assay.
        #
        assayNumber += 1

        #
        # Write a record to the probe prep file.
        #
        fpProbePrepFile.write(str(assayNumber) + '\t' +
                          probeID + '\t' +
                          probPrepType + '\t' +
                          hybridization + '\t' +
                          labelledWith + '\t' +
                          visualizedWith + '\n')

        #
        # Write a record to the assay file.
        #
        fpAssayFile.write(str(assayNumber) + '\t' +
                      markerID + '\t' +
                      jNumber + '\t' +
                      assayType + '\t' +
                      reportedGene + '\t' +
                      assayNote + '\t' +
                      createdBy + '\n')
        #
        # Maintain a lookup for finding all assay numbers that are used
        # for each marker ID. This will be needed when processing
        # input file 2 to create the specimen and results files.
        #
        if assayLookup.has_key(markerID):
            list = assayLookup[markerID]
        else:
            list = []
        list.append(assayNumber)
        assayLookup[markerID] = list

        line = fp1.readline()

    #
    # Process the records in input file 2 to create the specimen and results
    # output files.
    #
    specimenNumber = 0

    line = fp2.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        markerID = tokens[0]
        specimenLabel = tokens[1]
        genotypeID = tokens[2]
        age = tokens[3]
        ageNote = tokens[4]
        sex = tokens[5]
        fixation = tokens[6]
        embeddingMethod = tokens[7]
        hybridization = tokens[8]
        specimenNote = tokens[9]

        resultNumber = 0
        results = []

        #
        # Create an index into the token list that points to the first set
        # of results for the specimen and work through each set of results
        # one at a time.
        #
        i = FIRST_RESULT_IDX
        while i <= LAST_RESULT_IDX:
            strength = tokens[i]
            pattern = tokens[i+1]
            structure = tokens[i+2]
            theilerStage = tokens[i+3]
            resultNote = tokens[i+4]

            #
            # If there is a strength, save a set of results for this specimen
            # in the results list.
            #
            if strength != '':
                resultNumber += 1
                results.append((resultNumber, strength, pattern, structure,
                                theilerStage, resultNote))

            #
            # Advance the index to the next set of results.
            #
            i += 5

        #
        # Use the next unique specimen number for this specimen.
        #
        specimenNumber += 1

        #
        # Get the list of assay numbers that were used for this marker ID.
        #
        assayNumberList = assayLookup[markerID]

        #
        # Any records written to the specimen and results files need to be
        # duplicated for each assay that was created for the same marker ID.
        #
        for assayNumber in assayNumberList:

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
            # Process each tuple in the results list.
            #
            for r in results:

                (resultNumber, strength, pattern, structure, theilerStage, resultNote) = r

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
                                '' + '\n')

        line = fp2.readline()

    return 0


#
# Main
#
initialize()
buildProbeLookup()

openFiles()
process()
closeFiles()

sys.exit(0)
