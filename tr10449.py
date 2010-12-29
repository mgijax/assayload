#!/usr/local/bin/python

#
#  tr10449.py
###########################################################################
#
#  Purpose:
#
#      This script will create the input files for the GXD in situ load.
#
#  Usage:
#
#      tr10449.py
#
#  Env Vars:
#
#      ASSAY_DATA
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
#          5) Visualized With (all are "Alkaline phosphatase")
#          6) Reference (all are "J:148410")
#          7) Assay Type (all are "RNA in situ")
#          8) Reporter Gene (all are null)
#          9) Assay Note
#          10) Created By (all are "cms")
#          11) Specimen Label
#          12) Genotype ("CD-1" or "Quackenbush")
#          13) Age (either embryonic day 12.5 or 13.5)
#          14) Age Note
#          15) Sex (all are "Male")
#          16) Fixation (all are "4% Paraformaldehyde")
#          17) Embedding Method (all are "Not Applicable")
#          18) Hybridization (all are "whole mount")
#          19) Specimen Note
#          20) Strength (for result 1)
#          21) Pattern (for result 1)
#          22) MGI Structure (for result 1)
#          23) Theiler Stage (for result 1)
#          24) Note (for result 1)
#          .
#          .  There are 7 sets of result data (repeat columns 20-24).
#          .
#          50) Strength (for result 7)
#          51) Pattern (for result 7)
#          52) MGI Structure (for result 7)
#          53) Theiler Stage (for result 7)
#          54) Note (for result 7)
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
FIRST_RESULT_IDX = 19     # Index where the first result starts
LAST_RESULT_IDX = 49      # Index where the last result starts

# Date the probes were created
PROBE_CREATION_DATE = "12/15/2010"

#
#  GLOBALS
#
loadFile = os.environ['ASSAY_DATA']

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
    global fp
    global fpProbePrepFile, fpAssayFile, fpSpecimenFile, fpResultsFile

    #
    # Open the input file.
    #
    try:
        fp = open(loadFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + loadFile + '\n')
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
    fp.close()
    fpProbePrepFile.close()
    fpAssayFile.close()
    fpSpecimenFile.close()
    fpResultsFile.close()

    return


#
# Purpose: Use the records in the input file to create the probe prep,
#          assay, specimen and results output files.
# Returns: 0 if successful, 1 for an error
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def process ():

    #
    # Initialize assay and specimen counters.
    #
    assayNumber = 0
    specimenNumber = 1   # Always 1 for this load.

    #
    # Process the input records one at a time.
    #
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        markerID = tokens[0]
        probePrepType = tokens[1]
        hybridization1 = tokens[2]
        labelledWith = tokens[3]
        visualizedWith = tokens[4]
        # skip reference field (tokens[5]) and use the config setting.
        assayType = tokens[6]
        reportedGene = tokens[7]
        assayNote = tokens[8]
        # skip createdBy field (tokens[9]) and use the config setting.
        specimenLabel = tokens[10]
        genotypeID = tokens[11]
        age = tokens[12]
        ageNote = tokens[13]
        sex = tokens[14]
        fixation = tokens[15]
        embeddingMethod = tokens[16]
        hybridization2 = tokens[17]
        specimenNote = tokens[18]

        #
        # Find the MGI ID of the probe that was create for this marker.
        #
        probeID = probeLookup[markerID]

        #
        # Use the next unique assay number for each output file.
        #
        assayNumber += 1

        #
        # Write a record to the probe prep file.
        #
        fpProbePrepFile.write(str(assayNumber) + '\t' +
                          probeID + '\t' +
                          probePrepType + '\t' +
                          hybridization1 + '\t' +
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
                         hybridization2  + '\t' +
                         specimenNote + '\n')

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

        line = fp.readline()

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
