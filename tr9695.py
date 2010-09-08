#!/usr/local/bin/python

#
#  tr9695.py
###########################################################################
#
#  Purpose:
#
#      This script will create the input files for the GXD in situ load.
#
#  Usage:
#
#      tr9695.py
#
#  Env Vars:
#
#      LOADFILE1
#      LOADFILE2
#      GENOTYPE_LOOKUP
#      STRENGTH_LOOKUP
#      PATTERN_LOOKUP
#      SPEC_NOTES_LOOKUP
#      RES_NOTES_LOOKUP
#      PROBEPREP_FILE
#      ASSAY_FILE
#      SPECIMEN_FILE
#      RESULTS_FILE
#      REFERENCE
#      CREATEDBY
#
#  Inputs:
#
#      LoadFile1.txt - Tab-delimited fields:
#
#          1) Probe MGI ID
#          2) Marker MGI ID
#          3) MGI Symbol       # Not used
#          4) Annotation ID
#          5) Assay Note
#
#      LoadFile2.txt - Tab-delimited fields:
#
#          1) Annotation ID
#          2) Specimen Label
#          3) Figure Label                # Not used - same as specimen label
#          4) Image File                  # Not used
#          5) Cmpt Name (for result 1)    # Not used
#          6) Stage (for result 1)
#          7) EMAP ID (for result 1)
#          8) Strength (for result 1)
#          9) Pattern (for result 1)
#          10) Note (for result 1)
#          .
#          .  There are up to 49 sets of result data (repeat columns 5-10).
#          .
#          293) Cmpt Name (for result 49)    # Not used
#          294) Stage (for result 49)
#          295) EMAP ID (for result 49)
#          296) Strength (for result 49)
#          297) Pattern (for result 49)
#          298) Note (for result 49)
#
#      GenotypeLookupTable.txt - Tab-delimited fields:
#
#          1) Annotation ID
#          2) Genotype MGI ID
#
#      StrengthLookupTable.txt - Tab-delimited fields:
#
#          1) Eurexpress Strength Term
#          2) GXD Strength Term
#
#      PatternLookupTable.txt - Tab-delimited fields:
#
#          1) Eurexpress Pattern Term
#          2) GXD Pattern Term
#
#      SpecimenNotesLookupTable.txt - Tab-delimited fields:
#
#          1) Specimen Label
#          2) Specimen Note
#
#      ResultNotesLookupTable.txt - Tab-delimited fields:
#
#          1) Specimen Label
#          2) EMAP ID
#          3) Result Note (blank if none)
#          4) MGI Structure Key (pipe-separated structure keys for any
#                                additional structures to add for this
#                                specimen - blank if none).
#          5) Strength (for any additional structures)
#          6) Pattern (for any additional structures)
#          7) Result Note (for any additional structures - blank if none)
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

# Column indexes in LoadFile1
F1_PROBE_MGI_ID_IDX = 0
F1_MARKER_MGI_ID_IDX = 1
F1_ANNOT_ID_IDX = 3
F1_ASSAY_NOTE_IDX = 4

# Column indexes in LoadFile2
F2_ANNOT_ID_IDX = 0
F2_SPECIMEN_LABEL_IDX = 1
F2_FIRST_RESULT_IDX = 4      # Index where the first result starts
F2_LAST_RESULT_IDX = 292     # Index where the last result starts

PROBE_PREP_TYPE = "RNA"
PROBE_PREP_HYBRIDIZATION = "Antisense"
PROBE_PREP_LABELLED_WITH = "digoxigenin"
PROBE_PREP_VISUALIZED_WITH = "Alkaline phosphatase"

ASSAY_TYPE = "RNA in situ"
REPORTER_GENE = ""

SPECIMEN_AGE = "embryonic day 14.5"
SPECIMEN_AGE_NOTE = ""
SPECIMEN_SEX = "Not Specified"
SPECIMEN_FIXATION = "Fresh Frozen"
SPECIMEN_EMBEDDING = "Cryosection"
SPECIMEN_HYBRIDIZATION = "section"

THEILER_STAGE = "23"

#
#  GLOBALS
#
loadFile1 = os.environ['LOADFILE1']
loadFile2 = os.environ['LOADFILE2']
genotypeFile = os.environ['GENOTYPE_LOOKUP']
strengthFile = os.environ['STRENGTH_LOOKUP']
patternFile = os.environ['PATTERN_LOOKUP']
specimenNoteFile = os.environ['SPEC_NOTES_LOOKUP']
resultNoteFile = os.environ['RES_NOTES_LOOKUP']

probePrepFile = os.environ['PROBEPREP_FILE']
assayFile = os.environ['ASSAY_FILE']
specimenFile = os.environ['SPECIMEN_FILE']
resultsFile = os.environ['RESULTS_FILE']

jNumber = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']


#
# Purpose: Use the genotype lookup file to build a dictionary for looking
#          up the genotype ID for a given annotation ID.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildGenotypeLookup ():
    global genotypeLookup

    genotypeLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(genotypeFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + genotypeFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining lines of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        annotID = tokens[0]
        genotypeID = tokens[1]

        genotypeLookup[annotID] = genotypeID
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Use the strength lookup file to build a dictionary for looking
#          up the GXD strength term for a given Eurexpress strength term.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildStrengthLookup ():
    global strengthLookup

    strengthLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(strengthFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + strengthFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining lines of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        eurexpressTerm = tokens[0].lower()
        gxdTerm = tokens[1]

        strengthLookup[eurexpressTerm] = gxdTerm
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Use the pattern lookup file to build a dictionary for looking
#          up the GXD pattern term for a given Eurexpress pattern term.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildPatternLookup ():
    global patternLookup

    patternLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(patternFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + patternFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining lines of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        eurexpressTerm = tokens[0].lower()
        gxdTerm = tokens[1]

        patternLookup[eurexpressTerm] = gxdTerm
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Use the specimen note lookup file to build a dictionary for looking
#          up the specimen note for a given specimen label.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildSpecimenNoteLookup ():
    global specimenNoteLookup

    specimenNoteLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(specimenNoteFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + specimenNoteFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining lines of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        specimenLabel = tokens[0]
        specimenNote = tokens[1]

        specimenNoteLookup[specimenLabel] = specimenNote
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Use the result note lookup file to build a dictionary for looking
#          up the result note and additional structure info for a given
#          specimen label and EMAP ID.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildResultNoteLookup ():
    global resultNoteLookup

    resultNoteLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(resultNoteFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + resultNoteFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining lines of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        specimenLabel = tokens[0]
        emapID = tokens[1]
        resultNote = tokens[2]
        if tokens[3] != '':
            addStructureKeyList = string.split(tokens[3], '|')
        else:
            addStructureKeyList = []
        addStrength = tokens[4]
        addPattern = tokens[5]
        addResultNote = tokens[6]

        key = (specimenLabel, emapID)
        value = (resultNote, addStructureKeyList, addStrength,
                 addPattern, addResultNote)

        resultNoteLookup[key] = value
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Build two structure lookups, one for looking up the structure name
#          for a given EMAP ID and one for looking up the structure name for
#          a given MGI structure key. 
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variables
# Throws: Nothing
#
def buildStructureLookups ():
    global emapStructureLookup, mgiStructureLookup

    emapStructureLookup = {}
    mgiStructureLookup = {}

    #
    # Get the structure names for each EMAP ID.
    #
    results = db.sql('select s.edinburghKey, s.printName ' + \
                     'from GXD_Structure s, GXD_TheilerStage t ' + \
                     'where s.edinburghKey is not null and ' + \
                           's._Stage_key = t._Stage_key and ' + \
                           't.stage = %s' % (THEILER_STAGE),'auto')

    #
    # Add each key/name to the dictionary.
    #
    for r in results:
        emapStructureLookup[str(r['edinburghKey'])] = r['printName']

    #
    # Get the structure names for each MGI structure key.
    #
    results = db.sql('select s._Structure_key, s.printName ' + \
                     'from GXD_Structure s, GXD_TheilerStage t ' + \
                     'where s._Stage_key = t._Stage_key and ' + \
                           't.stage = %s' % (THEILER_STAGE),'auto')

    #
    # Add each key/name to the dictionary.
    #
    for r in results:
        mgiStructureLookup[str(r['_Structure_key'])] = r['printName']

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
        fp1 = open(loadFile1, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + loadFile1 + '\n')
        sys.exit(1)

    try:
        fp2 = open(loadFile2, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + loadFile2 + '\n')
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
# Purpose: Use the records in LoadFile1 to create the probe prep and assay
#          output files. Use the records in LoadFile2 to create the specimen
#          and results output files.
# Returns: 0 if successful, 1 for an error
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def process ():

    #
    # Process the records in LoadFile1 to create the probe prep and assay
    # output files.
    #
    assayNumber = 0
    assayLookup = {}

    #
    # Skip the header record and process the remaining records.
    #
    line = fp1.readline()
    line = fp1.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        probeMGIID = tokens[F1_PROBE_MGI_ID_IDX]
        markerMGIIDList = string.split(tokens[F1_MARKER_MGI_ID_IDX], '|')
        annotID = tokens[F1_ANNOT_ID_IDX]
        assayNote = tokens[F1_ASSAY_NOTE_IDX]

        #
        # There should be a new assay created for each marker.
        #
        for markerMGIID in markerMGIIDList:

            #
            # Use the next unique assay number for this assay.
            #
            assayNumber += 1

            #
            # Write a record to the probe prep file.
            #
            fpProbePrepFile.write(str(assayNumber) + '\t' +
                              probeMGIID + '\t' +
                              PROBE_PREP_TYPE + '\t' +
                              PROBE_PREP_HYBRIDIZATION + '\t' +
                              PROBE_PREP_LABELLED_WITH + '\t' +
                              PROBE_PREP_VISUALIZED_WITH + '\n')

            #
            # Write a record to the assay file.
            #
            fpAssayFile.write(str(assayNumber) + '\t' +
                          markerMGIID + '\t' +
                          jNumber + '\t' +
                          ASSAY_TYPE + '\t' +
                          REPORTER_GENE + '\t' +
                          assayNote + '\t' +
                          createdBy + '\n')
            #
            # Maintain a lookup for finding all assay numbers that are used
            # for each annotation ID. This will be needed when processing
            # LoadFile2.
            #
            if assayLookup.has_key(annotID):
                list = assayLookup[annotID]
            else:
                list = []
            list.append(assayNumber)
            assayLookup[annotID] = list

        line = fp1.readline()

    #
    # Process the records in LoadFile2 to create the specimen and results
    # output files.
    #
    specimenNumber = 0

    #
    # Skip the header record and process the remaining records.
    #
    line = fp2.readline()
    line = fp2.readline()
    while line:
        tokens = string.split(line[:-1], '\t')
        annotID = tokens[F2_ANNOT_ID_IDX]
        specimenLabel = tokens[F2_SPECIMEN_LABEL_IDX]

        resultNumber = 0
        results = []

        #
        # Look up the genotype ID for the annotation ID.
        #
        if genotypeLookup.has_key(annotID):
            genotypeID = genotypeLookup[annotID]
        else:
            sys.stderr.write('Genotype lookup failed:' +
                             ' annotID = ' + annotID + '\n')
            sys.exit(1)

        #
        # Look up the specimen note for the specimen.
        #
        if specimenNoteLookup.has_key(specimenLabel):
            specimenNote = specimenNoteLookup[specimenLabel]
        else:
            specimenNote = ''

        #
        # Create an index into the token list that points to the first set
        # of results for the specimen and work through each set of results
        # one at a time.
        #
        i = F2_FIRST_RESULT_IDX
        while i <= F2_LAST_RESULT_IDX:
            emapID = tokens[i+2]
            strength = tokens[i+3]
            pattern = tokens[i+4]
            resultNote = tokens[i+5]

            #
            # If there is no EMAP ID, there must not be any more results
            # for this specimen.
            #
            if emapID == '':
                break

            #
            # If the strength is "not examined", do not load this structure.
            # Just advance the index to the next set of results.
            #
            if strength == 'not examined':
                i += 6
                continue

            #
            # Look up the strength for the structure.
            #
            if strengthLookup.has_key(strength.lower()):
                strength = strengthLookup[strength.lower()]
            else:
                sys.stderr.write('Strength lookup failed:' +
                                 ' strength = ' + strength + '\n')
                sys.exit(1)

            #
            # Look up the pattern for the structure.
            #
            if patternLookup.has_key(pattern.lower()):
                pattern = patternLookup[pattern.lower()]
            else:
                sys.stderr.write('Pattern lookup failed:' +
                                 ' pattern = ' + pattern + '\n')
                sys.exit(1)

            #
            # Look up the structure name for the EMAP ID.
            #
            if emapStructureLookup.has_key(emapID):
                structure = emapStructureLookup[emapID]
            else:
                sys.stderr.write('Structure lookup failed:' +
                                 ' EMAP ID = ' + emapID + '\n')
                sys.exit(1)

            #
            # If there is a result note, look up the note that should be
            # used, along with any other structures that should be added.
            # If the note is missing from the lookup, the note file needs
            # to be corrected.
            #
            if resultNote != '':
                if resultNoteLookup.has_key((specimenLabel, emapID)):
                    (resultNote, addStructureKeyList, addStrength, addPattern, addResultNote) = resultNoteLookup[(specimenLabel, emapID)]
                else:
                    sys.stderr.write('Result note lookup failed:' +
                                     ' specimen label = ' + specimenLabel +
                                     ', EMAP ID = ' + emapID + '\n')
                    sys.exit(1)

            #
            # If there is no result note, there won't be any additional
            # structures to add either.
            #
            else:
                addStructureKeyList = []

            #
            # Use the next unique result number for this assay result.
            #
            resultNumber += 1

            #
            # Add an entry to the results list.
            #
            results.append((resultNumber, strength, pattern, structure,
                            resultNote))

            #
            # If there are any additional structures that need to be added,
            # look up the strength and pattern. These will be the same for
            # each additional structure. 
            #
            if len(addStructureKeyList) > 0:

                #
                # Look up the strength for the structure.
                #
                if strengthLookup.has_key(addStrength.lower()):
                    strength = strengthLookup[addStrength.lower()]
                else:
                    sys.stderr.write('Strength lookup failed:' +
                                     ' strength = ' + addStrength + '\n')
                    sys.exit(1)

                #
                # Look up the pattern for the structure.
                #
                if patternLookup.has_key(addPattern.lower()):
                    pattern = patternLookup[addPattern.lower()]
                else:
                    sys.stderr.write('Pattern lookup failed:' +
                                     ' pattern = ' + addPattern + '\n')
                    sys.exit(1)

                #
                # Add an entry to the results list for each additional
                # structure.
                #
                for structureKey in addStructureKeyList:

                    #
                    # Look up the structure name for the MGI structure key.
                    #
                    if mgiStructureLookup.has_key(structureKey):
                        structure = mgiStructureLookup[structureKey]
                    else:
                        sys.stderr.write('Structure lookup failed:' +
                                         ' MGI Key = ' + structureKey + '\n')
                        sys.exit(1)

                    #
                    # Use the next unique result number for this assay result.
                    #
                    resultNumber += 1

                    #
                    # Add an entry to the results list.
                    #
                    results.append((resultNumber, strength, pattern, structure, addResultNote))

            #
            # Advance the index to the next set of results.
            #
            i += 6

        #
        # Use the next unique specimen number for this specimen.
        #
        specimenNumber += 1

        #
        # Get the list of assay numbers for this annotID.
        #
        assayNumberList = assayLookup[annotID]

        #
        # Any records written to the specimen and results files need to be
        # duplicated for each assay that was created for the same
        # annotation ID.
        #
        for assayNumber in assayNumberList:

            #
            # Write a record to the specimen file.
            #
            fpSpecimenFile.write(str(assayNumber) + '\t' +
                             str(specimenNumber) + '\t' +
                             specimenLabel + '\t' +
                             genotypeID + '\t' +
                             SPECIMEN_AGE + '\t' +
                             SPECIMEN_AGE_NOTE + '\t' +
                             SPECIMEN_SEX + '\t' +
                             SPECIMEN_FIXATION + '\t' +
                             SPECIMEN_EMBEDDING + '\t' +
                             SPECIMEN_HYBRIDIZATION  + '\t' +
                             specimenNote + '\n')

            #
            # Process each tuple in the results list.
            #
            for r in results:

                (resultNumber, strength, pattern, structure, resultNote) = r

                #
                # Write a record to the results file.
                #
                fpResultsFile.write(str(assayNumber) + '\t' +
                                str(specimenNumber) + '\t' +
                                str(resultNumber) + '\t' +
                                strength + '\t' +
                                pattern + '\t' +
                                structure + '\t' +
                                THEILER_STAGE + '\t' +
                                resultNote + '\t' +
                                specimenLabel + '\n')

        line = fp2.readline()

    return 0


#
# Main
#
buildGenotypeLookup()
buildStrengthLookup()
buildPatternLookup()
buildSpecimenNoteLookup()
buildResultNoteLookup()
buildStructureLookups()

openFiles()
process()
closeFiles()

sys.exit(0)
