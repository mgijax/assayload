#!/usr/local/bin/python

#
#  tr8270.py
###########################################################################
#
#  Purpose:
#
#      This script will create the input files for the GXD in situ load.
#
#  Usage:
#
#      tr8270.py
#
#  Env Vars:
#
#      STRUCTURE_TRANS_FILE
#      IMAGE_LIST_FIG_FILE
#      BEST_IMAGE_FILE
#      STR_PATT_TRANS_FILE
#      ASSAY_FILE
#      PROBEPREP_FILE
#      SPECIMEN_FILE
#      RESULTS_FILE
#      REFERENCE
#      CREATEDBY
#
#  Inputs:
#
#      StrengthPatternTrans.txt - Tab-delimited fields:
#
#          1) Marker MGI ID
#          2) Marker Symbol
#          3) Analysis ID
#          4) Probe ID
#          5) Probe Name
#          6) Specimen ID
#          7) Primer Name
#          8) Forward Primer
#          9) Reverse Primer
#          10) Strain ID
#          11) Strain MGI ID
#          12) Method ID
#          13) Accession Number
#          14) Expression level for structure 1
#          15) Pattern for structure 1
#          .
#          .  Repeat expression/pattern for each of the 100 structures.
#          .
#          .  NOTE: Not all structure numbers are found in the input file.
#          .
#          212) Expression level for structure 126
#          213) Pattern for structure 126
#          214) Probe Sequence
#
#      ImageListFigLabels.txt - Tab-delimited fields:
#
#          1) Marker MGI ID
#          2) Marker Symbol
#          3) Analysis ID
#          4) Probe ID
#          5) Probe Name
#          6) Specimen ID
#          7) Accession Number
#          8) Set number for image 1
#          9) Slide/section number for image 1
#          10) Image file name for image 1
#          11) Figure label for image 1
#          .
#          .  Repeat image data for a maximum of 24 images.
#          .
#          100) Set number for image 24
#          101) Slide/section number for image 24
#          102) Image file name for image 24
#          103) Figure label for image 24
#
#      BestImage.txt - Tab-delimited fields:
#
#          1) Marker MGI ID
#          2) Marker Symbol
#          3) Analysis ID
#          4) Probe ID
#          5) Probe Name
#          6) Specimen ID
#          7) Accession Number
#          8) Best image slide/section for structure 1
#          .
#          .  Repeat best image for each of the 100 structures.
#          .
#          .  NOTE: Not all structures have a best image.
#          .
#          107) Best image slide/section for structure 100
#
#      StructureTranslation.txt - Tab-delimited fields:
#
#          1) Eichele Structure Name
#          2) Eichele Structure Code
#          3) MGI Structure Name
#          4) MGI Structure Key
#          5) Note
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
#          6) Label Coverage
#          7) Visualized With
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
#
#  Modification History:
#
#  Date        SE   Change Description
#  ----------  ---  -------------------------------------------------------
#
#  09/17/2007  DBM  Initial development
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

# Column indexes in StructureTranslation.txt
STRUCTURE_CODE_INDEX = 1
MGI_STRUCTURE_NAME_INDEX = 2
EXPRESSION_NOTE_INDEX = 4

# Column indexes in ImageListFigLabels.txt
FIRST_SLIDE_SECTION_INDEX = 8

# Column indexes in BestImage.txt
FIRST_STRUCTURE_CODE_INDEX = 7

# Column indexes in StrengthPatternTrans.txt
MARKER_MGI_ID_INDEX = 0
ANALYSIS_ID_INDEX = 2
PROBE_NAME_INDEX = 4
SPECIMEN_ID_INDEX = 5
STRAIN_MGI_ID_INDEX = 10
FIRST_EXPRESSION_INDEX = 13
LAST_EXPRESSION_INDEX = 211

ASSAY_TYPE = "RNA in situ"
REPORTER_GENE = ""
ASSAY_NOTE = '''Cryosections of fresh frozen material were fixed in 4% paraformaldehyde for 20 min. before further processing. A tyramide-biotin/streptavidin amplification step was included in the in situ hybridization procedure.'''

PROBE_PREP_TYPE = "RNA"
PROBE_PREP_HYBRIDIZATION = "Antisense"
PROBE_PREP_LABELLED_WITH = "Digoxigenin"
PROBE_PREP_LABEL_COVERAGE = "Not Specified"
PROBE_PREP_VISUALIZED_WITH = "Alkaline phosphatase"

SPECIMEN_NUMBER = "1"
SPECIMEN_AGE = "embryonic day 14.5"
SPECIMEN_AGE_NOTE = ""
SPECIMEN_SEX = "Not Specified"
SPECIMEN_FIXATION = "Fresh Frozen"
SPECIMEN_EMBEDDING = "Cryosection"
SPECIMEN_HYBRIDIZATION = "section"
SPECIMEN_NOTE = ""

THEILER_STAGE = "22"
RESULT_NOTE = "Expression was %s in %s."

#
#  GLOBALS
#
strTransFile = os.environ['STRUCTURE_TRANS_FILE']
imageListFile = os.environ['IMAGE_LIST_FIG_FILE']
bestImageFile = os.environ['BEST_IMAGE_FILE']
strPattTransFile = os.environ['STR_PATT_TRANS_FILE']

assayFile = os.environ['ASSAY_FILE']
probePrepFile = os.environ['PROBEPREP_FILE']
specimenFile = os.environ['SPECIMEN_FILE']
resultsFile = os.environ['RESULTS_FILE']

jNumber = os.environ['REFERENCE']
createdBy = os.environ['CREATEDBY']


#
# Purpose: Create a dictionary for looking up the MGI structure name and
#          expression note for an Eichele structure code.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildStructureLookup ():
    global structureLookup

    structureLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(strTransFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + strTransFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining line of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = re.split('\t', line[:-1])
        eicheleStructureCode = tokens[STRUCTURE_CODE_INDEX]
        structureName = tokens[MGI_STRUCTURE_NAME_INDEX]
        expNote = tokens[EXPRESSION_NOTE_INDEX]

        #
        # If there is a structure name, add an entry to the dictionary.
        #
        if structureName != '':
            dict = {}
            dict['name'] = structureName
            dict['note'] = expNote
            structureLookup[eicheleStructureCode] = dict
            #print eicheleStructureCode + " = " + str(dict)
        line = fp.readline()

    fp.close()

    return


#
# Purpose: Create a dictionary for looking up the figure label for an
#          analysis ID and slide/section number.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildFigureLabelLookup ():
    global figureLabelLookup

    figureLabelLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(imageListFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + imageListFile + '\n')
        sys.exit(1)

    #
    # Skip the heading record and process the remaining line of the file.
    #
    line = fp.readline()
    line = fp.readline()
    while line:
        tokens = re.split('\t', line[:-1])
        analysisID = tokens[ANALYSIS_ID_INDEX]
        i = FIRST_SLIDE_SECTION_INDEX

        #
        # Check each set of image data from the input line.
        #
        while i < len(tokens):
            slideSection = tokens[i]
            figureLabel = tokens[i+2]

            #
            # If there is no figure label, skip ahead to the next one.
            #
            if figureLabel == '':
                i += 4
                continue

            #
            # Add the figure label to the dictionary for the analysis ID and
            # current slide/section.
            #
            figureLabelLookup[(analysisID,slideSection)] = figureLabel
            #print "|" + analysisID + "," + slideSection + "| = " + figureLabel
            i += 4

        line = fp.readline()

    fp.close()

    return


#
# Purpose: Create a dictionary for looking up the slide/section number for an
#          analysis ID and Eichele structure code.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildSlideSectionLookup ():
    global slideSectionLookup

    slideSectionLookup = {}

    #
    # Open the input file.
    #
    try:
        fp = open(bestImageFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + bestImageFile + '\n')
        sys.exit(1)

    #
    # Save the fields of the header record because they contain the Eichele
    # structure codes that correspond to the slide/section numbers in the
    # subsequent lines of data.
    #
    line = fp.readline()
    header = re.split('\t', line[:-1])

    line = fp.readline()
    while line:
        tokens = re.split('\t', line[:-1])
        analysisID = tokens[ANALYSIS_ID_INDEX]
        i = FIRST_STRUCTURE_CODE_INDEX

        #
        # Check each slide/section number from the input line.
        #
        while i < len(tokens):
            eicheleStructureCode = header[i]
            slideSection = tokens[i]

            #
            # If there is no slide/section number, skip ahead to the next one.
            #
            if slideSection == '':
                i += 1
                continue

            #
            # Add the figure label to the dictionary for the analysis ID and
            # current slide/section.
            #
            slideSectionLookup[(analysisID,eicheleStructureCode)] = slideSection
            #print "|" + analysisID + "," + eicheleStructureCode + "| = " + slideSection
            i += 1

        line = fp.readline()

    fp.close()

    return


#
# Purpose: Create a dictionary for looking up the MGI ID for a probe that
#          has been previously added for this project.
# Returns: Nothing
# Assumes: Nothing
# Effects: Sets global variable
# Throws: Nothing
#
def buildProbeLookup ():
    global probeLookup

    probeLookup = {}

    #
    # Get the segment type key for 'cDNA' probes. This will allow the query
    # below to exclude primers that have been added with the probes.
    #
    results = db.sql('select t._Term_key ' + \
                     'from VOC_Vocab v, VOC_Term t ' + \
                     'where v.name = "Segment Type" and ' + \
                           'v._Vocab_key = t._Vocab_key and ' + \
                           't.term = "cDNA"', 'auto')
    segmentTypeKey = results[0]['_Term_key']

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
    # Get all the probes and their MGI IDs.
    #
    results = db.sql('select p.name, a.accID ' + \
                     'from PRB_Probe p, PRB_Reference r, ACC_Accession a ' + \
                     'where p._SegmentType_key = ' + str(segmentTypeKey) + ' and ' + \
                           'p._Probe_key = r._Probe_key and ' + \
                           'r._Refs_key = ' + str(refKey) + ' and ' + \
                           'p._Probe_key = a._Object_key and ' + \
                           'a._LogicalDB_key = 1 and ' + \
                           'a._MGIType_key = 3', 'auto')
    for r in results:
        probeLookup[r['name']] = r['accID']
        #print r['name'] + ' = ' + r['accID']

    return


#
# Purpose: Open the files.
# Returns: Nothing
# Assumes: The names of the files are set in the environment.
# Effects: Sets global variables
# Throws: Nothing
#
def openFiles ():
    global fpStrPattTrans
    global fpProbePrepFile, fpAssayFile, fpSpecimenFile, fpResultsFile

    #
    # Open the input file.
    #
    try:
        fpStrPattTrans = open(strPattTransFile, 'r')
    except:
        sys.stderr.write('Cannot open input file: ' + strPattTransFile + '\n')
        sys.exit(1)

    #
    # Open the output files.
    #
    try:
        fpAssayFile = open(assayFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + assayFile + '\n')
        sys.exit(1)

    try:
        fpProbePrepFile = open(probePrepFile, 'w')
    except:
        sys.stderr.write('Cannot open output file: ' + probePrepFile + '\n')
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
    fpStrPattTrans.close()
    fpAssayFile.close()
    fpProbePrepFile.close()
    fpSpecimenFile.close()
    fpResultsFile.close()

    return


#
# Purpose: Read the expression/pattern data from the input file and create
#          the output files.
# Returns: 0 if successful, 1 for an error
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing
#
def process ():

    assayNumber = 0

    #
    # Save the fields of the first header record because they contain the
    # Eichele structure codes that correspond to the expression/pattern values
    # in the subsequent lines of data.
    #
    line = fpStrPattTrans.readline()
    header = re.split('\t', line[:-1])

    #
    # Skip the second header record and process the remaining records.
    #
    line = fpStrPattTrans.readline()
    line = fpStrPattTrans.readline()
    while line:
        tokens = re.split('\t', line[:-1])
        markerMGIID = tokens[MARKER_MGI_ID_INDEX]
        analysisID = tokens[ANALYSIS_ID_INDEX]
        probeName = tokens[PROBE_NAME_INDEX]
        specimenID = tokens[SPECIMEN_ID_INDEX]
        strainMGIID = tokens[STRAIN_MGI_ID_INDEX]

        #
        # Look up the probe MGI ID for the current probe name.
        #
        probeMGIID = probeLookup[probeName]

        #
        # Use the next unique assay number for this assay.
        #
        assayNumber += 1

        #
        # Write a record to the assay file.
        #
        fpAssayFile.write(str(assayNumber) + '\t' +
                      markerMGIID + '\t' +
                      jNumber + '\t' +
                      ASSAY_TYPE + '\t' +
                      REPORTER_GENE + '\t' +
                      ASSAY_NOTE + '\t' +
                      createdBy + '\n')

        #
        # Write a record to the probe prep file.
        #
        fpProbePrepFile.write(str(assayNumber) + '\t' +
                          probeMGIID + '\t' +
                          PROBE_PREP_TYPE + '\t' +
                          PROBE_PREP_HYBRIDIZATION + '\t' +
                          PROBE_PREP_LABELLED_WITH + '\t' +
                          PROBE_PREP_LABEL_COVERAGE + '\t' +
                          PROBE_PREP_VISUALIZED_WITH + '\n')

        #
        # Write a record to the specimen file.
        #
        fpSpecimenFile.write(str(assayNumber) + '\t' +
                         SPECIMEN_NUMBER + '\t' +
                         specimenID + '\t' +
                         strainMGIID + '\t' +
                         SPECIMEN_AGE + '\t' +
                         SPECIMEN_AGE_NOTE + '\t' +
                         SPECIMEN_SEX + '\t' +
                         SPECIMEN_FIXATION + '\t' +
                         SPECIMEN_EMBEDDING + '\t' +
                         SPECIMEN_HYBRIDIZATION  + '\t' +
                         SPECIMEN_NOTE + '\n')

        #
        # Create an index into the token list that points to the first
        # expression value and work through each set of expression/pattern
        # values.
        #
        resultsNumber = 0
        i = FIRST_EXPRESSION_INDEX
        while i <= LAST_EXPRESSION_INDEX:
            expression = tokens[i]
            pattern = tokens[i+1]
            eicheleStructureCode = header[i]
            eicheleStructureName = header[i+1]

            #
            # If there is no expression for the current structure, skip to
            # the next structure.
            #
            if expression == '':
                i += 2
                continue

            #
            # If the Eichele structure can be translated to a MGI structure,
            # generate the expression note using the expression and Eichele
            # structure name. Otherwise, skip to the next structure.
            #
            if structureLookup.has_key(eicheleStructureCode):
                structure = structureLookup[eicheleStructureCode]
                structureName = structure['name']
                expNote = structure['note']
                if expNote != '':
                    expNote = RESULT_NOTE % (expression.lower(), eicheleStructureName)
            else:
                i += 2
                continue

            #
            # If there is a best image slide/section for the current
            # structure, use it to look up the figure label that identifies
            # the image for the result.
            #
            if slideSectionLookup.has_key((analysisID, eicheleStructureCode)):
                slideSection = slideSectionLookup[(analysisID, eicheleStructureCode)]
                figureLabel = figureLabelLookup[(analysisID, slideSection)]
            else:
                figureLabel = ''

            #
            # Use the next unique result number for this assay result.
            #
            resultsNumber += 1

            #
            # Write a record to the results file.
            #
            fpResultsFile.write(str(assayNumber) + '\t' +
                            SPECIMEN_NUMBER + '\t' +
                            str(resultsNumber) + '\t' +
                            expression + '\t' +
                            pattern + '\t' +
                            structureName + '\t' +
                            THEILER_STAGE + '\t' +
                            expNote + '\t' +
                            figureLabel + '\n')

            #
            # Advance the index to the next set of expression/pattern values.
            #
            i += 2

        line = fpStrPattTrans.readline()

    return 0


#
# Main
#
buildStructureLookup()
buildFigureLabelLookup()
buildSlideSectionLookup()
buildProbeLookup()
openFiles()
process()
closeFiles()

sys.exit(0)
