#!/usr/local/bin/python

#
#  tr6739insitu.py
###########################################################################
#
#  Purpose:
#
#      This script will use the insitu data in the retina_data input file
#      and the translations in the tissue_translation input file to create
#      a set of output files that will be used as input to the insituload.
#
#  Usage:
#
#      tr6739insitu.py
#
#  Env Vars:  None
#
#  Inputs:
#
#      The input file Retina_data.txt contains the following
#      tab-delimited fields:
#
#      1) SAGE_tag
#      2) Gene_author_description
#      3) MGI_Symbol
#      4) Probe_Accession_Number
#      5) MGI_Probe_Accession_ID
#      6) Probe_Name
#      7) Clone_Library
#      8 - 15) Fields containing insitu data for E12
#      16 - 24) Fields containing insitu data for E14
#      25 - 33) Fields containing insitu data for E16
#      34 - 45) Fields containing insitu data for E18
#      46 - 57) Fields containing insitu data for P0-P1
#      58 - 69) Fields containing insitu data for P2-P3
#      70 - 81) Fields containing insitu data for P4-P5
#      82 - 96) Fields containing insitu data for P6-P7
#      97 - 111) Fields containing insitu data for P8-P10
#      112 - 126) Fields containing insitu data for Adult
#      127 - 144) Fields containing insitu data for CNS
#      145 - 162) Fields containing insitu data for Embryo
#
#      The input file Tissue_Translation.txt contains 24 tab-delimited
#      fields. There are 6 ages (E12, E14, E16, E18, P0-P5 and P6-Adult)
#      and the following fields exist for each age:
#
#      1) The code word used by the author for a tissue at the specific age.
#      2) The structure from the GXD anatomical dictionary that should
#         be used to annotate the data for the author's tissue code.
#      3) An indicator as to whether a note is required in the retina
#         data file for the given author's tissue code. It will be "yes"
#         for required or blank if not required.
#      4) The GXD pattern of expression to be coded for the given author's
#         tissue code.
#
#  Outputs:
#
#      4 tab-delimited files to be used as input to the insituload:
#
#      1) In_Situ_probeprep.txt
#      2) In_Situ_assay.txt
#      3) In_Situ_specimen.txt
#      4) In_Situ_results.txt
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
#  08/18/2005  DBM  Initial development
#
###########################################################################

import sys
import os
import string
import db
import reportlib
import accessionlib
import tr6739lib

#
#  CONSTANTS
#

TAB = '\t'
CRT = '\n'
NULL = ''

JNUMBER = 'J:93300'

# Probe prep constants
PROBE_PREP_TYPE = 'RNA'
SENSE = 'Antisense'
LABELED_WITH = 'Digoxigenin'
LABEL_COVERAGE = 'Not Specified'
VISUALIZATED_WITH = 'Alkaline phosphatase'

# Assay constants
ASSAY_TYPE = 'RNA in situ'
REPORTER_GENE = NULL
ASSAY_NOTE = 'This data is from Supplemental Table 5, Gene column, ' + \
             'row: %s. After cryosectioning, tissue was fixed for ' + \
             '10\' in 4%% paraformaldehyde prior to processing.'

# Specimen constants
GENOTYPE1 = 'MGI:2166522'    # Adult
GENOTYPE2 = 'MGI:2166348'    # All other
AGE_NOTE = NULL
SEX1 = 'Male'             # Adult
SEX2 = 'Not Specified'    # All other
FIXATION = "Fresh Frozen"
EMBEDDING_METHOD = "Cryosection"
HYBRIDIZATION = "section"
SPECIMEN_NOTE = NULL


#
#  FUNCTIONS
#

# Purpose: Perform initialization for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def initialize():
    global createdBy
    global fpProbePrep, fpAssay, fpSpecimen, fpResults

    print 'Perform initialization'
    sys.stdout.flush()

    #
    #  Initialize global variables.
    #
    dataDir = os.environ['INSITUDATADIR']
    createdBy = os.environ['CREATEDBY']
    dbServer = os.environ['DBSERVER']
    dbName = os.environ['DBNAME']
    dbUser = os.environ['DBUSER']
    dbPasswordFile = os.environ['DBPASSWORDFILE']
    dbPassword = string.strip(open(dbPasswordFile,'r').readline())

    #
    #  Set up a connection to the database.
    #
    db.useOneConnection(1)
    db.set_sqlLogin(dbUser, dbPassword, dbServer, dbName)

    #
    #  Open the 4 output files.
    #
    probePrepFile = dataDir + '/In_Situ_probeprep.txt'
    assayFile = dataDir + '/In_Situ_assay.txt'
    specimenFile = dataDir + '/In_Situ_specimen.txt'
    resultsFile = dataDir + '/In_Situ_results.txt'

    try:
        fpProbePrep = open(probePrepFile,'w')
    except:
        sys.stderr.write('Could not open probe prep file: ' + probePrepFile)
        exit(1)

    try:
        fpAssay = open(assayFile,'w')
    except:
        sys.stderr.write('Could not open assay file: ' + assayFile)
        exit(1)

    try:
        fpSpecimen = open(specimenFile,'w')
    except:
        sys.stderr.write('Could not open specimen file: ' + specimenFile)
        exit(1)

    try:
        fpResults = open(resultsFile,'w')
    except:
        sys.stderr.write('Could not open results file: ' + resultsFile)
        exit(1)

    return


# Purpose: Perform cleanup steps for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def finalize():
    global fpProbePrep, fpAssay, fpSpecimen, fpResults

    db.useOneConnection(0)

    #
    #  Close the 4 output files.
    #
    fpProbePrep.close()
    fpAssay.close()
    fpSpecimen.close()
    fpResults.close()

    return


# Purpose: Create a dictionary of the line numbers that have QC errors, so
#          they can be skipped during processing of the input file.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def buildQCLineNumDict():
    global qcLineNum

    print 'Create a dictionary of input line numbers with discrepancies'
    sys.stdout.flush()

    results = db.sql('select lineNum from tempdb..TMP_QC', 'auto')

    qcLineNum = {}
    for r in results:
        qcLineNum[r['lineNum']] = r['lineNum']

    return


# Purpose: Create a dictionary of the probes and markers to be used when
#          processing each line of the input file.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def buildProbeMarkerDict():
    global probeMarker

    print 'Create a dictionary of probe and marker MGI IDs'
    sys.stdout.flush()

    results = db.sql('select lineNum, markerAccID, probeAccID ' + \
                     'from tempdb..TMP_ProbeMarker ' + \
                     'where probeAccID is not null ' + \
                     'union ' + \
                     'select tp.lineNum, tp.markerAccID, ' + \
                            'a.accID "probeAccID" ' + \
                     'from tempdb..TMP_ProbeMarker tp, ' + \
                          'PRB_Probe p, ACC_Accession a ' + \
                     'where tp.probeName is not null and ' + \
                           'tp.probeName = p.name and ' + \
                           'p._Probe_key = a._Object_key and ' + \
                           'a._MGIType_key = 3 and ' + \
                           'a._LogicalDB_key = 1 and ' + \
                           'a.prefixPart = "MGI:" and ' + \
                           'a.preferred = 1', 'auto')

    probeMarker = {}
    for r in results:
        probeMarker[r['lineNum']] = {'marker':r['markerAccID'],
                                     'probe':r['probeAccID']}

    return


# Purpose: Process all of the data from the retina data input file by
#          using the DataParser class to parse the file.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def processData():
    global createdBy
    global fpProbePrep, fpAssay, fpSpecimen, fpResults
    global qcLineNum, probeMarker

    print 'Process each input record that does not have a discrepancy'
    sys.stdout.flush()

    assayKey = 1
    specimenKey = 0
    resultKey = 0
    lastLabel = ''
    createdSpecimen = 0

    #
    #  Get a parser object that knows how to return tissue data from the
    #  input file.
    #
    parser = tr6739lib.DataParser()

    #
    #  Get the first set of tissue data.
    #
    data, translation, last = parser.nextTissue()

    #
    #  Continue processing tissues until there are none left.
    #
    while data != None:
        lineNum = data['lineNum']

        #
        #  If the current line number is one that was flagged with a QC
        #  error, skip over all of the tissues in the current record.
        #
        if qcLineNum.has_key(lineNum):
            while last != 1:
                data, translation, last = parser.nextTissue()
            data, translation, last = parser.nextTissue()
            continue

        #
        #  Extract the data from the parser.
        #
        geneDesc = parser.getGeneDesc()
        symbol = parser.getSymbol()
        tissue = translation['gxd_tissue']
        note = data['note']
        strength = data['strength']
        pattern = translation['pattern']
        label = data['label']
        age = data['age']
        stage = data['stage']
        images = data['images']

        if pattern == '':
            pattern = 'Not Specified'

        #
        #  If the tissue translates to 'no tissue', set the flag to
        #  skip to the next one.
        #
        if tissue == 'no tissue':
            skipTissue = 1
        else:
            skipTissue = 0

        # cycle specimen key from 1-12
        if label != lastLabel:
            specimenKey = specimenKey + 1
            createdSpecimen = 0
            lastLabel = label

        # cycle result key from 1-40
        resultKey = resultKey + 1

        #
        #  If the tissue is not being skipped, determine what needs to
        #  be written of the output files.
        #
        if skipTissue == 0:
            #
            #  Create a specimen record if it hasn't been done yet.
            #
            if createdSpecimen == 0:
                #
                #  Set values needed for the specimen record.
                #
                if label == 'Adult':
                    genotype = GENOTYPE1
                    sex = SEX1
                else:
                    genotype = GENOTYPE2
                    sex = SEX2

                #
                #  Write a record to the specimen file.
                #
                fpSpecimen.write(str(assayKey) + TAB +
                                 str(specimenKey) + TAB +
                                 label + TAB +
                                 genotype + TAB +
                                 age + TAB +
                                 AGE_NOTE + TAB +
                                 sex + TAB +
                                 FIXATION + TAB +
                                 EMBEDDING_METHOD + TAB +
                                 HYBRIDIZATION + TAB +
                                 SPECIMEN_NOTE + CRT)

                createdSpecimen = 1

            #
            #  Set values needed for the results record.
            #
            strengthValue = string.atof(strength)
            if strengthValue == 0.0:
                strength = 'Absent'
                pattern = 'Not Applicable'
            elif strengthValue > 0.0 and strengthValue <= 0.5:
                strength = 'Weak'
            elif strengthValue > 0.5 and strengthValue < 5.0:
                strength = 'Present'
            else:
                strength = 'Strong'

            #
            #  Write a record to the results file.
            #
            fpResults.write(str(assayKey) + TAB +
                            str(specimenKey) + TAB +
                            str(resultKey) + TAB +
                            strength + TAB +
                            pattern + TAB +
                            tissue + TAB +
                            stage + TAB +
                            note + TAB +
                            images + CRT)

        #
        #  If there are no more tissues to process for this input record,
        #  write a record to the probe prep and assay files.
        #
        if last == 1:
            fpProbePrep.write(str(assayKey) + TAB +
                              probeMarker[lineNum]['probe'] + TAB +
                              PROBE_PREP_TYPE + TAB +
                              SENSE + TAB +
                              LABELED_WITH + TAB +
                              LABEL_COVERAGE + TAB +
                              VISUALIZATED_WITH + CRT)

            fpAssay.write(str(assayKey) + TAB +
                          probeMarker[lineNum]['marker'] + TAB +
                          JNUMBER + TAB +
                          ASSAY_TYPE + TAB +
                          REPORTER_GENE + TAB +
                          ASSAY_NOTE % (geneDesc) + TAB +
                          createdBy + CRT)

            #
            #  Set the flags for the next input record.
            #
            assayKey = assayKey + 1
            specimenKey = 0
            resultKey = 0
            lastLabel = ''

        #
        #  Get the next set of tissue data.
        #
        data, translation, last = parser.nextTissue()

    #
    #  Destroy the parser object.
    #
    parser = None

    return


#
#  MAIN
#
initialize()
buildQCLineNumDict()
buildProbeMarkerDict()
processData()
finalize()

