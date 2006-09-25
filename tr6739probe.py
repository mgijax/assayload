#!/usr/local/bin/python

#
#  tr6739probe.py
###########################################################################
#
#  Purpose:
#
#      This script will use the insitu data in the retina_data input file
#      to update existing MGI probes and create an output file that will be
#      used as input to the probeload product to create new MGI probes.
#      The output file will contain data for 2 types of probes: those that
#      have known clone libraries (type 1) and those that have no source
#      information (type 2).
#
#  Usage:
#
#      tr6739probe.py
#
#  Env Vars:  None
#
#  Inputs:
#
#      The input file (Retina_data.txt) that contains the following
#      tab-delimited fields:
#
#      1) SAGE_tag
#      2) Gene_author_description
#      3) MGI_Symbol
#      4) Probe_Accession_Number
#      5) MGI_Probe_Accession_ID
#      6) Probe_Name
#      7) Clone_Library
#      8 - 162) These are additional fields that not used by this script.
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
#      A tab-delimited file (probe.txt) with the following fields to be
#      used as input to the probeload:
#
#      1) Probe Name
#      2) Reference (J:93300)
#      3) Source Name
#      4) Organism
#      5) Strain
#      6) Tissue
#      7) Gender
#      8) Cell Line
#      9) Age
#      10) Vector Type
#      11) Segment Type
#      12) Region Covered
#      13) Insert Site
#      14) Insert Size
#      15) MGI Marker
#      16) Relationship
#      17) Sequence ID
#      18) Notes
#      19) Created By
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

J_NUMBER = 'J:93300'
ORGANISM = 'mouse, laboratory'
STRAIN = 'Not Specified'
TISSUE = 'Not Specified'
GENDER = 'Not Specified'
CELLLINE = 'Not Specified'
AGE = 'Not Specified'
VECTOR_TYPE = 'Not Specified'
SEGMENT_TYPE = 'cDNA'
INSERT_SITE = 'Not Specified'
RELATIONSHIP = 'H'
LOGICAL_DB = 'Sequence DB'


#
#  FUNCTIONS
#

# Purpose: Perform initialization for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def initialize():
    global createdBy, refKey
    global dbServer, dbUser, dbPasswordFile
    global tmpPrbFile, tmpPrbMrkFile
    global fpQC, fpTmpPrb, fpTmpPrbMrk, fpPrbLoad

    print 'Perform initialization'
    sys.stdout.flush()

    #
    #  Initialize global variables.
    #
    dataDir = os.environ['DATADIR']
    createdBy = os.environ['CREATEDBY']
    refKey = accessionlib.get_Object_key(J_NUMBER, 'Reference')
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
    #  Open the QC report file.
    #
    fpQC = reportlib.init(outputfile='tr6739probe.rpt',
                          outputdir=os.environ['RPTDIR'], sqlLogging = 0)

    #
    #  Open the temp probe file.
    #
    tmpPrbFile = dataDir + '/TMP_Probe.bcp'

    try:
        fpTmpPrb = open(tmpPrbFile,'w')
    except:
        sys.stderr.write('Could not open temp probe file: ' + tmpPrbFile)
        exit(1)

    #
    #  Open the temp probe/marker file.
    #
    tmpPrbMrkFile = dataDir + '/TMP_ProbeMarker.bcp'

    try:
        fpTmpPrbMrk = open(tmpPrbMrkFile,'w')
    except:
        sys.stderr.write('Could not open temp probe/marker file: ' + \
                         tmpPrbMrkFile)
        exit(1)

    #
    #  Open the probe load file.
    #
    prbLoadFile = dataDir + '/probe.txt'

    try:
        fpPrbLoad = open(prbLoadFile,'w')
    except:
        sys.stderr.write('Could not open probe load file: ' + prbLoadFile)
        exit(1)

    return


# Purpose: Perform cleanup steps for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def finalize():
    global fpQC, fpPrbLoad

    db.useOneConnection(0)

    #
    #  Close the QC report file.
    #
    fpQC.close()

    #
    #  Close the probe load file.
    #
    fpPrbLoad.close()

    return


# Purpose: Preprocess the retina data input file to report any records that
#          have discrepancies and load probe-related information from the
#          remaining records into the TMP_Probe table.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def preprocess():
    global dbServer, dbUser, dbPasswordFile
    global tmpPrbFile
    global fpQC, fpTmpPrb
    global rptNoSymbol, rptNoStrength, rptBadStrength, rptNoTrans
    global rptNoteAAA, rptNoNote, rptNoTissues

    print 'Preprocess the retina data input file'
    sys.stdout.flush()

    rptNoSymbol = []
    rptNoTrans = []
    rptNoStrength = []
    rptBadStrength = []
    rptNoNote = []
    rptNoteAAA = []
    rptNoTissues = []
    qcLines = []
    hasTissue = 0
    lineError = 0
    goodLines = 0

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
        #
        #  Determine if the marker symbol is missing.
        #
        if parser.getSymbol() == '':
            if lineError == 0:
                rptNoSymbol.append({'lineNum':data['lineNum']})
                lineError = 1

        #
        #  Determine if a tissue cannot be translated.
        #
        if translation == None:
            rptNoTrans.append({'lineNum':data['lineNum'],
                               'label':data['label'],
                               'tissue':data['tissue']})
            lineError = 1
        else:
            #
            #  Determine if the tissue translates to something other than
            #  "no tissue".
            #
            if translation['gxd_tissue'] != 'no tissue':
                hasTissue = 1

                #
                #  Determine if a tissue has a missing strength.
                #
                if data['strength'] == '':
                    rptNoStrength.append({'lineNum':data['lineNum'],
                                          'label':data['label'],
                                          'tissue':data['tissue']})
                    lineError = 1

                #
                #  Determine if a tissue has an invalid strength.
                #
                else:
                    try:
                        value = string.atof(data['strength'])
                    except:
                        rptBadStrength.append({'lineNum':data['lineNum'],
                                               'label':data['label'],
                                               'tissue':data['tissue'],
                                               'strength':data['strength']})
                        lineError = 1

                #
                #  Determine if a tissue has a missing required note.
                #
                if translation['required'] == 'yes' and data['note'] == '':
                    rptNoNote.append({'lineNum':data['lineNum'],
                                      'label':data['label'],
                                      'tissue':data['tissue']})
                    lineError = 1

        #
        #  Determine if a note begins with "AAA".
        #
        if data['note'][0:3] == 'AAA':
            rptNoteAAA.append({'lineNum':data['lineNum'],
                               'label':data['label'],
                               'tissue':data['tissue'],
                               'note':data['note']})
            lineError = 1

        #
        #  If this is the last tissue for the input record, handle the
        #  reporting for this record.
        #
        if last == 1:
            #
            #  Make sure at least one tissue from the input record
            #  translated to something other than "no tissue".
            #
            if hasTissue == 0:
                rptNoTissues.append({'lineNum':data['lineNum']})
                lineError = 1

            #
            #  If there was any error found for the input record, save the
            #  line number.  Otherwise, write the probe/marker information
            #  to the output file so it can be loaded into the Tmp_Probe
            #  table.
            #
            if lineError == 1:
                qcLines.append(data['lineNum'])
            else:
                fpTmpPrb.write(str(data['lineNum']) + TAB +
                               parser.getSymbol() + TAB +
                               parser.getSeqID() + TAB +
                               parser.getProbeAccID() + TAB +
                               parser.getProbeName() + TAB +
                               parser.getCloneLib() + CRT)
                goodLines = goodLines + 1

            hasTissue = 0
            lineError = 0

        #
        #  Get the next set of tissue data.
        #
        data, translation, last = parser.nextTissue()

    #
    #  Destroy the parser object.
    #
    parser = None

    #
    #  Close the temp probe file.
    #
    fpTmpPrb.close()

    #
    #  Load the TMP_Probe table with the temp probe file.
    #
    print 'Load the TMP_Probe table with records to be processed'
    sys.stdout.flush()

    bcpCmd = 'cat ' + dbPasswordFile + ' | bcp tempdb..TMP_Probe in ' + \
             tmpPrbFile + ' -c -S' + dbServer + ' -U' + dbUser
    os.system(bcpCmd)

    print 'Records to be processed: ' + str(goodLines)
    sys.stdout.flush()

    print 'Write discrepancies to the QC report'
    sys.stdout.flush()

    #
    #  Report lines where the marker symbol is missing.
    #
    fpQC.write('Input lines where the marker symbol is missing.' + 2*CRT)
    fpQC.write('Line' + CRT)
    fpQC.write('----' + CRT)

    for r in rptNoSymbol:
        fpQC.write("%-4d" % (r['lineNum']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoSymbol)) + CRT)

    #
    #  Report cases where a tissue cannot be translated.
    #
    fpQC.write(3*CRT)
    fpQC.write('Tissues that could not be translated.' + 2*CRT)
    fpQC.write('Line Label    Tissue' + CRT)
    fpQC.write('---- -------- --------------------' + CRT)

    for r in rptNoTrans:
        fpQC.write("%-4d %-8s %-20s" % (r['lineNum'], r['label'], r['tissue']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoTrans)) + CRT)

    #
    #  Report cases where a tissue has a missing strength.
    #
    fpQC.write(3*CRT)
    fpQC.write('Tissues that are missing a strength.' + 2*CRT)
    fpQC.write('Line Label    Tissue' + CRT)
    fpQC.write('---- -------- --------------------' + CRT)

    for r in rptNoStrength:
        fpQC.write("%-4d %-8s %-20s" % (r['lineNum'], r['label'], r['tissue']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoStrength)) + CRT)

    #
    #  Report cases where a tissue has an invalid strength.
    #
    fpQC.write(3*CRT)
    fpQC.write('Tissues that have an invalid strength.' + 2*CRT)
    fpQC.write('Line Label    Tissue               Strength' + CRT)
    fpQC.write('---- -------- -------------------- ---------------' + CRT)

    for r in rptBadStrength:
        fpQC.write("%-4d %-8s %-20s %-15s" % (r['lineNum'], r['label'], r['tissue'], r['strength']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptBadStrength)) + CRT)

    #
    #  Report cases where a tissue has a missing required note.
    #
    fpQC.write(3*CRT)
    fpQC.write('Tissues that are missing a required note.' + 2*CRT)
    fpQC.write('Line Label    Tissue' + CRT)
    fpQC.write('---- -------- --------------------' + CRT)

    for r in rptNoNote:
        fpQC.write("%-4d %-8s %-20s" % (r['lineNum'], r['label'], r['tissue']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoNote)) + CRT)

    #
    #  Report cases where a note begins with "AAA".
    #
    fpQC.write(3*CRT)
    fpQC.write('Notes that begin with "AAA".' + 2*CRT)
    fpQC.write('Line Label    Tissue               Note' + CRT)
    fpQC.write('---- -------- -------------------- ------------------------------------------------------------' + CRT)

    for r in rptNoteAAA:
        fpQC.write("%-4d %-8s %-20s %-s" % (r['lineNum'], r['label'],
                                            r['tissue'], r['note']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoteAAA)) + CRT)

    #
    #  Report lines where all tissues are blank or translate to "no tissue".
    #
    fpQC.write(3*CRT)
    fpQC.write('Input lines where all tissues are blank or translate to "no tissue".' + 2*CRT)
    fpQC.write('Line' + CRT)
    fpQC.write('----' + CRT)

    for r in rptNoTissues:
        fpQC.write("%-4d" % (r['lineNum']))
        fpQC.write(CRT)
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(rptNoTissues)) + CRT)

    #
    #  For each input line number that had a QC error, insert it into
    #  the TMP_QC table.
    #
    print 'Load the TMP_QC table with line numbers that have discrepancies'
    sys.stdout.flush()

    for lineNum in qcLines:
        db.sql('insert into tempdb..TMP_QC values (' + \
               str(lineNum) + ',"probe")', 'auto')

    print 'Records with discrepancies: ' + str(len(qcLines))
    sys.stdout.flush()

    return


# Purpose: Write discrepancy conditions to the QC report.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def reportDiscrepancies():
    global fpQC

    print 'Identify input records with discrepancies in MGD'
    sys.stdout.flush()

    discrepCount = 0
    qcLines = []
    cmd = []

    #
    #  Find probes that are associated with a different marker than the
    #  one in the input file.
    #
    cmd.append('select tp.lineNum, tp.symbol, tp.seqID, tp.probeAccID, ' + \
                      'm1.symbol "db_symbol", pm1.relationship ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'ACC_Accession a, ' + \
                    'PRB_Marker pm1, ' + \
                    'MRK_Marker m1 ' + \
               'where tp.symbol is not null and ' + \
                     'tp.probeAccID = a.accID and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a._MGIType_key = 3 and ' + \
                     'a.preferred = 1 and ' + \
                     'a.prefixPart= "MGI:" and ' + \
                     'a._Object_key = pm1._Probe_key and ' + \
                     'pm1._Marker_key = m1._Marker_key and ' + \
                     'm1._Organism_key = 1 and ' + \
                     'not exists (select 1 ' + \
                                 'from PRB_Marker pm2, ' + \
                                      'MRK_Marker m2 ' + \
                                 'where a._Object_key = pm2._Probe_key and ' + \
                                       'pm2._Marker_key = m2._Marker_key and ' + \
                                       'm2._Organism_key = 1 and ' + \
                                       'm2.symbol = tp.symbol) and ' + \
                     'tp.lineNum not in (select lineNum from tempdb..TMP_QC) ' + \
               'order by tp.symbol')

    #
    #  Find probes that have a marker in the input file, but that marker
    #  does not exist in the database, so an association cannot be made.
    #
    cmd.append('select tp.lineNum, tp.symbol, tp.seqID, tp.probeAccID ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'ACC_Accession a ' + \
               'where tp.symbol is not null and ' + \
                     'tp.probeAccID = a.accID and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a._MGIType_key = 3 and ' + \
                     'a.preferred = 1 and ' + \
                     'a.prefixPart= "MGI:" and ' + \
                     'not exists (select 1 ' + \
                                 'from PRB_Marker pm ' + \
                                 'where a._Object_key = pm._Probe_key) and ' + \
                     'not exists (select 1 ' + \
                                 'from MRK_Marker m ' + \
                                 'where tp.symbol = m.symbol) and ' + \
                     'tp.lineNum not in (select lineNum from tempdb..TMP_QC) ' + \
               'order by tp.symbol')

    #
    #  Find sequence IDs in the input file that need to be associated
    #  with the given probe.
    #
    cmd.append('select tp.lineNum, tp.symbol, tp.seqID, tp.probeAccID ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'ACC_Accession a1 ' + \
               'where tp.symbol is not null and ' + \
                     'tp.seqID is not null and ' + \
                     'tp.probeAccID = a1.accID and ' + \
                     'a1._LogicalDB_key = 1 and ' + \
                     'a1._MGIType_key = 3 and ' + \
                     'a1.preferred = 1 and ' + \
                     'a1.prefixPart= "MGI:" and ' + \
                     'not exists (select 1 ' + \
                                 'from ACC_Accession a2 ' + \
                                 'where a2._LogicalDB_key = 9 and ' + \
                                       'a2._MGIType_key = 3 and ' + \
                                       'a2._Object_key = a1._Object_key and ' + \
                                       'a2.accID = tp.seqID) and ' + \
                     'tp.lineNum not in (select lineNum from tempdb..TMP_QC) ' + \
               'order by tp.symbol')

    #
    #  Find input records where the probe cannot be created because the
    #  probe name is missing.
    #
    cmd.append('select lineNum, symbol, seqID, cloneLibrary ' + \
               'from tempdb..TMP_Probe ' + \
               'where probeAccID is null and ' + \
                     'symbol is not null and ' + \
                     'probeName is null and ' + \
                     'lineNum not in (select lineNum from tempdb..TMP_QC) ' + \
               'order by symbol')

    #
    #  Find symbols that have no MGI ID.
    #
    cmd.append('select tp.lineNum, tp.symbol, tp.seqID, tp.probeAccID, ' + \
                      'tp.probeName, tp.cloneLibrary ' + \
               'from tempdb..TMP_Probe tp ' + \
               'where tp.symbol is not null and ' + \
                     'not exists (select 1 ' + \
                                 'from MRK_Marker m3, ' + \
                                      'ACC_Accession a ' + \
                                 'where tp.symbol = m3.symbol and ' + \
                                       'm3._Organism_key = 1 and ' + \
                                       'm3._Marker_key = a._Object_key and ' + \
                                       'a._MGIType_key = 2 and ' + \
                                       'a._LogicalDB_key = 1 and ' + \
                                       'a.prefixPart = "MGI:" and ' + \
                                       'a.preferred = 1) and ' + \
                     'tp.lineNum not in (select lineNum from tempdb..TMP_QC) ' + \
               'order by tp.symbol')

    results = db.sql(cmd, 'auto')

    print 'Write discrepancies to the QC report'
    sys.stdout.flush()

    fpQC.write(3*CRT)
    fpQC.write('Probes that are associated with a different marker than ' +
               'what is in the input file.' + 2*CRT)
    fpQC.write('Line Symbol         Seq ID     Probe ID     DB Symbol      Rel' + CRT)
    fpQC.write('---- -------------- ---------- ------------ -------------- ---' + CRT)

    for r in results[0]:
        if r['seqID'] == None:
            seqID = ""
        else:
            seqID = r['seqID']
        fpQC.write("%-4d %-14s %-10s %-12s %-14s %-3s" %
                   (r['lineNum'], r['symbol'], seqID, r['probeAccID'],
                   r['db_symbol'], r['relationship']))
        fpQC.write(CRT)
        qcLines.append(r['lineNum'])
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(results[0])) + CRT)
    discrepCount = discrepCount + len(results[0])

    fpQC.write(3*CRT)
    fpQC.write('Probes that have a marker in the input file, but that ' +
               'marker does not exist in the database.' + 2*CRT)
    fpQC.write('Line Symbol         Seq ID     Probe ID' + CRT)
    fpQC.write('---- -------------- ---------- ------------' + CRT)

    for r in results[1]:
        if r['seqID'] == None:
            seqID = ""
        else:
            seqID = r['seqID']
        fpQC.write("%-4d %-14s %-10s %-12s" %
                   (r['lineNum'], r['symbol'], seqID, r['probeAccID']))
        fpQC.write(CRT)
        qcLines.append(r['lineNum'])
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(results[1])) + CRT)
    discrepCount = discrepCount + len(results[1])

    fpQC.write(3*CRT)
    fpQC.write('Sequence IDs in the input file that need to be associated ' +
               'with the given probe.' + 2*CRT)
    fpQC.write('Line Symbol         Seq ID     Probe ID' + CRT)
    fpQC.write('---- -------------- ---------- ------------' + CRT)

    for r in results[2]:
        if r['seqID'] == None:
            seqID = ""
        else:
            seqID = r['seqID']
        fpQC.write("%-4d %-14s %-10s %-12s" %
                   (r['lineNum'], r['symbol'], seqID, r['probeAccID']))
        fpQC.write(CRT)
        qcLines.append(r['lineNum'])
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(results[2])) + CRT)
    discrepCount = discrepCount + len(results[2])

    fpQC.write(3*CRT)
    fpQC.write('Probes that cannot be created because the probe name ' +
               'is missing.' + 2*CRT)
    fpQC.write('Line Symbol         Seq ID     Clone Library' + CRT)
    fpQC.write('---- -------------- ---------- ----------------' + CRT)

    for r in results[3]:
        if r['seqID'] == None:
            seqID = ""
        else:
            seqID = r['seqID']
        if r['cloneLibrary'] == None:
            cloneLibrary = ""
        else:
            cloneLibrary = r['cloneLibrary']
        fpQC.write("%-4d %-14s %-10s %-16s" %
                   (r['lineNum'], r['symbol'], seqID, cloneLibrary))
        fpQC.write(CRT)
        qcLines.append(r['lineNum'])
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(results[3])) + CRT)
    discrepCount = discrepCount + len(results[3])

    fpQC.write(3*CRT)
    fpQC.write('Symbols that have no MGI ID.' + 2*CRT)
    fpQC.write('Line Symbol         Seq ID     Probe ID     Probe Name                   Clone Library' + CRT)
    fpQC.write('---- -------------- ---------- ------------ ---------------------------- ----------------' + CRT)

    for r in results[4]:
        if r['seqID'] == None:
            seqID = ""
        else:
            seqID = r['seqID']
        if r['probeAccID'] == None:
            probeAccID = ""
        else:
            probeAccID = r['probeAccID']
        if r['cloneLibrary'] == None:
            cloneLibrary = ""
        else:
            cloneLibrary = r['cloneLibrary']
        fpQC.write("%-4d %-14s %-10s %-12s %-28s %-16s" %
                   (r['lineNum'], r['symbol'], seqID, probeAccID,
                    r['probeName'], cloneLibrary))
        fpQC.write(CRT)
        qcLines.append(r['lineNum'])
    fpQC.write(CRT + 'Number of Discrepancies: ' + str(len(results[4])) + CRT)
    discrepCount = discrepCount + len(results[4])

    print 'Number of discrepancies: ' + str(discrepCount)
    sys.stdout.flush()

    #
    #  For each input line number that had a QC error, insert it into
    #  the TMP_QC table.
    #
    print 'Load the TMP_QC table with line numbers that have discrepancies'
    sys.stdout.flush()

    for lineNum in qcLines:
        db.sql('insert into tempdb..TMP_QC values (' + \
               str(lineNum) + ',"probe")', 'auto')

    print 'Records with discrepancies: ' + str(len(qcLines))
    sys.stdout.flush()

    #
    #  Delete any records from the TMP_Probe table that have been found
    #  to have discrepancies.
    #
    print 'Delete records from the TMP_Probe table that have discrepancies'
    sys.stdout.flush()

    db.sql('delete from tempdb..TMP_Probe ' + \
           'where lineNum in (select lineNum from tempdb..TMP_QC)', 'auto')

    return


# Purpose: If an input record has an MGI ID, it represents an existing
#          probe. Determine if the probe's relationship to the given
#          marker needs to be updated or added.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def updateExistingProbes():

    print 'Update probe/marker relationships for existing probes'
    sys.stdout.flush()

    cmd = []

    #
    #  Find probes from the input file that exist in MGD and are already
    #  associated with the given marker. If the probe/marker
    #  relationship is "P", update it to be "H".
    #
    cmd.append('update PRB_Marker ' + \
               'set _Refs_key = ' + str(refKey) + ', ' + \
                   'relationship = "H", ' + \
                   '_ModifiedBy_key = u._User_key, ' + \
                   'modification_date = getdate() ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'ACC_Accession a, ' + \
                    'PRB_Marker pm, ' + \
                    'MRK_Marker m, ' + \
                    'MGI_User u ' + \
               'where tp.symbol is not null and ' + \
                     'tp.probeAccID = a.accID and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a._MGIType_key = 3 and ' + \
                     'a.preferred = 1 and ' + \
                     'a.prefixPart= "MGI:" and ' + \
                     'a._Object_key = pm._Probe_key and ' + \
                     'pm.relationship = "P" and ' + \
                     'pm._Marker_key = m._Marker_key and ' + \
                     'm._Organism_key = 1 and ' + \
                     'm.symbol = tp.symbol and ' + \
                     'u.login = "' + createdBy + '"')

    #
    #  Find probes from the input file that exist in MGD, but are not
    #  associated with the given marker.
    #
    cmd.append('select a._Object_key, m._Marker_key, ' + \
                       str(refKey) + '"_Refs_key", ' + \
                      '"H" "relationship", ' + \
                      'u._User_key "_CreatedBy_key", ' + \
                      'u._User_key "_ModifiedBy_key", ' + \
                      'getdate() "creation_date", ' + \
                      'getdate() "modification_date" ' + \
               'into #NewProbeMarker ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'ACC_Accession a, ' + \
                    'MRK_Marker m, ' + \
                    'MGI_User u ' + \
               'where tp.symbol is not null and ' + \
                     'tp.probeAccID = a.accID and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a._MGIType_key = 3 and ' + \
                     'a.preferred = 1 and ' + \
                     'a.prefixPart= "MGI:" and ' + \
                     'not exists (select 1 ' + \
                                 'from PRB_Marker pm ' + \
                                 'where a._Object_key = pm._Probe_key) and ' + \
                     'tp.symbol = m.symbol and ' + \
                     'm._Organism_key = 1 and ' + \
                     'u.login = "' + createdBy + '"')

    cmd.append('select count(*) "count" from #NewProbeMarker')

    results = db.sql(cmd, 'auto')

    #
    #  If there are any cases where the probe/marker association is
    #  missing, add them.
    #
    print 'Add probe/marker relationships for existing probes'
    sys.stdout.flush()

    if results[2][0]['count'] > 0:
        cmd = []
        cmd.append('insert into PRB_Marker ' + \
                   'select _Object_key, _Marker_key, _Refs_key, ' + \
                          'relationship, _CreatedBy_key, _ModifiedBy_key, ' + \
                          'creation_date, modification_date ' + \
                   'from #NewProbeMarker')

        db.sql(cmd, 'auto')

    return


# Purpose: Create an input file for the probeload product that contains
#          attributes needed to create new probes.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def prepareNewProbes():
    global fpTmpPrbMrk, fpPrbLoad

    print 'Create the "probe.txt" input file for the probeload'
    sys.stdout.flush()

    cmd = []

    #  
    #  Get the vector type for each distinct library in the input file.
    #  
    cmd.append('select distinct s.name, t.term ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'PRB_Source s, ' + \
                    'VOC_Term t, ' + \
                    'VOC_Vocab v ' + \
               'where tp.cloneLibrary = s.name and ' + \
                     's._Vector_key = t._Term_key and ' + \
                     't._Vocab_key = v._Vocab_key and ' + \
                     'v.name = "Segment Vector Type"')

    #
    #  For probes that need to be created with a named library, retrieve
    #  the data necessary to write to the probe file.
    #
    cmd.append('select tp.lineNum, tp.probeName, tp.cloneLibrary, ' + \
                      'a.accID, tp.seqID ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'MRK_Marker m, ' + \
                    'ACC_Accession a ' + \
               'where tp.probeAccID is null and ' + \
                     'tp.symbol is not null and ' + \
                     'tp.probeName is not null and ' + \
                     'tp.cloneLibrary is not null and ' + \
                     'tp.symbol = m.symbol and ' + \
                     'm._Organism_key = 1 and ' + \
                     'm._Marker_key = a._Object_key and ' + \
                     'a._MGIType_key = 2 and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a.prefixPart = "MGI:" and ' + \
                     'a.preferred = 1 ' + \
               'order by tp.symbol')

    #  
    #  Build a temp table with all of the probe names that currently
    #  exist in the database that follow the "symbol cDNA#" pattern
    #  for sequentially naming probes.
    #  
    cmd.append('select lower(name) "name" ' + \
               'into #Probes ' + \
               'from PRB_Probe ' + \
               'where name like "% cdna[0-9]%"')

    #  
    #  Find all of the probe names in the temp table that begin with
    #  any of the symbol name that need to have new probe records created.
    #  
    cmd.append('select tp.symbol, p.name ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    '#Probes p ' + \
               'where tp.probeAccID is null and ' + \
                     'tp.cloneLibrary is null and ' + \
                     'tp.probeName is not null and ' + \
                     'charindex(lower(tp.symbol) + " cdna",p.name) > 0 ' + \
               'order by tp.symbol, p.name')

    #
    #  For probes that need to be created with an anonymous library,
    #  retrieve the data necessary to write to the probe file.
    #
    cmd.append('select tp.lineNum, tp.symbol, a.accID, tp.seqID ' + \
               'from tempdb..TMP_Probe tp, ' + \
                    'MRK_Marker m, ' + \
                    'ACC_Accession a ' + \
               'where tp.probeAccID is null and ' + \
                     'tp.symbol is not null and ' + \
                     'tp.probeName is not null and ' + \
                     'tp.cloneLibrary is null and ' + \
                     'tp.symbol = m.symbol and ' + \
                     'm._Organism_key = 1 and ' + \
                     'm._Marker_key = a._Object_key and ' + \
                     'a._MGIType_key = 2 and ' + \
                     'a._LogicalDB_key = 1 and ' + \
                     'a.prefixPart = "MGI:" and ' + \
                     'a.preferred = 1 ' + \
               'order by tp.symbol')

    results = db.sql(cmd, 'auto')

    #
    #  Build a dictionary of library names and their vector types.
    #
    libVector = {}
    for r in results[0]:
        libVector[r['name']] = r['term']

    #
    #  Process each record for a name library.
    #
    for r in results[1]:
        lineNum = r['lineNum']
        probeName = r['probeName']
        library = r['cloneLibrary']
        mgiID = r['accID']
        seqID = r['seqID']
        vectorType = libVector[library]

        if seqID == None:
            seqID = ''
        else:
            seqID = LOGICAL_DB + ":" + seqID + "|"

        fpPrbLoad.write(probeName + TAB + J_NUMBER + TAB + library + TAB +
                        NULL + TAB + NULL + TAB + NULL + TAB + NULL + TAB +
                        NULL + TAB + NULL + TAB + vectorType + TAB +
                        SEGMENT_TYPE + TAB + NULL + TAB + INSERT_SITE + TAB +
                        NULL + TAB + mgiID + TAB + RELATIONSHIP + TAB +
                        seqID + TAB + NULL + TAB + createdBy + CRT)

        fpTmpPrbMrk.write(str(lineNum) + TAB + mgiID + TAB + TAB +
                          probeName + CRT)

    #
    #  Build a dictionary that keeps track of the largest sequential
    #  number used in conjunction with a symbol for creating the name
    #  of a probe that follows the "symbol cDNA#" pattern.  For example,
    #  if a symbol is "Wnt1" and a probe is found with the name
    #  "Wnt1 cDNA6", the dictionary entry for "Wtn1" would be 6.
    #
    symbolNumber = {}
    for r in results[3]:
        symbol = r['symbol']
        name = r['name']

        #
        #  Extract the sequential number from the probe name.
        #
        s = string.splitfields(name,string.lower(symbol)+' cdna')
        number = string.atoi(s[1])

        #
        #  The the symbol already exists and the new number is less than
        #  the one in the dictionary, skip it.
        #
        if symbolNumber.has_key(symbol):
            if number < symbolNumber[symbol]:
                continue
        #
        #  Store the new number for the symbol.
        #
        symbolNumber[symbol] = number

    #
    #  Process each record for an anonymous library.
    #
    for r in results[4]:
        lineNum = r['lineNum']
        symbol = r['symbol']
        mgiID = r['accID']
        seqID = r['seqID']

        if seqID == None:
            seqID = ''
        else:
            seqID = LOGICAL_DB + ":" + seqID + "|"

        if symbolNumber.has_key(symbol):
            number = symbolNumber[symbol] + 1
        else:
            number = 1
        symbolNumber[symbol] = number

        probeName = symbol + ' cDNA' + str(number)

        fpPrbLoad.write(probeName + TAB + J_NUMBER + TAB + NULL + TAB +
                        ORGANISM + TAB + STRAIN + TAB + TISSUE + TAB +
                        GENDER + TAB + CELLLINE + TAB + AGE + TAB +
                        VECTOR_TYPE + TAB + SEGMENT_TYPE + TAB +
                        NULL + TAB + NULL + TAB + NULL + TAB +
                        mgiID + TAB + RELATIONSHIP + TAB +
                        seqID + TAB + NULL + TAB + createdBy + CRT)

        fpTmpPrbMrk.write(str(lineNum) + TAB + mgiID + TAB + TAB +
                          probeName + CRT)

    return


# Purpose: Load the TMP_ProbeMarker table.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def loadTmpPrbMrkTable():
    global dbServer, dbUser, dbPasswordFile
    global tmpPrbMrkFile
    global fpTmpPrbMrk

    print 'Load the TMP_ProbeMarker table'
    sys.stdout.flush()

    #
    #  Close the temp probe/marker file.
    #
    fpTmpPrbMrk.close()

    #
    #  Load the TMP_ProbeMarker table with the temp probe/marker file.
    #  This will load records that represent probe/marker relationships
    #  for new probes that will be created by the probeload.
    #
    bcpCmd = 'cat ' + dbPasswordFile + ' | bcp tempdb..TMP_ProbeMarker in ' + \
             tmpPrbMrkFile + ' -c -S' + dbServer + ' -U' + dbUser
    os.system(bcpCmd)

    #
    #  Load the TMP_ProbeMarker table with additional records that
    #  represent probe/marker relationships for probes that already
    #  exist and have a probe accession ID in the input file.
    #
    db.sql('insert into tempdb..TMP_ProbeMarker ' + \
           'select tp.lineNum, a.accID, tp.probeAccID, null ' + \
           'from tempdb..TMP_Probe tp, ' + \
                'MRK_Marker m, ' + \
                'ACC_Accession a ' + \
           'where tp.probeAccID is not null and ' + \
                 'tp.symbol = m.symbol and ' + \
                 'm._Organism_key = 1 and ' + \
                 'm._Marker_key = a._Object_key and ' + \
                 'a._MGIType_key = 2 and ' + \
                 'a._LogicalDB_key = 1 and ' + \
                 'a.prefixPart = "MGI:" and ' + \
                 'a.preferred = 1',
           'auto')

    return


#
#  MAIN
#
initialize()
preprocess()
reportDiscrepancies()
updateExistingProbes()
prepareNewProbes()
loadTmpPrbMrkTable()
finalize()

