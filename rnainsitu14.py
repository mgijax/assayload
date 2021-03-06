
#
# Program: rnainsitu14.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To translate RNA InSitu 14.5 input files into input files
#	for the insituload.py program.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
#	rnainsitu14.py
#
# Envvars:
#
# Inputs:
#
#       In_Situ.txt, a tab-delimited file in the format:
#		field 1: Human Gene
#		field 2: Mouse Gene Symbol
#               field 3: MGI Marker Accession ID
#		field 4: ISH_number
#		field 5: Specimen
#		field 6: Tissue Quality
#		field 7: Overall Expression
#		field 8-50: Results
#		field 51: Image 1
#		field 52: Image 2
#
#	In_Situ_Tissues.txt, a tab-delimited file in the format:
#		field 1: Reported Tissue
#		field 2: MGI Tissue
#		field 3: Theiler Stage
#
#	probe_table.txt, a tab-delimited file in the format:
#		field 1: Human Gene
#		field 2: Mouse Gene
#               field 3: MGI Marker Accession ID
#		field 4: Clone Name
#		field 5: Clone Origin
#		field 6: Clone Acc ID
#		field 7: Clone Mapping
#		field 8: Clone MGI ID
#		field 9: Clone Library ID
#		field 10: Clone Library Name
#
# Outputs:
#
#       4 tab-delimited files:
#
#	In_Situ_probeprep.txt
#	In_Situ_assay.txt
#	In_Situ_specimen.txt
#	In_Situ_results.txt
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

inInSituFile = ''	# file descriptor
inTissueFile = ''	# file descriptor
inProbeFile = ''	# file descriptor
prepFile = ''		# file descriptor
assayFile = ''          # file descriptor
specimenFile = ''       # file descriptor
resultsFile = ''        # file descriptor

datadir = os.environ['INSITU14DATADIR']

inInSituFileName = datadir + '/tr4800/14.5_In_situ.txt'
inTissueFileName = datadir + '/tr4800/14.5_In_Situ_tissues.txt'
inProbeFileName = datadir + '/tr4800/probe_table.txt'
prepFileName = datadir + '/In_Situ_probeprep.txt'
assayFileName = datadir + '/In_Situ_assay.txt'
specimenFileName = datadir + '/In_Situ_specimen.txt'
resultsFileName = datadir + '/In_Situ_results.txt'

# constants for probe prep
prepType = 'RNA'
hybridization = 'Antisense'
labelledWith = 'Digoxigenin'
labelCoverage = 'continuously'
visualizedWith = 'Alkaline phosphatase'

# constants for assay
reference = 'J:80502'
assayType = 'RNA In Situ'
createdBy = os.environ['CREATEDBY']

# constants for specimen
specimenLabel = '%s 14.5dpc'
genotype = 'MGI:2166653'
age = 'embryonic day 14.5'
ageNote = 'Age of embryo at noon of plug day not specified in reference.'
sex = 'Not Specified'
fixation = '4% Paraformaldehyde'
embedding = 'Cryosection'
specimenHybridization = 'section'
specimenNote = NULL
resultNote = NULL
epiTissue = 'epithalamus'
epiNote = 'Expression was restricted to the pineal gland primordium.'
ubiExpression = 'ubiquitous'
ubiPattern = 'U'

# translation of patterns

patternTrans = {'U':'Homogeneous', 'R':'Regionally restricted', 'I':'Regionally restricted', \
   'U,R':'Regionally restricted', 'R,U':'Regionally restricted', 'i':'Regionally restricted'}

ABSENT = 'Absent'
NA = 'Not Applicable'

# translation of input file strengths and MGI strengths

strengthTrans = {'*':'Weak', '**':'Moderate', '***':'Strong', '':'Absent', \
    'U,R':'Present', 'R,U':'Present', '-':'Absent', '+':'Weak'}
presentStrength = ['U,R', 'R,U']

# Purpose: prints error message and exits
# Returns: nothing
# Assumes: nothing
# Effects: exits with exit status
# Throws: nothing

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (str.
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
    global inInSituFile, inTissueFile, inProbeFile, prepFile, assayFile, specimenFile, resultsFile
 
    try:
        inInSituFile = open(inInSituFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inInSituFileName)

    try:
        inTissueFile = open(inTissueFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inTissueFileName)

    try:
        inProbeFile = open(inProbeFileName, 'r')
    except:
        exit(1, 'Could not open file %s\n' % inProbeFileName)

    try:
        prepFile = open(prepFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % prepFileName)

    try:
        assayFile = open(assayFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % assayFileName)

    try:
        specimenFile = open(specimenFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % specimenFileName)

    try:
        resultsFile = open(resultsFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % resultsFileName)

    return

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  writes data to output files
# Throws:   nothing

def process():

    global resultNote

    tissueTrans = {}	# maps input tissue to MGI tissue and Theiler Stage
    probeTrans = {}	# maps probe to MGI Gene

    for line in inTissueFile.readlines():
        tokens = str.split(line[:-1], TAB)
        badTissue = tokens[0]
        goodTissue = tokens[1]
        theilerStage = tokens[2]

        key = badTissue
        value = goodTissue + '|' + theilerStage
        if key not in tissueTrans:
            tissueTrans[key] = []
        tissueTrans[key].append(value)

    for line in inProbeFile.readlines():
        tokens = str.split(line[:-1], TAB)
        mgiID = tokens[2]
        probeID = tokens[7]

        if len(mgiID) == 0:
            continue

        key = mgiID
        value = probeID
        if key not in probeTrans:
            probeTrans[key] = []
        probeTrans[key].append(value)

    assay = 0	# unique Assay ID

    # For each line in the input file

    for line in inInSituFile.readlines():

        # Split the line into tokens
        tokens = str.split(line[:-1], TAB)

        # processing first line (header)
        # grab the Tissue headings

        if assay == 0:
            tissueLabels = tokens[7:50]
            assay = assay + 1
            continue

        # else process an actual data line

        try:
            humanGene = tokens[0]
            mouseGene = tokens[1]
            accID = str.strip(tokens[2])
            ishNumber = tokens[3]
            specimen = tokens[4]
            tissueQuality = tokens[5]
            overallExpression = tokens[6]
            results = tokens[7:50]
            imageFileName1 = tokens[50]
            imageFileName2= tokens[51]

        except:
            print('Invalid Line (%d): %s\n' % (assay, line))

        if len(mouseGene) == 0:
            continue

        # create one assay per probe for given marker

        if accID not in probeTrans:
            print('Cannot find MGI ID in Probe file: %s\n' % (accID))
            continue

        for probeID in probeTrans[accID]:

            # write the probe prep information

            prepFile.write(str(assay) + TAB + \
                probeID + TAB + \
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
                TAB + \
                createdBy + CRT)

            # write the specimen (one for each Assay)

            specimen = 1

            specimenFile.write(str(assay) + TAB + \
                    str(specimen) + TAB + \
                    specimenLabel % (mouseGene) + TAB + \
                    genotype + TAB + \
                    age + TAB + \
                    ageNote + TAB + \
                    sex + TAB + \
                    fixation + TAB + \
                    embedding + TAB + \
                    specimenHybridization + TAB + \
                    specimenNote + CRT)

            # set defaults
            defaultStrength = ABSENT
            defaultPattern = NA

            # if overall expression is ubiquitous, then set the default Strength to
            # the strength value specified in this column

            if len(overallExpression) > 0:
                try:
                    [oExpression, oStrength] = str.split(overallExpression, ' ')

                    if oExpression == ubiExpression:
                        defaultStrength = strengthTrans[oStrength]
                        defaultPattern = patternTrans[ubiPattern]

                except:
                    pass

            # one result for each Tissue 

            result = 1
            for i in range(len(tissueLabels)):

                # Translate the Tissue into a Tissue and Age

                for t in tissueTrans[tissueLabels[i]]:

                    [tissue, theilerStage] = str.split(t, '|')
                    resultNote = NULL

                    if len(results) < i + 1:
                        # not every tissue was assayed; no results
                        continue

                    # if there are results...

                    if len(results[i]) > 1:

                        [inStrength, inPattern] = str.split(results[i], ' ')

                        if inPattern in presentStrength:
                            strength = strengthTrans[inPattern]
                            pattern = patternTrans[inPattern]
                        else:
                            strength = strengthTrans[inStrength]
                            pattern = patternTrans[inPattern]

                        if tissue == epiTissue:
                           resultNote = epiNote
                           pattern = patternTrans['R']

                    # no strength or pattern given

                    else:
                        strength = defaultStrength
                        pattern = defaultPattern

                    resultsFile.write(str(assay) + TAB + \
                        str(specimen) + TAB + \
                        str(result) + TAB + \
                        strength + TAB + \
                        pattern + TAB + \
                        tissue + TAB + \
                        theilerStage + TAB + \
                        resultNote + CRT)

                result = result + 1

            specimen = specimen + 1
            assay = assay + 1

        # end of for probeID in probeTrans[accID]
    # end of "for line in inInSituFile.readlines():"

#
# Main
#

init()
process()
exit(0)
