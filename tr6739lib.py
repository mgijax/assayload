#!/usr/local/bin/python
#
#  tr6739lib.py
###########################################################################
#
#  Purpose:
#
#      This module contains classes used by the assay load for tr6739.
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

#
#  CONSTANTS
#

TAB = '\t'

# List of specimen labels, one for each of the 12 sets of columns that
# represent an age.
#
LABEL_LIST = [ 'E12','E14','E16','E18','P0/1','P2/3',
               'P4/5','P6/7','P8/10','Adult','CNS','Embryo' ]

# List of age specifications, one for each of the 12 sets of columns that
# represent an age.
#
AGE_LIST = [ 'embryonic day 12', 'embryonic day 14',
             'embryonic day 16', 'embryonic day 18',
             'postnatal day 0-1', 'postnatal day 2-3',
             'postnatal day 4-5', 'postnatal day 6-7',
             'postnatal day 8-10', 'postnatal adult',
             'postnatal day 6', 'embryonic day 16' ]

# List of Theiler stages, one for each of the 12 sets of columns that
# represent an age.
#
THEILER_STAGE_LIST = [ '20', '22', '24', '26', '28', '28',
                       '28', '28', '28', '28', '28', '24' ]

# Column numbers within the Retina_data.txt input file.
#
GENE_DESC_COLUMN = 2
SYMBOL_COLUMN = 3
SEQ_ID_COLUMN = 4
PROBE_ACC_ID_COLUMN = 5
PROBE_NAME_COLUMN = 6
CLONE_LIB_COLUMN = 7

# A list of values that are used when parsing the Retina_data.txt
# input file. Each entry in the list corresponds to a set of
# tissue, note and strength columns from an input record. Each entry
# contains the following values for a specific age/tissue:
# 1) Column # containing the tissue
# 2) Column # containing the image(s) (if any)
# 3) A number (1-12) that represents one of the groups of columns in the
#    Retina_data.txt input file for a given age.
# 4) A number (1-6) that represents one of the groups of columns in the
#    Tissue_Translation.txt input file where the given tissue can be
#    translated.
#
PARSER_VALUES = [ [10,0,1,1],
                  [13,0,1,1],
                  [19,18,2,2],
                  [22,18,2,2],
                  [28,27,3,3],
                  [31,27,3,3],
                  [37,36,4,4],
                  [40,36,4,4],
                  [43,36,4,4],
                  [49,48,5,5],
                  [52,48,5,5],
                  [55,48,5,5],
                  [61,60,6,5],
                  [64,60,6,5],
                  [67,60,6,5],
                  [73,72,7,5],
                  [76,72,7,5],
                  [79,72,7,5],
                  [85,84,8,6],
                  [88,84,8,6],
                  [91,84,8,6],
                  [94,84,8,6],
                  [100,99,9,6],
                  [103,99,9,6],
                  [106,99,9,6],
                  [109,99,9,6],
                  [115,114,10,6],
                  [118,114,10,6],
                  [121,114,10,6],
                  [124,114,10,6],
                  [130,0,11,6],
                  [133,0,11,6],
                  [136,0,11,6],
                  [139,0,11,6],
                  [142,0,11,6],
                  [148,0,12,3],
                  [151,0,12,3],
                  [154,0,12,3],
                  [157,0,12,3],
                  [160,0,12,3] ]

# Number of fields in the Retina_data.txt input file.
MAX_FIELDS = 162

# Number of tissues that must be parsed from the Retina_data.txt input file.
MAX_TISSUES = len(PARSER_VALUES)


#
#  CLASSES
#

# IS: An object that knows how to look up an author tissue in the
#     tissue translation file to find the GXD tissue name along with
#     the required note indicator and pattern.
# HAS: Dictionaries that are loaded from the input file for each
#      of the 6 age groups.
# DOES: Performs lookups
#
class Translator:

    # Purpose: Get the tissue translation information from the input file and
    #          load it into a list of dictionaries for lookups.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def __init__ (self):
        self.translation = {},{},{},{},{},{}
        lineNum = 0

        #
        #  Set the name of the tissue translation file.
        #
        tissueTransFile = os.environ['DATADIR'] + '/Tissue_Translation.txt'

        #
        #  Open the tissue translation file.
        #
        try:
            fpTrans = open(tissueTransFile,'r')
        except:
            sys.stderr.write('Could not open input file: ' + tissueTransFile)
            exit(1)

        #
        #  Read the tissue translation file and create a list of translation
        #  values for each age.
        #
        for line in fpTrans.readlines():
            lineNum = lineNum + 1
            if lineNum == 1:
                continue

            tokens = string.split(line[:-1],TAB)

            i = 0
            while i < len(tokens):
                if tokens[i] != '':
                    data = {}
                    data['gxd_tissue'] = tokens[i+1]
                    data['required'] = tokens[i+2]
                    data['pattern'] = tokens[i+3]

                    self.translation[i / 4][string.strip(string.lower(tokens[i]))] = data

                i = i + 4

        #
        #  Close the tissue translation file.
        #
        fpTrans.close()


    # Purpose: Get the translation values for a given age and author
    #          tissue.
    # Returns: A dictionary containing the values.
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getValues (self, age, tissue):
        if tissue == '':
            tissue = 'blank'
        else:
            tissue = string.strip(string.lower(tissue))
        try:
            d = self.translation[age][tissue]
            return d
        except:
            return None


    # Purpose: Get the GXD tissue for a given age and author tissue.
    # Returns: The tissue or None (if not found)
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getGXDTissue (self, age, tissue):
        if tissue == '':
            tissue = 'blank'
        else:
            tissue = string.strip(string.lower(tissue))
        try:
            value = self.translation[age][tissue]['gxd_tissue']
            return value
        except:
            return None


    # Purpose: Determine is a note field is required for a given age
    #          and author tissue.
    # Returns: 1 (required) or 0 (not required)
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def isRequired (self, age, tissue):
        if tissue == '':
            tissue = 'blank'
        else:
            tissue = string.strip(string.lower(tissue))
        try:
            value = self.translation[age][tissue]['required']
            if string.lower(value) == "yes":
                return 1
            else:
                return 0
        except:
            return 0


    # Purpose: Get the pattern for a given age and author tissue.
    # Returns: The pattern or None (if not found)
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getPattern (self, age, tissue):
        if tissue == '':
            tissue = 'blank'
        else:
            tissue = string.strip(string.lower(tissue))
        try:
            value = self.translation[age][tissue]['pattern']
            return value
        except:
            return None


    # Purpose: Print the contents of the translator for debugging.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def debug (self):
        for trans in self.translation:
            for tissue in trans.keys():
                print "|%s|%s|%s|%s|" % (tissue,
                    trans[tissue]['gxd_tissue'],
                    trans[tissue]['required'],
                    trans[tissue]['pattern'])
            print
            sys.stdout.flush()

        return


# IS: An object that knows how to parse the retina data input file
#     and return each set of tissue data to the caller.
# HAS: A Translator object to translate each tissue that is read from
#      the input file.
# DOES: Returns each set of tissue data.
#
class DataParser:

    # Purpose: Open the retina data input file to prepare it for parsing.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def __init__ (self):
        self.fields = []
        self.tissueNum = 0
        self.foundTissue = 0
        self.geneDesc = ''
        self.symbol = ''
        self.seqID = ''
        self.probeAccID = ''
        self.probeName = ''
        self.cloneLib = ''
        self.translator = Translator()

        #
        #  Set the name of the retina data file.
        #
        retinaDataFile = os.environ['DATADIR'] + '/Retina_data.txt'

        #
        #  Open the retina data file.
        #
        try:
            self.fpData = open(retinaDataFile,'r')
        except:
            sys.stderr.write('Could not open input file: ' + retinaDataFile)
            exit(1)

        #
        #  Read one line to skip over the header.
        #
        self.fpData.readline()
        self.lineNum = 1


    # Purpose: Close the retina data input file.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def __del__ (self):
        self.translator = None
        self.fpData.close()


    # Purpose: Get the next set of tissue data from the input file.
    # Returns: A dictionary of values pertaining to the current tissue
    #          from the retina data input file, a dictionary of translated
    #          values from the tissue translation file and an indicator
    #          to tell if this is the last tissue for the current input
    #          record.
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def nextTissue (self):

        #
        #  If the first tissue has not been read yet or there are no more
        #  tissues in the current input record, a new record must be read.
        #
        if self.tissueNum == 0 or self.tissueNum == MAX_TISSUES:
            self.line = self.fpData.readline()
            if self.line == '':
                return None,None,None

            self.lineNum = self.lineNum + 1
            self.fields = string.split(self.line[:-1],TAB)

            #
            #  Expand the field list (if necessary) so that it is the
            #  proper length for parsing.
            #
            for i in range(len(self.fields),MAX_FIELDS):
                self.fields.append('')

            #
            #  Set the instance variables using the fields from the
            #  input record.
            #
            self.geneDesc = self.fields[GENE_DESC_COLUMN-1]
            self.symbol = self.fields[SYMBOL_COLUMN-1]
            self.seqID = self.fields[SEQ_ID_COLUMN-1]
            self.probeAccID = self.fields[PROBE_ACC_ID_COLUMN-1]
            self.probeName = self.fields[PROBE_NAME_COLUMN-1]
            self.cloneLib = self.fields[CLONE_LIB_COLUMN-1]

            self.tissueNum = 1
            self.foundTissue = 0

        #
        #  If there are more tissues to process from the current input
        #  record, advance the counter.
        #
        else:
            self.tissueNum = self.tissueNum + 1

        #
        #  Get the set of parse values for the current tissue.
        #
        tissueColumn = PARSER_VALUES[self.tissueNum-1][0]
        imageColumn = PARSER_VALUES[self.tissueNum-1][1]
        dataAgeGroup = PARSER_VALUES[self.tissueNum-1][2]
        transAgeGroup = PARSER_VALUES[self.tissueNum-1][3]

        #
        #  Get the tissue, note and strength from the appropriate fields.
        #
        tissue = self.fields[tissueColumn-1]
        note = self.fields[tissueColumn]
        strength = self.fields[tissueColumn+1]

        #
        #  Build a dictionary of values pertaining to the current tissue.
        #
        data = {}
        data['lineNum'] = self.lineNum
        data['tissue'] = tissue
        data['note'] = note
        data['strength'] = strength
        data['label'] = LABEL_LIST[dataAgeGroup-1]
        data['age'] = AGE_LIST[dataAgeGroup-1]
        data['stage'] = THEILER_STAGE_LIST[dataAgeGroup-1]

        #
        #  If there are images for the current tissue, strip off any
        #  leading/trailing spaces and create a comma-separated list to
        #  add to the dictionary.
        #
        images = ''
        if imageColumn > 0:
            for i in string.split(self.fields[imageColumn-1],','):
                if string.strip(i) != '':
                    if len(images) == 0:
                        images = string.strip(i)
                    else:
                        images = images + ',' + string.strip(i)
        data['images'] = images

        #
        #  Get the translation values for the current tissue.
        #
        translation = self.translator.getValues(transAgeGroup-1,tissue)

        if self.tissueNum == MAX_TISSUES:
            last = 1
        else:
            last = 0

        return data,translation,last


    # Purpose: Get the gene description from the current input record.
    # Returns: The gene description
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getGeneDesc (self):
        return self.geneDesc


    # Purpose: Get the symbol from the current input record.
    # Returns: The symbol
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getSymbol (self):
        return self.symbol


    # Purpose: Get the sequence ID from the current input record.
    # Returns: The sequence ID
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getSeqID (self):
        return self.seqID


    # Purpose: Get the probe accession ID from the current input record.
    # Returns: The MGI ID of the probe
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getProbeAccID (self):
        return self.probeAccID


    # Purpose: Get the probe name from the current input record.
    # Returns: The probe name
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getProbeName (self):
        return self.probeName


    # Purpose: Get the clone library from the current input record.
    # Returns: The clone library name
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing
    #
    def getCloneLib (self):
        return self.cloneLib

