
#
# Program: gxdloadlib.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	Provide a set of library functions that any GXD data load can use.
#
# Requirements Satisfied by This Program:
#
# Usage:
#
# Envvars:
#
# Inputs:
#
# Outputs:
#
# Exit Codes:
#
# Assumes:
#
# Bugs:
#
# Implementation:
#
# 03/05/2014	lec
#	- TR11471/verifyPrepAntibody/verifyPrepSecondary
#
# 01/20/2010	lec
#	- TR9560; remove verifyPrepCoverage
#

import sys
import os
import accessionlib
import db

#globals

antibodyDict = {}       # antibody
assayTypeDict = {}      # assay type
coverageDict = {}       # probe coverage
embeddingDict = {}      # embedding method
fixationDict = {}       # fixation
gelRNATypeDict = {}     # gel rna types
gelControlDict = {}     # gel control
gelUnitsDict = {}       # gel units
gelStrengthDict = {}    # gel strength
genotypeDict = {}       # genotype
idxassayDict = {}	# index assay
idxpriorityDict = {}	# index priority
idxstageDict = {}	# index stage
labelDict = {}          # probe label
patternDict = {}        # pattern
reporterGeneDict = {}   # reporter gene
secondaryDict = {}      # antibody/secondary antibody
senseDict = {}          # probe sense
strengthDict = {}       # strength
visualDict = {}         # probe visualization

prepTypeList = ['DNA', 'RNA', 'Not Specified']  # probe prep types
hybridizationList = ['section', 'whole mount', 'section from whole mount', 'Not Specified']

# Purpose:  verify Antibody Accession ID
# Returns:  Antibody Key if Antibody is valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Antibody exists either in the Antibody dictionary or the database
#	writes to the error file if the Antibody is invalid
#	adds the Antibody id and key to the Antibody dictionary if the Antibody is valid
# Throws:  nothing

def verifyAntibody(
    antibodyID,	# Accession ID of the Antibody (str.
    lineNum,	# line number (integer)
    errorFile   # error file (file descriptor)
    ):

    global antibodyDict

    antibodyKey = 0

    if antibodyID in antibodyDict:
        return antibodyDict[antibodyID]
    else:
        results = db.sql('select _Object_key from ACC_Accession where _MGIType_key = 6 and accID = \'%s\' ' % (antibodyID), 'auto')

        for r in results:
            if r['_Object_key'] is None:
                if errorFile != None:
                    errorFile.write('Invalid Mouse Probe (%d) %s\n' % (lineNum, antibodyID))
                antibodyKey = 0
            else:
                antibodyKey = r['_Object_key']
                antibodyDict[antibodyID] = antibodyKey

    return antibodyKey

# Purpose:  verify Assay Type
# Returns:  Assay Type key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Assay Type exists in the Assay Type dictionary
#	writes to the error file if the Assay Type is invalid
# Throws:  nothing

def verifyAssayType(
    assayType, 	# Assay Type value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global assayTypeDict

    assayTypeKey = 0

    if len(assayTypeDict) == 0:
        results = db.sql('select _AssayType_key, assayType from GXD_AssayType', 'auto')
        for r in results:
            assayTypeDict[r['assayType']] = r['_AssayType_key']

    if assayType in assayTypeDict:
        assayTypeKey = assayTypeDict[assayType]
    else:
        if errorFile != None:
            errorFile.write('Invalid Assay Type (%d): %s\n' % (lineNum, assayType))

    return assayTypeKey

# Purpose:  verify Embedding Method
# Returns:  Embedding Method key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Embedding Method exists in the Embedding Method dictionary
#	writes to the error file if the Embedding Method is invalid
# Throws:  nothing

def verifyEmbeddingMethod(
    embedding, 	# Embedding Method value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global embeddingDict

    embeddingKey = 0

    if len(embeddingDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 155', 'auto')
        for r in results:
            embeddingDict[r['term']] = r['_Term_key']

    if embedding in embeddingDict:
        embeddingKey = embeddingDict[embedding]
    else:
        if errorFile != None:
            errorFile.write('Invalid Embedding Method (%d): %s\n' % (lineNum, embedding))

    return embeddingKey

# Purpose:  verify Fixation Method
# Returns:  Fixation key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Fixation Method exists in the Fixation dictionary
#	writes to the error file if the Fixation Method is invalid
# Throws:  nothing

def verifyFixationMethod(
    fixation, 	# Fixation Method value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global fixationDict

    fixationKey = 0

    if len(fixationDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 156', 'auto')
        for r in results:
            fixationDict[r['term']] = r['_Term_key']

    if fixation in fixationDict:
        fixationKey = fixationDict[fixation]
    else:
        if errorFile != None:
            errorFile.write('Invalid Fixation (%d): %s\n' % (lineNum, fixation))

    return fixationKey

# Purpose:  verify Gel RNA Type
# Returns:  Gel RNA Type key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel RNA Type exists in the Gel RNA Type dictionary
#	writes to the error file if the Gel RNA Type is invalid
# Throws:  nothing

def verifyGelRNAType(
    gelRNAType, 	# Gel RNA Type value (str.
    lineNum,	# line number (integer)
    errorFile	# error file (file descriptor)
    ):

    global gelRNATypeDict

    gelRNATypeKey = 0

    if len(gelRNATypeDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 172', 'auto')
        for r in results:
            gelRNATypeDict[r['term']] = r['_Term_key']

    if gelRNAType in gelRNATypeDict:
        gelRNATypeKey = gelRNATypeDict[gelRNAType]
    else:
        if errorFile != None:
            errorFile.write('Invalid Gel RNA Type (%d): %s\n' % (lineNum, gelRNAType))

    return gelRNATypeKey

# Purpose:  verify Gel Control
# Returns:  Gel Control key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Control exists in the Gel Control dictionary
#	writes to the error file if the Gel Control is invalid
# Throws:  nothing

def verifyGelControl(
    gelControl, 	# Gel Control value (str.
    lineNum,	# line number (integer)
    errorFile	# error file (file descriptor)
    ):

    global gelControlDict

    gelControlKey = 0

    if len(gelControlDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 154', 'auto')
        for r in results:
            gelControlDict[r['term']] = r['_Term_key']

    if gelControl in gelControlDict:
        gelControlKey = gelControlDict[gelControl]
    else:
        if errorFile != None:
            errorFile.write('Invalid Gel Control (%d): %s\n' % (lineNum, gelControl))

    return gelControlKey

# Purpose:  verify Gel Units
# Returns:  Gel Units key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Units exists in the Gel Units dictionary
#	writes to the error file if the Gel Units is invalid
# Throws:  nothing

def verifyGelUnits(
    gelUnits, 	# Gel Units value (str.
    lineNum,	# line number (integer)
    errorFile	# error file (file descriptor)
    ):

    global gelUnitsDict

    gelUnitsKey = 0

    if len(gelUnitsDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 173', 'auto')
        for r in results:
            gelUnitsDict[r['term']] = r['_Term_key']

    if gelUnits in gelUnitsDict:
        gelUnitsKey = gelUnitsDict[gelUnits]
    else:
        if errorFile != None:
            errorFile.write('Invalid Gel Units (%d): %s\n' % (lineNum, gelUnits))

    return gelUnitsKey

# Purpose:  verifies the genotype
# Returns:  the primary key of the genotype or 0 if invalid
# Assumes:  nothing
# Effects:  verifies that the Genotype exists by checking the genotypeDict
#	dictionary for the Genotype ID or the database.
#	writes to the error file if the Genotype is invalid.
#	adds the Genotype ID/Key to the global genotypeDict dictionary if the
#	genotype is valid.
# Throws:

def verifyGenotype(
    genotypeID,          # genotype accession ID; MGI:#### (str.
    lineNum,		 # line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global genotypeDict

    if genotypeID in genotypeDict:
        genotypeKey = genotypeDict[genotypeID]
    else:
        genotypeKey = accessionlib.get_Object_key(genotypeID, 'Genotype')
        if genotypeKey is None:
            if errorFile != None:
                errorFile.write('Invalid Genotype (%d): %s\n' % (lineNum, genotypeID))
            genotypeKey = 0
        else:
            genotypeDict[genotypeID] = genotypeKey

    return(genotypeKey)

# Purpose:  verify Hybridization Type
# Returns:  1 if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Hybridization exists in the hybridization list
#	writes to the error file if the Hybridization is invalid
# Throws:  nothing

def verifyHybridization(
    hybridization, # value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global hybridizationList

    if hybridization in hybridizationList:
        return 1
    else:
        if errorFile != None:
            errorFile.write('Invalid Hybridization (%d): %s\n' % (lineNum, hybridization))
        return 0

# Purpose:  verify Index Assay
# Returns:  Assay key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Assay exists in the idxassay dictionary
#	writes to the error file if the Assay is invalid
# Throws:  nothing

def verifyIdxAssay(
    idxassay, 	# Assay value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global idxassayDict

    idxassayKey = 0

    if len(idxassayDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 12', 'auto')
        for r in results:
            idxassayDict[r['term']] = r['_Term_key']

    if idxassay in idxassayDict:
        idxassayKey = idxassayDict[idxassay]
    else:
        if errorFile != None:
            errorFile.write('Invalid Index Assay (%d): %s\n' % (lineNum, idxassay))

    return idxassayKey

# Purpose:  verify Index Priority
# Returns:  Priority key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Priority exists in the idxpriority dictionary
#	writes to the error file if the Priority is invalid
# Throws:  nothing

def verifyIdxPriority(
    idxpriority, 	# Priority value (str.
    lineNum,		# line number (integer)
    errorFile	   	# error file (file descriptor)
    ):

    global idxpriorityDict

    idxpriorityKey = 0

    if len(idxpriorityDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 11', 'auto')
        for r in results:
            idxpriorityDict[r['term']] = r['_Term_key']

    if idxpriority in idxpriorityDict:
        idxpriorityKey = idxpriorityDict[idxpriority]
    else:
        if errorFile != None:
            errorFile.write('Invalid Index Priority (%d): %s\n' % (lineNum, idxpriority))

    return idxpriorityKey

# Purpose:  verify Index Stage
# Returns:  Stage key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Stage exists in the idxstage dictionary
#	writes to the error file if the Stage is invalid
# Throws:  nothing

def verifyIdxStage(
    idxstage, 	# Stage value (str.
    lineNum,		# line number (integer)
    errorFile	   	# error file (file descriptor)
    ):

    global idxstageDict

    idxstageKey = 0

    if len(idxstageDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 13', 'auto')
        for r in results:
            idxstageDict[r['term']] = r['_Term_key']

    if idxstage in idxstageDict:
        idxstageKey = idxstageDict[idxstage]
    else:
        if errorFile != None:
            errorFile.write('Invalid Index Stage (%d): %s\n' % (lineNum, idxstage))

    return idxstageKey

# Purpose:  verify Probe Prep Label
# Returns:  Probe Prep Label key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Label exists in the Label dictionary
#	writes to the error file if the Prep Label is invalid
# Throws:  nothing

def verifyPrepLabel(
    label, 	# Label value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global labelDict

    labelKey = 0

    if len(labelDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 152', 'auto')
        for r in results:
            labelDict[r['term']] = r['_Term_key']

    if label in labelDict:
        labelKey = labelDict[label]
    else:
        if errorFile != None:
            errorFile.write('Invalid Prep Label (%d): %s\n' % (lineNum, label))

    return labelKey

# Purpose:  verify Antibody Prep Secondary
# Returns:  Antibody Prep Secondary key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Secondary exists in the secondary dictionary
#	writes to the error file if the Prep Secondary is invalid
# Throws:  nothing

def verifyPrepSecondary(
    secondary, 	# Secondary value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global secondaryDict

    secondaryKey = 0

    if len(secondaryDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 160', 'auto')
        for r in results:
            secondaryDict[r['term']] = r['_Term_key']

    if secondary in secondaryDict:
        secondaryKey = secondaryDict[secondary]
    else:
        if errorFile != None:
            errorFile.write('Invalid Prep Secondary (%d): %s\n' % (lineNum, secondary))

    return secondaryKey

# Purpose:  verify Probe Prep Sense
# Returns:  Probe Prep Sense key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Sense exists in the sense dictionary
#	writes to the error file if the Prep Sense is invalid
# Throws:  nothing

def verifyPrepSense(
    sense, 	# Sense value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global senseDict

    senseKey = 0

    if len(senseDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 159', 'auto')
        for r in results:
            senseDict[r['term']] = r['_Term_key']

    if sense in senseDict:
        senseKey = senseDict[sense]
    else:
        if errorFile != None:
            errorFile.write('Invalid Prep Sense (%d): %s\n' % (lineNum, sense))

    return senseKey

# Purpose:  verify Probe Prep Type
# Returns:  1 if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Type exists in the probeType list
#	writes to the error file if the Prep Type is invalid
# Throws:  nothing

def verifyPrepType(
    prepType, 	# Type value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global prepTypeList

    if prepType in prepTypeList:
        return 1
    else:
        if errorFile != None:
            errorFile.write('Invalid Prep Type (%d): %s\n' % (lineNum, prepType))
        return 0

# Purpose:  verify Probe Prep Visualization
# Returns:  Probe Prep Visualization key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Prep Visualization exists in the Visualization dictionary
#	writes to the error file if the Prep Visualization is invalid
# Throws:  nothing

def verifyPrepVisualization(
    visualization, 	# Visualization value (str.
    lineNum,		# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global visualDict

    visualKey = 0

    if len(visualDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 157', 'auto')
        for r in results:
            visualDict[r['term']] = r['_Term_key']

    if visualization in visualDict:
        visualKey = visualDict[visualization]
    else:
        if errorFile != None:
            errorFile.write('Invalid Prep Visualization (%d): %s\n' % (lineNum, visualization))

    return visualKey

# Purpose:  verify Gel Strength
# Returns:  Gel Strength key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Gel Strength exists in the Gel Strength dictionary
#	writes to the error file if the Gel Strength is invalid
# Throws:  nothing

def verifyGelStrength(
    gelStrength, 	# Gel Strength value (str.
    lineNum,	# line number (integer)
    errorFile	# error file (file descriptor)
    ):

    global gelStrengthDict

    gelStrengthKey = 0

    if len(gelStrengthDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 163', 'auto')
        for r in results:
            gelStrengthDict[r['term']] = r['_Term_key']

    if gelStrength in gelStrengthDict:
        gelStrengthKey = gelStrengthDict[gelStrength]
    else:
        if errorFile != None:
            errorFile.write('Invalid Gel Strength (%d): %s\n' % (lineNum, gelStrength))

    return gelStrengthKey

# Purpose:  verify Strength
# Returns:  Strength key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Strength exists in the Strength dictionary
#	writes to the error file if the Strength is invalid
# Throws:  nothing

def verifyStrength(
    strength, 	# Strength value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global strengthDict

    strengthKey = 0

    if len(strengthDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 163', 'auto')
        for r in results:
            strengthDict[r['term']] = r['_Term_key']

    if strength in strengthDict:
        strengthKey = strengthDict[strength]
    else:
        if errorFile != None:
            errorFile.write('Invalid Strength (%d): %s\n' % (lineNum, strength))

    return strengthKey

# Purpose:  verify Pattern
# Returns:  Pattern key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Pattern exists in the Pattern dictionary
#	writes to the error file if the Pattern is invalid
# Throws:  nothing

def verifyPattern(
    pattern, 	# Pattern value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global patternDict

    patternKey = 0

    if len(patternDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 153', 'auto')
        for r in results:
            patternDict[r['term']] = r['_Term_key']

    if pattern in patternDict:
        patternKey = patternDict[pattern]
    else:
        if errorFile != None:
            errorFile.write('Invalid Pattern (%d): %s\n' % (lineNum, pattern))

    return patternKey

# Purpose:  verify Reporter Gene
# Returns:  Reporter Gene key if valid, else 0
# Assumes:  nothing
# Effects:  verifies that the Reporter Gene exists in the Reporter Gene dictionary
#	writes to the error file if the Reporter Gene is invalid
# Throws:  nothing

def verifyReporterGene(
    reporterGene, # Reporter Gene value (str.
    lineNum,	# line number (integer)
    errorFile	   # error file (file descriptor)
    ):

    global reporterGeneDict

    reporterGeneKey = 0

    if len(reporterGeneDict) == 0:
        results = db.sql('select _Term_key, term from VOC_Term where _Vocab_key = 14', 'auto')
        for r in results:
            reporterGeneDict[r['term']] = r['_Term_key']

    if reporterGene in reporterGeneDict:
        reporterGeneKey = reporterGeneDict[reporterGene]
    else:
        if errorFile != None:
            errorFile.write('Invalid Reporter Gene (%d): %s\n' % (lineNum, reporterGene))

    return reporterGeneKey
