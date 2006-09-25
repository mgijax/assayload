#!/bin/csh -f

#
# TR 5154
#
# Wrapper script for loading TR 5154 data into GXD Assay
#
# Processing:
#	1. Load Probes (manual)
#	2. Load RT-PCR (includes primers)
#	3. Load InSitus
#	4. Generate Index record
#	5. Load Images (gxdimageload)
#	6. Load InSitu/Image associations (gxdimageload)
#

cd `dirname $0` && source ./Configuration

setenv REFERENCE J:80501

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
cd ${RTPCRDATADIR}

${ASSAYLOAD}/J80501rtpcr.py
${ASSAYLOAD}/gelload.py >>& $LOG

cd ${INSITUDATADIR}

${ASSAYLOAD}/J80501rnainsitu.py
${ASSAYLOAD}/insituload.py >>& $LOG

cd `dirname $0`

${ASSAYLOAD}/indexload.py >>& $LOG
${ASSAYLOAD}/createPrbReference.csh J:80501

date >> $LOG

