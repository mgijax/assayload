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

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
cd ${RTPCRDATADIR}

${ASSAYLOADINSTALLDIR}/J80501rtpcr.py
${ASSAYLOADINSTALLDIR}/gelload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG

cd ${INSITUDATADIR}

${ASSAYLOADINSTALLDIR}/J80501rnainsitu.py
${ASSAYLOADINSTALLDIR}/insituload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG

cd `dirname $0`

${ASSAYLOADINSTALLDIR}/indexload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} -RJ:80501 >>& $LOG
${ASSAYLOADINSTALLDIR}/createPrbReference.csh J:80501

date >> $LOG

