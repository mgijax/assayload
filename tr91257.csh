#!/bin/csh -f

#
# TR 91257
#
# Wrapper script for loading TR 6118 data into GXD Assay
#
# Processing:
#	1. Load Probes
#	2. Load InSitus
#	3. Generate Index record
#	4. Load Images (gxdimageload)
#	5. Load InSitu/Image associations (gxdimageload)
#

cd `dirname $0` && source ./Configuration

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
cd ${INSITUDATADIR}

#${ASSAYLOADINSTALLDIR}/J91257probe.py
#${PROBELOAD} -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} -Iprobe.txt >>& $LOG

${ASSAYLOADINSTALLDIR}/J91257insitu.py
${ASSAYLOADINSTALLDIR}/insituload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG

cd `dirname $0`

#${ASSAYLOADINSTALLDIR}/indexload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} -RJ:80501 >>& $LOG

date >> $LOG

