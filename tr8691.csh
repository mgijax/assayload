#!/bin/csh -f

#
# TR 8691
#
# Wrapper script for loading TR 8691 data into GXD Index
#

cd `dirname $0`

# Freeman

source ./tr8691-freeman.config
cd ${LOGDIR}

setenv LOG ${LOGDIR}/$0.log
rm -rf ${LOG}
touch ${LOG}
 
date > ${LOG}
 
${ASSAYLOAD}/indexload.py >>& ${LOG}

# Deltagen

source ./tr8691-deltagen.config
cd ${LOGDIR}

setenv LOG ${LOGDIR}/$0.log
rm -rf ${LOG}
touch ${LOG}
 
date > ${LOG}
 
${ASSAYLOAD}/indexload.py >>& ${LOG}

# Marker Cache re-load

${MRKCACHELOAD}/mrkref.csh >>& ${LOG}

date >> ${LOG}

