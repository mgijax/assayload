#!/bin/csh -f

#
# TR 9365
#
# Wrapper script for loading data into GXD Assay
#
# Processing:
#	1. Load Probes
#	2. Load RTPCR
#	3. Generate Index record
#

cd `dirname $0` && source ./Configuration

source ./tr9365.config

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
cat - <<EOSQL | doisql.csh $MGD_DBSERVER $MGD_DBNAME $0 | tee -a $LOG

use ${MGD_DBNAME}
go

update BIB_Refs set journal = 'MGI Direct Data Submission' where _Refs_key = 141558
go

checkpoint
go

quit

EOSQL

# load primers
${PROBELOAD}/primerload.csh tr9365.config | tee -a $LOG

# put this into gel load format
${ASSAYLOAD}/tr9365rtpcr.py | tee -a $LOG

# run gelload
cd ${RTPCRDATADIR}
${ASSAYLOAD}/gelload.py | tee -a $LOG

# load index
cd `dirname $0`
${ASSAYLOAD}/indexload.py | tee -a $LOG

# update marker cache
#${MRKCACHELOAD}/mrkref.csh | tee -a $LOG

date >> $LOG

