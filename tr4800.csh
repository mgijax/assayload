#!/bin/csh -f

#
# TR 4800
#
# Wrapper script for loading TR 4800 data into GXD Assay
#
# Processing:
#	1. Load Probes (manual)
#	2. Load RT-PCR (includes primers)
#	3. Load InSitus
#	4. Load Images
#	5. Load InSitu/Image associations
#

cd `dirname $0` && source ./Configuration

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
rtpcr.py
gelload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG
exit 0

setenv INSITUDATALOAD $INSITU10DATALOAD
setenv INSITUDATA 10
rnainsitu.py
insituload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG
#setenv INSITUDATALOAD $INSITU14DATALOAD
#setenv INSITUDATA 14
#rnainsitu.py
#rnainsitu.py
#insituload.py -S${DBSERVER} -D${DBNAME} -U${DBUSER} -P${DBPASSWORDFILE} -M${LOADMODE} >>& $LOG

date >> $LOG

