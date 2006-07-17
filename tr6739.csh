#!/bin/csh -f

#
# TR 6739
#
# Wrapper script for loading TR 6739 data into GXD
#
# Processing:
#	1. Load Probes
#	2. Load InSitus
#	3. Generate Index record
#	4. Reload MRK_Reference table
#

cd `dirname $0` && source ./Configuration

setenv LOG ${LOGDIR}/`basename $0`.log
rm -rf $LOG
 
date > $LOG

cd ${DATADIR}

#
#  Create the temp tables.
#
echo "" >>& $LOG
echo "Create the temp tables" >>& $LOG
cat - <<EOSQL | doisql.csh ${MGD_DBSERVER} ${MGD_DBNAME} $0 >> $LOG

use tempdb
go

create table TMP_Probe (lineNum int not null,
                        symbol varchar(50) null,
                        seqID varchar(30) null,
                        probeAccID varchar(30) null,
                        probeName varchar(40) null,
                        cloneLibrary varchar(255) null)
go
grant all on TMP_Probe to public
go

create table TMP_ProbeMarker (lineNum int not null,
                              markerAccID varchar(30) not null,
                              probeAccID varchar(30) null,
                              probeName varchar(40) null)
go
grant all on TMP_ProbeMarker to public
go

create table TMP_QC (lineNum int not null,
                     process varchar(20) not null)
go
grant all on TMP_QC to public
go

quit
EOSQL
 

#cat - <<EOSQL | doisql.csh $0 >> $LOG
#
#use tempdb
#go
#
#truncate table TMP_Probe
#go
#truncate table TMP_ProbeMarker
#go
#truncate table TMP_QC
#go
#
#quit
#EOSQL
 

#
#  Create the input file for the probe load.
#
echo "" >>& $LOG
echo "Call tr6739probe.py" >>& $LOG
${ASSAYLOADINSTALLDIR}/tr6739probe.py >>& $LOG


#
#  Call the probe load.
#
echo "" >>& $LOG
echo "Call probeload.py" >>& $LOG
${PROBELOAD} -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${LOADMODE} -Iprobe.txt >>& $LOG


#
#  Create the input files for the in situ load.
#
echo "" >>& $LOG
echo "Call tr6739insitu.py" >>& $LOG
${ASSAYLOADINSTALLDIR}/tr6739insitu.py >>& $LOG


#
#  Call the in situ load.
#
echo "" >>& $LOG
echo "Call insituload.py" >>& $LOG
${ASSAYLOADINSTALLDIR}/insituload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${LOADMODE} >>& $LOG


#
#  Drop the temp tables.
#
echo "" >>& $LOG
echo "Drop the temp tables" >>& $LOG
cat - <<EOSQL | doisql.csh $0 >> $LOG

use tempdb
go

drop table TMP_Probe
go

drop table TMP_ProbeMarker
go

drop table TMP_QC
go

quit
EOSQL


#
#  Call the index load.
#
echo "" >>& $LOG
echo "Call indexload.py" >>& $LOG
${ASSAYLOADINSTALLDIR}/indexload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${LOADMODE} -RJ:93300 >>& $LOG


#
#  Reload the MRK_Reference table.
#
echo "" >>& $LOG
echo "Reload the MRK_Reference table" >>& $LOG
${MRKREFLOAD} >>& ${LOG}

date >> $LOG
