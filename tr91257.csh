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

setenv REFERENCE J:91257
setenv PROBELOADINPUT probe.txt

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date > $LOG
 
#cat - <<EOSQL | doisql.csh ${MGD_DBSERVER} ${MGD_DBNAME} $0 >> $LOG

#use ${MGD_DBNAME}
#go

#update BIB_Refs set journal = 'MGI Direct Data Submission' where _Refs_key = 92242
#go

#declare @key integer
#select @key = max(_Source_key) + 1 from PRB_Source
#insert into PRB_Source 
#values(@key,63468,316367,1,433,450,315167,316335,null,null,null,'embryonic day 13.5',13.5,13.5,0,1000,1000,getdate(),getdate())
#go

#declare @key integer
#select @key = max(_Source_key) + 1 from PRB_Source
#insert into PRB_Source 
#values(@key,63468,316370,1,-1,-1,315167,316335,null,null,null,'Not Specified',-1.0,-1.0,0,1000,1000,getdate(),getdate())
#go

#checkpoint
#go
#
#quit

#EOSQL

cd ${INSITUDATADIR}

# probes

${ASSAYLOAD}/J91257probe.py
${PROBELOAD}/probeload.py >>& $LOG

# 1st set (section)

${ASSAYLOAD}/J91257insitu.py
${ASSAYLOAD}/insituload.py >>& $LOG

# 2nd set (whole mount)

cd ${WMDATADIR}
setenv INSITUDATADIR ${WMDATADIR}

${ASSAYLOAD}/J91257insituWM.py
${ASSAYLOAD}/insituload.py >>& $LOG

cd `dirname $0`
${ASSAYLOAD}/indexload.py >>& $LOG

${MRKCACHELOAD}/mrkref.csh >>& ${LOG}

date >> $LOG

