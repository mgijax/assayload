#!/bin/csh -f

#
# Create Probe Reference records for new Reference
#

cd `dirname $0` && source ./Configuration

setenv JNUM $1

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date >> $LOG
 
# must be converted to postgres

cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

declare
v_maxRefKey := select max(_Reference_key) + 1 from PRB_Reference ;

CREATE TEMP TABLE newpref ON COMMIT DROP
AS SELECT row_number() over (ORDER BY _Probe_key) as seq,
from PRB_Probe p, GXD_ProbePrep pp, GXD_Assay a
where a._Refs_key = (select _Object_key from BIB_Acc_View where accID = '${JNUM}')
and a._ProbePrep_key = pp._ProbePrep_key
and pp._Probe_key = p._Probe_key
and not exists (select 1 from PRB_Reference r where p._Probe_key = r._Probe_key and r._Refs_key = v_refsKey)
;

insert into PRB_Reference
select seq + v_maxRefKey, 
_Probe_key, 
(select _Refs_key from BIB_Citation_Cache where jnumID = '${JNUM}'), null, 0, 0, now(), now()
from newpref
;

EOSQL

date >> $LOG

