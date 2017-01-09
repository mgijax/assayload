#!/bin/csh -f

#
# Update Putative Relationships to Encodes
# between Probes/Markers for given GXD Assays (via Probe Prep).
#

cd `dirname $0` && source ./Configuration

setenv JNUM $1

setenv LOG $0.log
rm -rf $LOG
touch $LOG
 
date >> $LOG
 
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a $LOG

update PRB_Marker
set relationship = 'E'
from PRB_Marker p, GXD_ProbePrep pp, GXD_Assay a
where a._Refs_key = (select _Refs_key from BIB_Citation_Cache where jnumID = '${JNUM}')
and a._ProbePrep_key = pp._ProbePrep_key
and pp._Probe_key = p._Probe_key
and a._Marker_key = p._Marker_key
and p.relationship = 'P'
;

EOSQL

date >> $LOG

