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
 
cat - <<EOSQL | doisql.csh ${MGD_DBSERVER} ${MGD_DBNAME} $0 >> $LOG

use ${MGD_DBNAME}
go

declare @refsKey integer
select @refsKey = _Object_key from BIB_Acc_View where accID = "${JNUM}"

update PRB_Marker
set relationship = "E"
from PRB_Marker p, GXD_ProbePrep pp, GXD_Assay a
where a._Refs_key = @refsKey
and a._ProbePrep_key = pp._ProbePrep_key
and pp._Probe_key = p._Probe_key
and a._Marker_key = p._Marker_key
and p.relationship = "P"
go

checkpoint
go

end

EOSQL

date >> $LOG

