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
 
cat - <<EOSQL | doisql.csh $0 >> $LOG

use $DBNAME
go

declare @refsKey integer
select @refsKey = _Object_key from BIB_Acc_View where accID = "${JNUM}"

select p._Probe_key, seq = identity(5)
into #newpref
from PRB_Probe p, GXD_ProbePrep pp, GXD_Assay a
where a._Refs_key = @refsKey
and a._ProbePrep_key = pp._ProbePrep_key
and pp._Probe_key = p._Probe_key
and not exists (select 1 from PRB_Reference r
where p._Probe_key = r._Probe_key
and r._Refs_key = @refsKey)

declare @maxRefKey integer
select @maxRefKey = max(_Reference_key) + 1 from PRB_Reference

insert into PRB_Reference
select seq + @maxRefKey, _Probe_key, @refsKey, null, 0, 0, getdate(), getdate()
from #newpref
go

checkpoint
go

end

EOSQL

date >> $LOG

