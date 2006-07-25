#!/bin/csh -f

#
# TR 7819
#
# Wrapper script for loading TR 7819 data into GXD Assay
#
# Processing:
#	1. Load InSitus
#	2. Generate Index record
#	3. Load Images (gxdimageload)
#	4. Load InSitu/Image associations (gxdimageload)
#

cd `dirname $0` && source ./tr7819.config

./tr7819.py
./insituload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${LOADMODE}
#./indexload.py -S${MGD_DBSERVER} -D${MGD_DBNAME} -U${MGD_DBUSER} -P${MGD_DBPASSWORDFILE} -M${LOADMODE} -RJ:101679

#${MRKREFLOAD}

