#!/bin/csh -f

#
# TR 7982
#
# Wrapper script for loading TR 7982 data into GXD Assay
#
# Processing:
#	1. Load InSitus
#	2. Generate Index record
#	3. Load Images (gxdimageload)
#	4. Load InSitu/Image associations (gxdimageload)
#

cd `dirname $0` && source ./tr7982.config

./tr7982.py
./insituload.py

${MRKCACHELOAD}/mrkref.csh

