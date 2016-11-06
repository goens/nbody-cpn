#!/bin/bash
if [ $# -ne 1 ]; then
	 echo "usage: $0 N"
	 exit 1
fi
i=$1
sed -i "s/PNgeneric(NUM_WORKERS,[0-9]*/PNgeneric(NUM_WORKERS,$i/" nbody.cpnm4
