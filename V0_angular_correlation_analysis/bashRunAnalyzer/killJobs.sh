#!/bin/bash
for i in {1..3600}
do
  qselect -u jdeclerc | xargs qdel
  sleep 2
  qstat -u jdeclerc 
  sleep 5
done
