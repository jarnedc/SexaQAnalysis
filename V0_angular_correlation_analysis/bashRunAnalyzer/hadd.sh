#!/bin/bash

array=($(ls -d ${INPUT_PATH}))

echo hadd -f  ${OUTPUTFILENAME} ${INPUT_PATH}*.root
