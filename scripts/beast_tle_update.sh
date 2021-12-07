#!/bin/bash

CURRENTDATE=`date +"%Y-%m-%d %T"`
GCOMMSG='git commit -C ~/myobs "TLE update '$CURRENTDATE'"'
eval $GCOMMSG
git push origin main
