#!/bin/bash
cd ~/myobs
CURRENTDATE=`date +"%Y-%m-%d %T"`
GCOMMSG='git commit "TLE update '$CURRENTDATE'"'
eval $GCOMMSG
git push origin main
