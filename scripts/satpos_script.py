#! /usr/bin/env python
from myobs import satcensus
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('tlefile', help="Name of tle file for which to generate satpos script.")
args = ap.parse_args()

satcensus.satpos_script(args.tlefile)
