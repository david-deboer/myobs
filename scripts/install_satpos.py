#! /usr/bin/env python
from os import path
import shutil


shutil.copy('src/satpos', path.dirname(__file__))
