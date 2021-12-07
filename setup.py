#! /usr/bin/env python
# -*- mode: python; coding: utf-8 -*-
# Copyright 2019 David DeBoer
# Licensed under the 2-clause BSD license.

from setuptools import setup
import glob

setup_args = {
    'name': "myobs",
    'description': "various routines for observing etc",
    'license': "BSD",
    'author': "David DeBoer",
    'author_email': "ddeboer@berkeley.edu",
    'version': '0.1',
    'scripts': glob.glob('scripts/*'),
    'packages': ['myobs'],
    # 'package_data': {"myobs": ["data/*"]}
}

if __name__ == '__main__':
    setup(**setup_args)
