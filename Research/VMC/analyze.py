#!/usr/bin/env python

import VMC_tools as H
import sys, os

where = sys.argv[1]
H.read_raw_data(where)
H.plot_data(where)

