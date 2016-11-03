#!/usr/bin/env python2

"""
Takes stdin or a filenameand parses paramemter assignments.  Outputs parameters required by SIZE.template.

To use:
    $ cat SIZE | parseSize.py
or
    $ parseSize.py SIZE

"""


import fileinput
import re
params = {}

for line in fileinput.input():
    if not line.startswith(('c', 'C')):
        for m in re.finditer(r'\b(?P<name>\w+)\b *= *(?P<value>(?:[a-z0-9\+\-\*/ \t]+|\([a-z0-9\+\-\*/ \t]+\))+)', line, re.I):
            params[m.group('name').lower().strip()] = m.group('value').lower().strip()

# print all
for x in sorted(params):
   print "{0} = {1}".format(x, params[x])

# print only variables in size.template
# template_params = [
#             'ldim',
#             'lx1',
#             'lxd',
#             'lx2',
#             'lx1m',
#             'lelg',
#             'lp',
#             'lelt',
#             'ldimt',
#             'lelx',
#             'lely',
#             'lelz',
#             'ax1',
#             'ax2',
#             'lbx1',
#             'lbx2',
#             'lbelt',
#             'lpx1',
#             'lpx2',
#             'lpelt',
#             'lpert',
#             'lelecmt',
#             'toteq',
#             'mxprev',
#             'lgmres',
#             'lorder',
#             'lhis',
#             'maxobj',
#             'maxmbr',
#             'nsessmax',
#             'nmaxl',
#             'nfldmax',
#             'nmaxcom',
# ]
#
# for x in template_params:
#     print "{0} = {1}".format(x, params.get(x, None))


