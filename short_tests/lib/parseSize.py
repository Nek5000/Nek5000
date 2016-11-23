#!/usr/bin/env python2

"""
Takes stdin or a filename and parses FORTRAN variable assignments.  Prints assigned variables to stdout.  This is intended to be used for

To use:
    $ parseSize.py location/of/SIZE
or
    $ cat location/of/SIZE | parseSize.py

"""
import fileinput
import re

def getSizeParams(files=None):
    """ Parse variable assignments from a Nek5000 SIZE file.

    This is intended to be used for Nek5000 SIZE files but may work for other FORTRAN source files.

    :param files: May be a list of filenames to parse sequentially.
            * If None, then takes filenames from sys.argv[1:].
            * If None and if sys.argv[1:] is None, then takes contents of stdin
    :return:  A dict of {variable : value} pairs.  All variables are forced to lowercase.  All values are strings.
    """
    params = {}
    for line in fileinput.input(files):
        if not line.startswith(('c', 'C')):
            for m in re.finditer(r'\b(?P<name>\w+)\b *= *(?P<value>(?:[a-z0-9\+\-\*/ \t]+|\([a-z0-9\+\-\*/ \t]+\))+)', line, re.I):
                params[m.group('name').lower().strip()] = m.group('value').lower().strip()
    return params

if __name__ == '__main__':
    params = getSizeParams()

    # Print all params in alphabetical order.
    # for x in sorted(params):
    #     print "{0} = {1}".format(x, params[x])

    # Print only variables in SIZE.template in the order that they appear in SIZE.template
    template_params = [
                'ldim',
                'lx1',
                'lxd',
                'lx2',
                'lx1m',
                'lelg',
                'lp',
                'lelt',
                'ldimt',
                'lelx',
                'lely',
                'lelz',
                'ax1',
                'ax2',
                'lbx1',
                'lbx2',
                'lbelt',
                'lpx1',
                'lpx2',
                'lpelt',
                'lpert',
                'lelecmt',
                'toteq',
                'lcvx1',
                'lcvelt',
                'mxprev',
                'lgmres',
                'lorder',
                'lhis',
                'maxobj',
                'maxmbr',
                'nsessmax',
                'nmaxl',
                'nfldmax',
                'nmaxcom',
    ]

    for x in template_params:
        print "{:<9}= '{}',".format(x, params.get(x, None))


