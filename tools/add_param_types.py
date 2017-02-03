#!/usr/bin/env python2

"""
Takes stdin or a filename and parses FORTRAN parameters.  Adds explict types to implicitly-typed parameters

To use:
    $ parseSize.py location/of/SIZE > destination/of/new/SIZE
or
    $ cat location/of/SIZE | parseSize.py > destination/of/new/SIZE

This usage is intended to be UNIX-like and to prevent users from accidentally clobbering their SIZE file.

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

def getTypedParams(files=None):
    """ Get names of explictly-typed variables from a Nek5000 SIZE file

    This may likely work for an arbitrary FORTRAN 77 source file.

    :param files: May be a list of filenames to parse sequentially.
            * If None, then takes filenames from sys.argv[1:].
            * If None and if sys.argv[1:] is None, then takes contents of stdin
    :return: A set of variable names.  Does not return values.
    """
    declaredVars = set()
    for line in fileinput.input(files):
        line = line.lstrip()
        for m in re.finditer(r'^(?P<type>(?:integer)|(?:real)|(?:double precision)|(?:complex)|(?:logical)|(?:character))[/ \t]*?(?P<var>(?:[a-z0-9/ \t,])*)', line, re.I):
            varList = [v.lower().strip() for v in m.group('var').split(',')]
            declaredVars.update(varList)
    return declaredVars


def writeParamTypes(files=None):
    """ Prepends explicit type declarations to a Nek5000 SIZE file

    Parameters that are already explicitly-typed will be ignored; this should avoid any compiler errors.

    The edited SIZE file will be printed to stdout.  The user may redirect it on the command line.

    This may likely work for an arbitray FORTRAN 77 source file.

    :param files: May be a list of filenames to parse sequentially.
            * If None, then takes filenames from sys.argv[1:].
            * If None and if sys.argv[1:] is None, then takes contents of stdin
    """

    # Get parameter declarations and explicitly typed params
    allParams   = getSizeParams(files)   # A dict
    typedParams = getTypedParams(files)  # A tuple

    # Lists of variables that must be explicitly typed
    ints  = []
    reals = []

    for p in sorted(allParams.iterkeys()):
        if p not in typedParams:
            if p.startswith(('i','j','k','l','m','I','J','K','L','M','N')):
                ints.append(p)
            else:
                reals.append(p)

    print 'c automatically prepended by add_param_types.py\n'
    for i in ints:
        print('      integer {0}'.format(i))
    for r in reals:
        print('      real {0}'.format(r))
    print '\n'
    for line in fileinput.input(files):
        print line.rstrip()


if __name__ == '__main__':
    writeParamTypes()
