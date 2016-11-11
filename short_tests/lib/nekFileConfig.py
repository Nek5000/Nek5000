import os, stat, re

def config_makenek(infile, outfile, source_root=None, f77=None, cc=None, ifmpi=None, pplist=None, usr_lflags=None):

    with open(infile, 'r') as f:
        lines = f.readlines()

    if source_root:
        lines = [re.sub(r'^SOURCE_ROOT=\"+.+?\"+', r'SOURCE_ROOT="{0}"'.format(source_root), l)
                 for l in lines]
    if f77:
        lines = [re.sub(r'^F77=\"+.+?\"+', r'F77="{0}"'.format(f77), l)
                 for l in lines]
    if cc:
        lines = [re.sub(r'^CC=\"+.+?\"+', r'CC="{0}"'.format(cc), l)
                 for l in lines]
    if ifmpi:
        lines = [re.sub(r'^#*IFMPI=\"+.+?\"+', r'IFMPI="{0}"'.format(ifmpi), l)
                 for l in lines]

    if pplist:
        lines = [re.sub(r'^#*PPLIST=\"+.+?\"+', r'PPLIST="{0}"'.format(pplist), l)
                 for l in lines]

    if usr_lflags:
        lines = [re.sub(r'^#*USR_LFLAGS=\"+.+?\"+', r'USR_LFLAGS="{0}"'.format(usr_lflags), l)
                 for l in lines]

    lines = [re.sub(r'(^source\s+\$SOURCE_ROOT/makenek.inc)', r'\g<1> >compiler.out', l)
             for l in lines]

    lines = [re.sub(r'(.+)2>&1\s+\|\s*tee\s+compiler.out', r'\g<1>', l)
             for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)
    os.chmod(outfile,
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH |
             stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |
             stat.S_IWUSR )


def config_maketools(infile, outfile, f77=None, cc=None, bigmem=None):

    with open(infile, 'r') as f:
        lines = f.readlines()

    if f77:
        lines = [re.sub(r'^F77=\"+.+?\"+', r'F77="{0}"'.format(f77), l)
                 for l in lines]
    if cc:
        lines = [re.sub(r'^CC=\"+.+?\"+', r'CC="{0}"'.format(cc), l)
                 for l in lines]
    if bigmem:
        lines = [re.sub(r'BIGMEM=\"+.+?\"+', r'BIGMEM="{0}"'.format(bigmem), l)
                 for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)
    os.chmod(outfile,
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH |
             stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |
             stat.S_IWUSR )


def config_basics_inc(infile, outfile, nelm):

    with open(infile, 'r') as f:
        lines = f.readlines()

    lines = [re.sub(r'(.*nelm *= *)[ 0-9]+(.*)', r'\g<1>{0}\g<2>'.format(nelm), l, flags=re.I)
             for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)


def config_size( infile, outfile, **kwargs ):

    with open(infile, 'r') as f:
        lines = f.readlines()

    # Substitute all the variables
    for key, value in kwargs.iteritems():
        if value:
            lines = [
                re.sub(
                    r'(.*\bparameter\b.*\b{0} *= *)\S+?( *[),])'.format(key),
                    r'\g<1>{0}\g<2>'.format(value), l, flags=re.I)
                for l in lines ]

    with open(outfile, 'w') as f:
        f.writelines(lines)


def config_parfile(infile, outfile, general=None, problemtype=None, mesh=None, pressure=None, velocity=None,
                   temperature=None, cvode=None):
    import ConfigParser

    parfile = ConfigParser.SafeConfigParser()
    parfile.read(infile)

    if general:
        for k,v in general.iteritems():
            parfile.set('GENERAL', k, v)

    if problemtype:
        for k,v in problemtype.iteritems():
            parfile.set('PROBLEMTYPE', k, v)

    if mesh:
        for k,v in mesh.iteritems():
            parfile.set('MESH', k, v)

    if pressure:
        for k,v in pressure.iteritems():
            parfile.set('PRESSURE', k, v)

    if velocity:
        for k,v in velocity.iteritems():
            parfile.set('VELOCITY', k, v)

    if temperature:
        for k,v in temperature.iteritems():
            parfile.set('TEMPERATURE', k, v)

    if cvode:
        for k,v in cvode.iteritems():
            parfile.set('CVODE', k, v)

    with open(outfile, 'w') as f:
        parfile.write(f)


