import os, stat, re


def config_makenek(infile, outfile, opts):
    with open(infile, 'r') as f:
        lines = f.readlines()

    for key, val in opts.iteritems():
        lines = [re.sub(r'^#*{0}=\"+.+?\"+'.format(key), r'{0}="{1}"'.format(key, val), l) for l in lines]

    lines = [re.sub(r'(^source\s+\$SOURCE_ROOT/makenek.inc)', r'\g<1> >compiler.out', l)
             for l in lines]

    lines = [re.sub(r'(.+)2>&1\s+\|\s*tee\s+compiler.out', r'\g<1>', l)
             for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)
    os.chmod(outfile,
             stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH |
             stat.S_IRUSR | stat.S_IRGRP | stat.S_IROTH |
             stat.S_IWUSR)


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
             stat.S_IWUSR)


def config_basics_inc(infile, outfile, nelm):
    with open(infile, 'r') as f:
        lines = f.readlines()

    lines = [re.sub(r'(.*nelm *= *)[ 0-9]+(.*)', r'\g<1>{0}\g<2>'.format(nelm), l, flags=re.I)
             for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)


def config_size(infile, outfile, opts):
    with open(infile, 'r') as f:
        lines = f.readlines()

    # Substitute all the variables
    for key, value in opts.iteritems():
        if value:
            lines = [
                re.sub(
                    r'(.*\bparameter\b.*\b{0} *= *)\S+?( *[),])'.format(key),
                    r'\g<1>{0}\g<2>'.format(value), l, flags=re.I)
                for l in lines]

    with open(outfile, 'w') as f:
        f.writelines(lines)


def config_parfile(infile, outfile, opts):
    """ Set values in a parfile using ConfigParser

    Given a path to infile, substitute the options & values in kwargs, then
    output to outfile.  The infile and outfile can be the same file.

    opts is interpreted as a nested dict of the form:
        {section: {optname: value, ...}, ...}
    where "optname = value" are set in [section]. If 'optname' is not set in
    infile, then it will be added to outfile.  If 'optname' is already set in
    infile, then it will be overridden in outfile.  If an option is listed in
    infile but is not listed in in 'opts', then it will be copied to outfile
    without modification..

    Args:
        infile (str): Path to input file
        outfile (str): Path to output file
        opts ({section: {optname : value}}): Set each "optname = value" in
            each "[section]"
    """
    import ConfigParser

    parfile = ConfigParser.SafeConfigParser()
    parfile.read(infile)

    for section, name_vals in opts.iteritems():
        for name, val in name_vals.iteritems():
            parfile.set(section, name, val)

    with open(outfile, 'w') as f:
        parfile.write(f)
