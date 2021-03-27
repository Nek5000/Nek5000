import os
import stat
import re


def config_size(params, infile, outfile):
    """Redefine parameters in a SIZE file

    Given a path to infile, redefine the variables & values in params, then
    output to outfile.  The infile and outfile can be the same file.

    Only the parameters already declared in infile may be set.  This function
    doesn't allow you to declare any new parameters.
    (TODO: Raise warning when attempting to declare a new parameter.)

    Args:
        params ({variable : value}): Each 'variable' will be set to 'value'
            in the output SIZE file
        infile (str): Path to the input SIZE file
        outfile (str): Path to output SIZE file

    """
    with open(infile, "r") as file:
        lines = file.readlines()

    # Substitute all the variables
    for key, value in list(params.items()):
        if value:
            lines = [
                re.sub(
                    fr"(.*\bparameter\b.*\b{key} *= *)\S+?( *[),])",
                    fr"\g<1>{value}\g<2>",
                    l,
                    flags=re.I,
                )
                for l in lines
            ]

    with open(outfile, "w") as file:
        file.writelines(lines)


def config_parfile(opts, infile, outfile):
    """Set values in a parfile using ConfigParser

    Given a path to infile, substitute the options & values in opts, then
    output to outfile.  The infile and outfile can be the same file.

    opts is interpreted as a nested dict of the form:
        {section: {optname: value, ...}, ...}
    where "optname = value" are set in [section]. If 'optname' is not set in
    infile, then it will be added to outfile.  If 'optname' is already set in
    infile, then it will be overridden in outfile.  If an option is listed in
    infile but is not listed in in 'opts', then it will be copied to outfile
    without modification..

    Args:
        opts ({section: {optname : value}}): Set each "optname = value" in
            each "[section]"
        infile (str): Path to input parfile
        outfile (str): Path to output parfile

    """
    import configparser

    parfile = configparser.SafeConfigParser()
    parfile.read(infile)

    for section, name_vals in list(opts.items()):
        for name, val in list(name_vals.items()):
            parfile.set(section, name, val)

    with open(outfile, "w") as file:
        parfile.write(file)
