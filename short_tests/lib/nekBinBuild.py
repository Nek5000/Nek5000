import os
import sys
from subprocess import call, check_call, Popen, PIPE, STDOUT

def build_tools(tools_root, tools_bin, f77=None, cc=None, bigmem=None,
                targets=('clean', 'all'), verbose=False):

    print('Compiling tools... ')
    print('    Using source directory "{0}"'.format(tools_root))
    print('    Using output directory "{0}"'.format(tools_bin))
    print('    Using FC "{0}"'.format(f77))
    print('    Using CC "{0}"'.format(cc))

    maketools_in  = os.path.join(tools_root, 'maketools')
    maketools_log = os.path.join(tools_root, 'maketools.out')

    my_env = os.environ.copy()
    if f77: my_env["FC"] = f77
    if cc:  my_env["CC"] = cc 
    my_env["bin_nek_tools"] = tools_bin 

    with open(maketools_log, 'w') as f:
        for t in targets:
            proc = Popen(
                [maketools_in, t],
                env=my_env,
                cwd=tools_root,
                stderr=STDOUT,
                stdout=PIPE
            )
            for line in proc.stdout:
                sys.stdout.write(line)
                f.write(line)

    proc.wait()

    if proc.returncode != 0:
        print('Could not compile tools! Check "{0}" for details.'.format(maketools_log))
        #exit(-1)
    else:
        print('Successfully compiled tools!')

def build_nek(source_root, usr_file, cwd=None, opts=None, verbose=False):

    if not opts:
        _opts = {}
    else:
        _opts = opts.copy()
    _opts.update(NEK_SOURCE_ROOT=source_root)

    print('Compiling nek5000...')
    print('    Using source directory "{0}"'.format(source_root))
    print('    Using working directory "{0}"'.format(cwd))
    print('    Using .usr file "{0}"'.format(usr_file))
    for key, val in _opts.iteritems():
        print('    Using {0}="{1}"'.format(key, val))

    my_env = os.environ.copy()
    if source_root      : my_env["NEK_SOURCE_ROOT"] = source_root
    if _opts.get('F77')   : my_env["FC"] = _opts.get('F77') 
    if _opts.get('CC')    : my_env["CC"] = _opts.get('CC')
    if _opts.get('PPLIST') : my_env["PPLIST"] = _opts.get('PPLIST') 

    makenek_in  = os.path.join(source_root, 'bin', 'makenek')
    logfile     = os.path.join(cwd, 'compiler.out')

    proc = Popen(
        [makenek_in, 'clean'], 
        cwd=cwd,
        env=my_env,
        stdin=PIPE)
    proc.communicate(input='Y\n')
    proc.wait()

    proc = Popen(
        [makenek_in, usr_file], 
        cwd=cwd,
        env=my_env,
        stdin=PIPE, 
        stderr=STDOUT, 
        stdout=PIPE)

    f = open(logfile, 'w')
    for line in proc.stdout:
        f.write(line)
        if verbose:
           sys.stdout.write(line)

    proc.wait()

    if proc.returncode != 0:
        print('Could not compile nek5000! Check "{0}" for details.'.format(logfile))
        exit(-1)
    else:
        print('Successfully compiled nek5000!')
