import unittest
import inspect
import os
from functools import wraps

###############################################################################
#  DECORATORS
###############################################################################

def pn_pn_serial(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self.mpi_procs = self.serial_procs
        self.log_suffix = '.pn_pn'
        if self.ifmpi:
            self.log_suffix += '.parallel'
        else:
            self.log_suffix += '.serial'
        method(self, *args, **kwargs)
    return wrapper

def pn_pn_2_serial(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self.mpi_procs = self.serial_procs
        self.log_suffix = '.pn_pn_2'
        if self.ifmpi:
            self.log_suffix += '.parallel'
        else:
            self.log_suffix= '.serial'
        method(self, *args, **kwargs)
    return wrapper

def pn_pn_parallel(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self.mpi_procs = self.parallel_procs
        if not self.ifmpi:
            self.skipTest("Skipping \"{0}\"; MPI is not enabled.".format(self.id()))
        else:
            self.log_suffix = '.pn_pn'
            if self.ifmpi:
                self.log_suffix += '.parallel'
            else:
                self.log_suffix += '.serial'
            method(self, *args, **kwargs)
    return wrapper

def pn_pn_2_parallel(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        self.mpi_procs = self.parallel_procs
        if not self.ifmpi:
            self.skipTest("Skipping \"{0}\"; MPI is not enabled.".format(self.id()))
        else:
            # Set number of mpi procs
            self.log_suffix = '.pn_pn_2'
            if self.ifmpi:
                self.log_suffix += '.parallel'
            else:
                self.log_suffix += '.serial'
            method(self, *args, **kwargs)
    return wrapper

###############################################################################
#  BASE TEST CASE
###############################################################################

class NekTestCase(unittest.TestCase):
    """ Base class for Nek unittests

    This defines a setUpClass method to:
        (a) get the relevant environment variables for compilers, directories
        (b) add env vars to maketools, makenek
        (b) build tools
    All subclassed TestCases will need to do these things.

    Class attributes:
        f77 (str):            The Fortran 77 compiler to use     [default: 'gfortran']
        cc (str):             The C compiler to use              [default: 'gcc']
        ifmpi (bool):         Perform compilation/tests with MPI [default: False]
        source_root (str):    Path to Nek source directory;overridden by $NEK_SOURCE_ROOT env variable
                              [default: '$HOME/nek5_svn/trunk/nek']
        tools_root (str):     Path to Nek tools directory; overridden by $TOOLS_ROOT env variable
                              [default: '$HOME/nek5_svn/trunk/tools']
        examples_root (str):  Path to Nek examples directory; overridden by $EXAMPLES_ROOT env variable
                              [default: '$HOME/nek5_svn/examples']
        makenek (str):        Path to makenek                    [default: source_root/makenek]
        tools_bin (str):      Directory to place compiled tools  [default: tools_root/bin]

    Subclass attributes:
        These aren't meaningful in the base class.  They're intended for a subclass that represents
        a particular example problem.
        example_subdir (str): The subdirectory for the subclass' example.  Assumed that it's in example_root
        rea_file (str):       The .rea file for the subclass' example, minus the '.rea' extension.  Assumed
                              that it's in example_root/example_dir
        size_file (str):      The SIZE file for the subclass' example.  Assuemed that it's in
                              example_root/example_subdir
    """
    # Defined in subclasses only; declared here to make syntax checker happy
    example_subdir      = ""
    case_name           = ""

    def __init__(self, *args, **kwargs):
        # These can be overridden by self.get_opts
        self.f77            = ""
        self.cc             = ""
        self.pplist         = ""
        self.usr_lflags     = ""
        self.ifmpi          = True

        self.source_root    = os.path.dirname(os.path.dirname(inspect.getabsfile(self.__class__)))
        self.examples_root  = os.path.dirname(inspect.getabsfile(self.__class__))
        self.makenek        = os.path.join(self.source_root, 'bin', 'makenek')
        self.tools_root     = os.path.join(self.source_root, 'tools')
        self.tools_bin      = os.path.join(self.source_root, 'bin')
        self.log_root       = ""
        self.verbose        = True
        self.serial_procs   = 1
        self.parallel_procs = 2
        self.size_params    = {}

        # These are overridden by method decorators (pn_pn_serial, pn_pn_parallel,
        # pn_pn_2_serial, and pn_pn_2_parallel)
        self.log_suffix = ""
        self.mpi_procs  = None

        # Empy list of delayed fails
        self._delayed_failures = []

        self.get_opts()

        unittest.TestCase.__init__(self, *args, **kwargs)


    def assertAlmostEqualDelayed(self, test_val, target_val, delta, label):
        if abs(test_val-target_val) <= delta:
            msg = '    SUCCESS: {0}: Test value {1} equals target value {2} +/- {3}'.format(label, test_val, target_val, delta)
        else:
            msg = '    FAILURE: {0}: Test value {1} exceeds target value {2} +/- {3}'.format(label, test_val, target_val, delta)
            self._delayed_failures.append(msg)
        print(msg)


    def assertIsNotNullDelayed(self, test_val, label):
        if test_val:
            msg = 'SUCCESS: Found phrase "{0}" in logfile.'.format(label)
        else:
            msg = 'FAILURE: Unexpectedly did not find phrase "{0}" in logfile'.format(label)
            self._delayed_failures.append(msg)
        print(msg)

    def assertIsNullDelayed(self, test_val, label):
        if test_val:
            msg = 'FAILURE: Found phrase "{0}" in logfile.'.format(label)
            self._delayed_failures.append(msg)
        else:
            msg = 'SUCCESS: Did not find phrase "{0}" in logfile'.format(label)
        print(msg)

    def assertDelayedFailures(self):
        if self._delayed_failures:
            report = [
                '\n\nFailed assertions:{0}\n'.format(len(self._delayed_failures))
            ]
            for i,failure in enumerate(self._delayed_failures, start=1):
                report.append('{0}: {1}'.format(i, failure))
            #self._delayed_failures = []
            self.fail('\n'.join(report))


    def get_opts(self):

        print("Getting setup options...")

        # Get compiler options from env
        self.f77            = os.environ.get('FC', self.f77)
        self.cc             = os.environ.get('CC', self.cc)
        self.pplist         = os.environ.get('PPLIST', self.pplist)
        self.usr_lflags     = os.environ.get('USR_LFLAGS', self.usr_lflags)
        self.ifmpi          = os.environ.get('MPI', self.ifmpi)

        # Get paths from env
        try:
            self.source_root = os.path.abspath(os.environ['NEK_SOURCE_ROOT'])
        except KeyError:
            pass
        else:
            self.makenek        = os.path.join(self.source_root, 'bin', 'makenek')
            self.tools_root     = os.path.join(self.source_root, 'tools')
            self.tools_bin      = os.path.join(self.source_root, 'bin')

        self.examples_root = os.path.abspath(os.environ.get('EXAMPLES_ROOT', self.examples_root))
        self.tools_root    = os.path.abspath(os.environ.get('TOOLS_ROOT', self.tools_root))
        self.tools_bin     = os.path.abspath(os.environ.get('TOOLS_BIN', self.tools_bin))

        try:
            self.log_root = os.path.abspath(os.environ['LOG_ROOT'])
        except KeyError:
            pass

        self.verbose        = str(os.environ.get('VERBOSE_TESTS', self.verbose)).lower() == 'true'
        self.parallel_procs = int(os.environ.get('PARALLEL_PROCS', self.parallel_procs))

        # Print everything out
        for varname, varval in (
                ('FC', self.f77),
                ('CC', self.cc),
                ('PPLIST', self.pplist),
                ('USR_LFLAGS', self.usr_lflags),
                ('IFMPI', self.ifmpi),
                ('NEK_SOURCE_ROOT', self.source_root),
                ('EXAMPLES_ROOT', self.examples_root),
                ('LOG_ROOT', self.log_root),
                ('TOOLS_ROOT', self.tools_root),
                ('TOOLS_BIN', self.tools_bin),
                ('VERBOSE_TESTS', self.verbose),
                ('PARALLEL_PROCS', self.parallel_procs)
        ):
            if varval:
                print('    Using {0:14} = "{1}"'.format(varname, varval))

        # Verify that pathnames are valid
        for varname, varval in (
                ('NEK_SOURCE_ROOT', self.source_root),
                ('EXAMPLES_ROOT', self.examples_root),
                ('LOG_ROOT', self.log_root),
                ('TOOLS_ROOT', self.tools_root),
                ('TOOLS_BIN', self.tools_bin),
        ):
            if varval and not os.path.isdir(varval):
                raise OSError('The {0} directory "{1}" does not exist. Please the env variable ${0} to a valid directory.'.format(varname, varval))

        print("Finished getting setup options!")

    def build_tools(self, targets=None, tools_root=None, tools_bin=None, f77=None, cc=None, bigmem=None, verbose=None):
        from lib.nekBinBuild import build_tools
        build_tools(
            targets    = targets    if targets    else ('clean', 'genmap'),
            tools_root = tools_root if tools_root else self.tools_root,
            tools_bin  = tools_bin  if tools_bin  else self.tools_bin,
            f77        = f77        ,
            cc         = cc         ,
            bigmem     = bigmem     if bigmem     else 'false',
            verbose    = verbose    if verbose    else self.verbose
        )

    def config_size(self, params=None, infile=None, outfile=None):
        from lib.nekFileConfig import config_size
        cls = self.__class__

        if not infile:
            infile = os.path.join(self.source_root, 'core', 'SIZE.template')
        if not outfile:
            outfile = os.path.join(self.examples_root, cls.example_subdir, 'SIZE')
        if not params:
            params = self.size_params

        config_size(params=params, infile=infile, outfile=outfile)

    def mkSIZE(self, case=None):
        cls = self.__class__

        if not case:
            case = cls.case_name

        workdir = os.path.join(self.examples_root, cls.example_subdir)
        os.system('cd ' + workdir + ' ; ' + self.source_root + '/bin/mkSIZE ' + case)

    def config_parfile(self, opts=None, infile=None, outfile=None):
        from lib.nekFileConfig import config_parfile
        cls = self.__class__

        if not infile:
            infile = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.par')
        if not outfile:
            outfile = infile
        if not opts:
            opts = {}

        config_parfile(opts=opts, infile=infile, outfile=outfile)

    def run_genmap(self, rea_file=None, tol='0.5'):

        from lib.nekBinRun import run_meshgen
        cls = self.__class__

        if not rea_file:
            rea_file = cls.case_name

        run_meshgen(
            command = os.path.join(self.tools_bin, 'genmap'),
            stdin   = [rea_file, tol],
            cwd     = os.path.join(self.examples_root, cls.example_subdir),
            verbose = self.verbose
        )

    def run_genbox(self, box_file=None):
        from lib.nekBinRun import run_meshgen
        if not box_file:
            box_file = self.__class__.case_name

        # Fix extension, in case user doesn't provide it
        root, ext = os.path.splitext(box_file)
        if ext != '.box':
            box_file = root + ext + '.box'

        run_meshgen(
            command = os.path.join(self.tools_bin, 'genbox'),
            stdin   = [box_file],
            cwd     = os.path.join(self.examples_root, self.__class__.example_subdir),
            verbose = self.verbose
        )

    def run_n2to3(self, stdin):
        from lib.nekBinRun import run_meshgen
        run_meshgen(
            command = os.path.join(self.tools_bin, 'n2to3'),
            stdin   = stdin,
            cwd     = os.path.join(self.examples_root, self.__class__.example_subdir),
        )

    def build_nek(self, opts=None, usr_file=None):
        from lib.nekBinBuild import build_nek
        cls = self.__class__

        if not usr_file:
            usr_file = cls.case_name

        all_opts = dict(
            FC = self.f77,
            CC = self.cc,
            PPLIST = self.pplist,
            USR_LFLAGS = self.usr_lflags,
            MPI = int(self.ifmpi),
        )
        if opts:
            all_opts.update(opts)

        build_nek(
            source_root = self.source_root,
            usr_file    = usr_file,
            cwd         = os.path.join(self.examples_root, cls.example_subdir),
            opts        = all_opts,
            verbose     = self.verbose,
        )

    def run_nek(self, rea_file=None, step_limit=None):
        from lib.nekBinRun import run_nek
        cls = self.__class__
        run_nek(
            cwd        = os.path.join(self.examples_root, cls.example_subdir),
            rea_file   = cls.case_name if not rea_file else rea_file,
            ifmpi      = self.ifmpi,
            log_suffix = self.log_suffix,
            n_procs    = self.mpi_procs,
            step_limit = step_limit,
            verbose    = self.verbose
        )

    def move_logs(self):
        cls = self.__class__
        if self.log_root:
            src_dir = os.path.join(self.examples_root, cls.example_subdir)
            dest_dir = os.path.join(self.log_root, cls.example_subdir)
            print("Moving logs...")
            if not os.path.isdir(dest_dir):
                print("    Making subdirectory {0}")
                os.makedirs(os.path.join(dest_dir))

            for f in os.listdir(src_dir):
                if f == 'compiler.out' or f == 'genmap.out' or 'log' in f:
                    src_file = os.path.join(src_dir, f)
                    dest_file = os.path.join(dest_dir, f)
                    try:
                        os.rename(src_file, dest_file)
                    except OSError as E:
                        # TODO: change to warnings.warning
                        print("    Could not move {0} to {1}: {2}".format(src_file, dest_file, E))
                    else:
                        print("    Moved {0} to {1}".format(src_file, dest_file))

    def mvn(self, src_prefix, dest_prefix):
        from lib.nekBinRun import mvn
        cls = self.__class__
        mvn(src_prefix, dest_prefix,
            cwd = os.path.join(self.examples_root, cls.example_subdir)
        )


    def get_value_from_log(self, label, column, row=0, logfile=None):
        cls = self.__class__
        if not logfile:
            logfile = os.path.join(
                self.examples_root,
                cls.example_subdir,
                '{0}.log.{1}{2}'.format(cls.case_name, self.mpi_procs, self.log_suffix)
            )
        # Get all lines with label
        with open(logfile, 'r') as f:
            line_list = [l for l in f if label in l]
        if not line_list:
            raise ValueError("Could not find label \"{0}\" in logfile \"{1}\".  The run may have failed.".format(label, logfile))
        try:
            value = float(line_list[row].split()[column])
        except ValueError:
            raise ValueError("Attempted to parse non-numerical value in logfile, \"{0}\".  Logfile may be malformatted or run may have failed".format(logfile))
        except IndexError:
            raise IndexError("Fewer rows and/or columns than expected in logfile, \"{0}\".  Logfile may be malformmated or run may have failed.".format(logfile))
        else:
            return value

    def get_phrase_from_log(self, label, logfile=None, row=0):
        cls = self.__class__
        if not logfile:
            logfile = os.path.join(
                self.examples_root,
                cls.example_subdir,
                '{0}.log.{1}{2}'.format(cls.case_name, self.mpi_procs, self.log_suffix)
            )

        with open(logfile, 'r') as f:
            line_list = [l for l in f if label in l]

        try:
            line = line_list[row]
        except IndexError:
            return None
        else:
            return line
