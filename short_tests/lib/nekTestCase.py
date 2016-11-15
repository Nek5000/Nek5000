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
        source_root (str):    Path to Nek source directory;overridden by $SOURCE_ROOT env variable
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
        self.f77            = 'mpif77'
        self.cc             = 'mpicc'
        self.ifmpi          = True
        self.verbose        = True
        self.source_root    = os.path.dirname(os.path.dirname(inspect.getabsfile(self.__class__)))
        self.examples_root  = os.path.dirname(inspect.getabsfile(self.__class__))
        self.tools_root     = ''
        self.tools_bin      = ''
        self.log_root       = ''
        self.makenek        = ''
        self.serial_procs   = 1
        self.parallel_procs = 4
        self.size_params    = {}
        self.cvode_dir      = ""

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

        # Get compilers from env, default to GNU
        # --------------------------------------
        self.f77     = os.environ.get('F77',   self.f77)
        self.cc      = os.environ.get('CC',    self.cc)
        self.ifmpi   = os.environ.get('IFMPI', self.ifmpi)
        self.verbose = os.environ.get('VERBOSE_TESTS', self.verbose)
        self.parallel_procs = int(os.environ.get('PARALLEL_PROCS', self.parallel_procs))

        # String/bool conversions
        self.ifmpi = str(self.ifmpi).lower()
        self.ifmpi = self.ifmpi == 'yes' or self.ifmpi == 'true'

        self.verbose = str(self.verbose).lower()
        self.verbose = self.verbose == 'yes' or self.verbose == 'true'

        for varname, varval in (
                ('F77', self.f77),
                ('CC', self.cc),
                ('IFMPI', str(self.ifmpi).lower()),
                ('VERBOSE_TESTS', str(self.verbose).lower()),
                ('PARALLEL_PROCS', self.parallel_procs)
        ):
            print('    Using {0}={1}'.format(varname, varval))

        # SOURCE_ROOT and EXAMPLES_ROOT must be defined.  Get from env and fail early if they don't exist
        # -----------------------------------------------------------------------------------------------
        self.source_root   = os.path.abspath(os.environ.get('SOURCE_ROOT',   self.source_root))
        self.examples_root = os.path.abspath(os.environ.get('EXAMPLES_ROOT', self.examples_root))

        for (varname, varval) in (('SOURCE_ROOT', self.source_root), ('EXAMPLES_ROOT', self.examples_root)):
            if os.path.isdir(varval):
                print('    Using {0} at {1}'.format(varname, varval))
            else:
                raise ValueError(
                    'The {0} directory "{1}" does not exist. Please provide a valid directory using the env variable {0} ROOT.'.format(varname, varval))

        # TOOLS_ROOT and TOOLS_BIN have default values, if not defined.  Raise error if they don't exist
        # ------------------------------------------------------------------------------------------------
        self.tools_root = os.environ.get('TOOLS_ROOT', self.tools_root)
        if self.tools_root:
            self.tools_root = os.path.abspath(self.tools_root)
        else:
            self.tools_root = os.path.abspath(os.path.join(self.source_root, 'tools'))

        if os.path.isdir(self.tools_root):
            print('    Using {0} at {1}'.format('TOOLS_ROOT', self.tools_root))
        else:
            raise ValueError(
                'The {0} directory "{1}" does not exist. Please provide a valid directory using the env variable {0} ROOT.'.format('TOOLS_ROOT', self.tools_root))

        # TOOLS_BIN has a default value.
        # -----------------------------
        self.tools_bin = os.environ.get('TOOLS_BIN', self.tools_bin)
        if self.tools_bin:
            self.tools_bin = os.path.abspath(self.tools_bin)
        else:
            self.tools_bin = os.path.abspath(os.path.join(self.source_root, 'bin'))

        # LOG_ROOT has no default value and can remain undefined
        # ------------------------------------------------------
        self.log_root = os.environ.get('LOG_ROOT', '')
        if self.log_root:
            self.log_root = os.path.abspath(self.log_root)

        # If TOOLS_BIN or LOG_ROOT don't exist, make them
        #------------------------------------------------
        for varval, varname in ((self.tools_bin, 'TOOLS_BIN'), (self.log_root,  'LOG_ROOT')):
            if varval:
                if os.path.isdir(varval):
                    print('    Using {0} at {1}'.format(varname, varval))
                else:
                    print('    The {0} directory, "{1}" does not exist.  It will be created'.format(varname, varval))
                    os.makedirs(varval)

        # CVODE_DIR doesn't need to be defined.  It defaults to ""
        #---------------------------------------------------------
        self.cvode_dir = os.environ.get('CVODE_DIR', self.cvode_dir)

        # Default destination of makenek
        # ------------------------------
        if not self.makenek:
            self.makenek   = os.path.join(self.source_root, 'core', 'makenek')

        print("Finished getting setup options!")

    def build_tools(self, targets=None, tools_root=None, tools_bin=None, f77=None, cc=None, bigmem=None, verbose=None):
        from lib.nekBinBuild import build_tools
        build_tools(
            targets    = targets    if targets    else ('clean', 'genmap'),
            tools_root = tools_root if tools_root else self.tools_root,
            tools_bin  = tools_bin  if tools_bin  else self.tools_bin,
            f77        = f77        if f77        else 'gfortran',
            cc         = cc         if cc         else 'gcc',
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
            F77   = self.f77,
            CC    = self.cc,
            IFMPI = str(self.ifmpi).lower(),
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
