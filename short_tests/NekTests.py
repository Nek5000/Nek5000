#!/usr/bin/env python
from lib.nekTestCase import *
from unittest import skip

import re

###############################################################################

class FsHydro(NekTestCase):
    example_subdir = 'fs_hydro'
    case_name       = 'fs_hydro'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '12',
            lxd       = '18',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '1',
            lhis      = '100',
            lelx      = '20',
            lely      = '60',
            lelz      = '1',
            lx1m      = 'lx1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=90., label='gmres')

        axhm = self.get_value_from_log('axhm', column=-3,)
        self.assertAlmostEqualDelayed(axhm, target_val=0., delta=19180., label='axhm')

        amp = self.get_value_from_log('AMP', column=-2, row=-1)
        self.assertAlmostEqualDelayed(amp, target_val=-5.2985368e-05, delta=4e-05, label='AMP')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

###############################################################################

class Axi(NekTestCase):
    example_subdir  = 'axi'
    case_name        = 'axi'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        errmx = self.get_value_from_log('err ', column=-2, row=-1)
        self.assertAlmostEqualDelayed(errmx, target_val=2.639593E-15, delta=1E-15, label='L2 err')

        errl2 = self.get_value_from_log('err ', column=-3, row=-1)
        self.assertAlmostEqualDelayed(errl2, target_val=1.043610E-14, delta=1E-15, label='MAX err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        errmx = self.get_value_from_log('err ', column=-2, row=-1)
        self.assertAlmostEqualDelayed(errmx, target_val=3.963589E-15, delta=1E-15, label='L2 err')

        errl2 = self.get_value_from_log('err ', column=-3, row=-1)
        self.assertAlmostEqualDelayed(errl2, target_val=1.687539E-14, delta=1E-15, label='MAX err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################

class Benard_Ray9(NekTestCase):
    example_subdir = 'benard'
    case_name = 'ray_9'

    def setUp(self):
        self.size_params = dict (
            ldim      = '2',
            lx1       = '12',
            lxd       = '16',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])
        self.run_genmap(rea_file='ray_9')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9')

        gmres = self.get_value_from_log('gmres ', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=17., label='gmres')

        rac = self.get_value_from_log('converged_rac', column=-2)
        self.assertAlmostEqualDelayed(rac, target_val=1707.79, delta=0.01, label='rac')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9')

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=11., label='gmres')

        rac = self.get_value_from_log('converged_rac', column=-2)
        self.assertAlmostEqualDelayed(rac, target_val=1707.79, delta=0.01, label='rac')

    def tearDown(self):
        self.move_logs()

#####################################################################

class Eddy_EddyUv(NekTestCase):
    example_subdir  = 'eddy'
    case_name        = 'eddy_uv'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])

        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*DIVERGENCE$', '      0.10000E-08', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        # Update SIZE parameters for PnPn
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=12., label='gmres')

        runtime = self.get_value_from_log('total solver time w/o IO ', column=-2,)
        self.assertAlmostEqualDelayed(runtime, target_val=0., delta=60., label='runtime')

        crsl = self.get_value_from_log('crsl ', column=-3,)
        self.assertAlmostEqualDelayed(crsl, target_val=0., delta=1500., label='crsl')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=1.448024E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=10., label='gmres')

        crsl = self.get_value_from_log('crsl ', column=-3,)
        self.assertAlmostEqualDelayed(crsl, target_val=0., delta=1500., label='crsl')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=1.448024E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=8., label='gmres')

        crsl = self.get_value_from_log('crsl ', column=-3,)
        self.assertAlmostEqualDelayed(crsl, target_val=0., delta=1050., label='crsl')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.759103E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=7.842019E-05, delta=1E-06, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=6.896211E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################

class Eddy_LegacySize(NekTestCase):
    example_subdir  = 'eddy'
    case_name       = 'eddy_uv'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            lx2       = '',
            ly2       = '',
            lz2       = '',
        )

        self.build_tools(['genmap'])

        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*DIVERGENCE$', '      0.10000E-08', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.size_params['ly2'] = 'ly1'
        self.size_params['lz2'] = 'lz1'
        self.config_size(infile=os.path.join(self.examples_root, self.__class__.example_subdir, 'SIZE.legacy'))
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=34., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=1.448024E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.size_params['ly2'] = 'ly1-2'
        self.size_params['lz2'] = 'lz1'
        self.config_size(infile=os.path.join(self.examples_root, self.__class__.example_subdir, 'SIZE.legacy'))
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=22., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.759103E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=7.842019E-05, delta=1E-06, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=6.896211E-05, delta=1E-06, label='P err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

######################################################################

class Eddy_Rich(NekTestCase):
    example_subdir  = 'eddy_rich'
    case_name        = 'eddy_rich'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])

        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*DIVERGENCE$', '      0.10000E-08', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=5.497982E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=8.064398E-05, delta=1E-06, label='Y err')

        perr = self.get_value_from_log('P err', column=-5, row=-1)
        self.assertAlmostEqualDelayed(perr, target_val=2.272926E-04, delta=1E-04, label='P err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################

class Eddy_Neknek(NekTestCase):
    example_subdir  = 'eddy_neknek'
    case_name       = 'eddy_uv'

    def setUp(self):

        self.size_params = dict(
            ldim='2',
            lx1='8',
            lxd='12',
            lx2='lx1-2',
            lelg='1000',
            lpert='1',
            nsessmax='2',
        )

        self.build_tools(['genmap'])

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        from lib.nekBinRun import run_neknek
        from re import sub

        cls = self.__class__
        cwd = os.path.join(self.examples_root, cls.example_subdir)

        # Tweak the .rea files and run genmap
        for rea_file in ('inside', 'outside'):
            rea_path = os.path.join(cwd, rea_file + '.rea')
            with open(rea_path, 'r') as f:
                lines = [sub(r'^.*DIVERGENCE$', '      1.0000000E-06     p21 DIVERGENCE', l) for l in f]
            with open(rea_path, 'w') as f:
                f.writelines(lines)
            self.run_genmap(os.path.join(cwd, rea_file),tol='0.2')

        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek(opts={'PPLIST':'NEKNEK'})
        run_neknek(
            cwd = cwd,
            inside = 'inside',
            outside = 'outside',
            np_inside = 1,
            np_outside = 1,
            step_limit = 1000,
            log_suffix = self.log_suffix,
            verbose = self.verbose,
        )

        logfile  = os.path.join(cwd, '{inside}{np_in}.{outside}{np_out}.log{sfx}'.format(
            inside = 'inside',
            outside = 'outside',
            np_in = 1,
            np_out = 1,
            sfx = self.log_suffix
        ))

        xerr_inside = self.get_value_from_log('X err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_inside, target_val=7.163001E-04, delta=1E-05, label='X err  inside')

        xerr_global = self.get_value_from_log('X err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_global, target_val=8.580050E-04, delta=1E-05, label='X err   global')

        xerr_outside = self.get_value_from_log('X err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_outside, target_val=8.580050E-04, delta=1E-05, label='X err  outside')

        yerr_inside = self.get_value_from_log('Y err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_inside, target_val=9.012947E-04, delta=1E-05, label='Y err  inside')

        yerr_global = self.get_value_from_log('Y err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_global, target_val=9.877146E-04, delta=1E-05, label='Y err   global')

        yerr_outside = self.get_value_from_log('Y err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_outside, target_val=9.877146E-04, delta=1E-05, label='Y err  outside')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        from lib.nekBinRun import run_neknek
        from re import sub

        cls = self.__class__
        cwd = os.path.join(self.examples_root, cls.example_subdir)

        # Tweak the .rea files and run genmap
        for rea_file in ('inside', 'outside'):
            rea_path = os.path.join(cwd, rea_file + '.rea')
            with open(rea_path, 'r') as f:
                lines = [sub(r'^.*DIVERGENCE$', '      1.0000000E-11     p21 DIVERGENCE', l) for l in f]
            with open(rea_path, 'w') as f:
                f.writelines(lines)
            self.run_genmap(os.path.join(cwd, rea_file),tol='0.2')

        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(opts={'PPLIST':'NEKNEK'})
        run_neknek(
            cwd = cwd,
            inside = 'inside',
            outside = 'outside',
            np_inside = 1,
            np_outside = 1,
            step_limit = 1000,
            log_suffix = self.log_suffix,
            verbose = self.verbose,
        )

        logfile  = os.path.join(cwd, '{inside}{np_in}.{outside}{np_out}.log{sfx}'.format(
            inside = 'inside',
            outside = 'outside',
            np_in = 1,
            np_out = 1,
            sfx = self.log_suffix
        ))

        xerr_inside = self.get_value_from_log('X err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_inside, target_val=7.431657E-04, delta=1E-04, label='X err  inside')

        xerr_global = self.get_value_from_log('X err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_global, target_val=8.696332E-04, delta=1E-04, label='X err   global')

        xerr_outside = self.get_value_from_log('X err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(xerr_outside, target_val=8.696332E-04, delta=1E-05, label='X err  outside')

        yerr_inside = self.get_value_from_log('Y err  inside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_inside, target_val=9.250194E-04, delta=1E-04, label='Y err  inside')

        yerr_global = self.get_value_from_log('Y err   global', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_global, target_val=9.878329E-04, delta=1E-04, label='Y err   global')

        yerr_outside = self.get_value_from_log('Y err  outside', logfile=logfile, column=-7, row=-1)
        self.assertAlmostEqualDelayed(yerr_outside, target_val=9.878329E-04, delta=1E-05, label='Y err  outside')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class KovStokes(NekTestCase):
    # Note: Legacy Analysis.py script only checked Pn-Pn-2 test cases
    example_subdir = 'kov_stokes'
    case_name = 'kov_stokes'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '14',
            lxd       = '20',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        err = self.get_value_from_log(label='err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err, target_val= 3.93249E-08, delta=6e-08, label='err')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1-0'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        err = self.get_value_from_log(label='err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err, target_val=5.05960E-13 , delta=1e-14, label='err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class Ethier(NekTestCase):
    example_subdir = 'ethier'
    case_name = 'ethier'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '50',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=1000)

        herr = self.get_value_from_log(label='hpts err', column=-1, row=-1)
        self.assertAlmostEqualDelayed(herr, target_val=1.3776e-08, delta=1e-08, label='hpts err')

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=7., label='gmres')

        vxerr = self.get_value_from_log(label='L2 err', column=-4, row=-1)
        self.assertAlmostEqualDelayed(vxerr, target_val=3.635317e-05, delta=1e-07, label='VX err')

        prerr = self.get_value_from_log(label='L2 err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(prerr, target_val=1.127384e-04, delta=1e-06, label='PR err')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()

        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.par')
        with open(rea_path, 'r') as f:
            lines = [sub(r'.*preconditioner.*', 'preconditioner = semg_amg', l, flags=re.I) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)

        self.run_nek(step_limit=1000)

        herr = self.get_value_from_log(label='hpts err', column=-1, row=-1)
        self.assertAlmostEqualDelayed(herr, target_val=1.3776e-08, delta=1e-08, label='hpts err')

        gmres = self.get_value_from_log('gmres ', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=14., label='gmres')

        vxerr = self.get_value_from_log(label='L2 err', column=-4, row=-1)
        self.assertAlmostEqualDelayed(vxerr, target_val=2.407549E-006, delta=1e-08, label='VX err')

        prerr = self.get_value_from_log(label='L2 err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(prerr, target_val=7.554325E-005, delta=1e-08, label='PR err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class LowMachTest(NekTestCase):
    example_subdir = 'lowMach_test'
    case_name       = 'lowMach_test'

    def setUp(self):
        self.size_params = dict (
            ldim      = '2',
            lx1       = '14',
            lxd       = '20',
            lx2       = 'lx1-0',
            lx1m      = 'lx1',
            lelg      = '500',
        )
        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=35, label='gmres')

        vx = self.get_value_from_log(label='ERROR VX', column=-5, row=-1)
        self.assertAlmostEqualDelayed(vx, target_val=2.6938e-09, delta=1e-10, label='VX')

        errt = self.get_value_from_log(label='ERROR T', column=-5, row=-1)
        self.assertAlmostEqualDelayed(errt, target_val=4.5532e-12, delta=1e-13, label='T')

        qtl = self.get_value_from_log(label='ERROR QTL', column=-5, row=-1)
        self.assertAlmostEqualDelayed(qtl, target_val=2.6557E-06, delta=1e-07, label='QTL')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class MvCylCvode(NekTestCase):
    example_subdir = 'mv_cyl'
    case_name = 'mv_cyl'

    def setUp(self):
        self.size_params = dict (
            ldim     = '2',
            lx1      = '8',
            lxd      = '12',
            lx2      = 'lx1-0',
            lx1m     = 'lx1',
            lelg     = '500',
            ldimt    = '10',
            lcvelt   = 'lelt',
        )
        self.config_size()
        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.log_suffix += '.steps_1e3'
        self.config_parfile({'GENERAL' : {'numSteps' : '1e3', 'dt' : '1e-3'}})
        self.build_nek(opts={'PPLIST':'CVODE'})
        self.run_nek()

        err3 = self.get_value_from_log('err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err3, target_val=4.128154e-05, delta=1e-6, label='err p0th')

        err2 = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=0.6348537E-06, delta=1e-7, label='err dp/dt')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################

class ConjHt(NekTestCase):
    example_subdir  = 'conj_ht'
    case_name        = 'conj_ht'

    def setUp(self):
        self.build_tools(['genmap'])
        self.run_genmap()
        self.size_params = dict (
            ldim     = '2',
            lx1      = '4',
            lxd      = '7',
            lx2      = 'lx1-0',
            lelg     = '100',
            ldimt    = '2',
            lcvelt   = 'lelt',
        )

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek(opts={'PPLIST':'CVODE'})
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=20., label='gmres')

        axhm = self.get_value_from_log('axhm', column=-3,)
        self.assertAlmostEqualDelayed(axhm, target_val=0., delta=25800., label='axhm')

        tmax = self.get_value_from_log('tmax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(tmax, target_val=13.109, delta=1E-03, label='tmax')

        terr = self.get_value_from_log('tmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(terr, target_val=3.88933e-06, delta=5E-07, label='terr')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(opts={'PPLIST':'CVODE'})
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=16., label='gmres')

        axhm = self.get_value_from_log('axhm', column=-3,)
        self.assertAlmostEqualDelayed(axhm, target_val=0., delta=30600.0, label='axhm')

        tmax = self.get_value_from_log('tmax', column=-3, row=-1)
        self.assertAlmostEqualDelayed(tmax, target_val=13.1208, delta=1E-03, label='tmax')

        terr = self.get_value_from_log('tmax', column=-2, row=-1)
        self.assertAlmostEqualDelayed(terr, target_val=8.53943e-06, delta=5E-07, label='terr')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class CmtInviscidVortex(NekTestCase):
    example_subdir = os.path.join('CMT', 'inviscid_vortex')
    case_name = 'pvort'

    def diff_l2norms(self):
        def get_line(filename, line_num=0):
            with open(filename) as f:
                line = f.readlines()[line_num]
            return [float(x) for x in line.split()[1:]]

        cls = self.__class__
        test_vals = get_line(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat'))
        ref_vals = get_line(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat.ref'))
        for t, r in zip(test_vals, ref_vals):
            self.assertAlmostEqual(t, r, delta=0.1*r,
                msg='FAILURE: Last line of l2norms.dat differed from reference values by > 10%\n  test vals:{0}\n  ref vals: {1}'.format(test_vals, ref_vals))
        print('SUCCESS: Last line of l2norms.dat was within 10% of reference values\n  test vals:{0}\n  ref vals: {1}'.format(test_vals, ref_vals))

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '25',
            lxd       = '36',
            lx2       = 'lx1-0',
            lelg      = '50',
            ldimt     = '3',
            toteq     = '5',
        )
        self.config_size()
        self.build_tools(['genmap'])
        self.run_genmap()

        cls = self.__class__
        try:
            os.remove(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat'))
        except OSError:
            pass

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.build_nek(opts={'PPLIST':'CMTNEK'})
        self.run_nek(step_limit=1000)
        self.diff_l2norms()

    def tearDown(self):
        self.move_logs()
        
####################################################################

class LinCav_Dir(NekTestCase):
    example_subdir = 'dfh_cav'
    case_name = 'lin_dfh_cav_dir'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpelt     = 'lelt',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_dfh_cav')
        self.run_nek(step_limit=None)  

        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-7.57304E-03, delta=1E-06, label='growth rate')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

class LinCav_Adj(NekTestCase):
    example_subdir = 'dfh_cav'
    case_name = 'lin_dfh_cav_adj'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '9',
            lxd       = '13',
            lx2       = 'lx1-2',
            lelg      = '500',
            lpelt     = 'lelt',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='lin_dfh_cav')
        self.run_nek(step_limit=None)  

        omega = self.get_value_from_log('Energy', column=-3, row=-1)
        self.assertAlmostEqualDelayed(omega, target_val=-7.57304E-03, delta=1E-06, label='growth rate')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

#####################################################################

class DoubleShear(NekTestCase):
    example_subdir  = 'double_shear'
    case_name        = 'double_shear'

    def setUp(self):

        # Default SIZE parameters. Can be overridden in test cases
        self.size_params = dict(
            ldim      = '2',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        en = self.get_value_from_log('L2 norm: ', column=-1)
        self.assertAlmostEqualDelayed(en, target_val=0.86147612200488166, delta=1e-5, label='Energy')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        en = self.get_value_from_log('L2 norm: ', column=-1)
        self.assertAlmostEqualDelayed(en, target_val=0.85973996057622892, delta=1e-5, label='Energy')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class IO_Test(NekTestCase):
    example_subdir = 'io_test'
    case_name = 'io_test'

    def setUp(self):
        self.size_params = dict(
            ldim      = '3',
            lx1       = '6',
            lxd       = '9',
            lx2       = 'lx1-2',
            lelg      = '100',
            ldimt     = '3',
        )

        self.build_tools(['genmap'])
        self.run_genmap()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek(usr_file='io_test')
        self.run_nek(step_limit=None)  
        
        phrase = self.get_phrase_from_log('All I/O tests PASSED')
        self.assertIsNotNullDelayed(phrase, label='All I/O tests PASSED')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################

class InclDef(NekTestCase):
    example_subdir = 'incl_def'
    case_name = 'incl_def'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1',
            lelg      = '100',
            ldimt     = '2',
            lhis      = '100',
            lpert     = '1',
            toteq     = '1',
            lelx      = '1',
            lely      = '1',
            lelz      = '1',
            lx1m      = '1',
            lbelt     = '1',
            lpelt     = '1',
            lcvelt    = '1',
        )
        self.build_tools(['genmap'])
        self.run_genmap(tol='0.01')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        phrase = self.get_phrase_from_log('All include files added with success')
        self.assertIsNotNullDelayed(phrase, label='All include files added with success')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

###############################################################

if __name__ == '__main__':
    import unittest, argparse, os

    # Get arguments from command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--f77", default='mpif77', help="The Fortran 77 compiler to use [default: mpif77]")
    parser.add_argument("--cc", default='mpicc',  help="The C compiler to use [default: mpicc]")
    parser.add_argument("--ifmpi", default='true', choices=['true', 'false'], help="Enable/disable parallel tests with MPI [default: true]")
    parser.add_argument("--nprocs", default='4', help="Number of processes to use for MPI tests [default: 4]")
    parser.add_argument("-v", "--verbose", action='store_true', help="Enable verbose output")
 
    args = parser.parse_args()

    # # Set environment
    os.environ['CC'] = args.cc
    os.environ['FC'] = args.f77
    os.environ['IFMPI'] = args.ifmpi
    os.environ['PARALLEL_PROCS'] = args.nprocs
    if args.verbose:
        os.environ['VERBOSE_TESTS'] = 'true'
        ut_verbose = 2
    else:
        os.environ['VERBOSE_TESTS'] = 'false'
        ut_verbose = 1

    testList = (
               FsHydro,
               Axi, 
               Eddy_Neknek,
               Eddy_EddyUv,
               Eddy_LegacySize, 
               Benard_Ray9, 
               KovStokes, 
               LowMachTest, 
               MvCylCvode, 
               CmtInviscidVortex,
               DoubleShear,
               Ethier,
               LinCav_Dir,
               LinCav_Adj,
               IO_Test,
               InclDef   
               ) 

    suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
