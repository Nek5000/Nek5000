#!/usr/bin/env python
from lib.nekTestCase import *
from unittest import skip

###############################################################################
#  axi: axi.rea
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

        self.build_tools(['genbox', 'genmap'])
        self.run_genbox()
        self.mvn('box', self.__class__.case_name)
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        pres = self.get_value_from_log('PRES ', column=-4)
        self.assertAlmostEqualDelayed(pres, target_val=0., delta=76., label='PRES')

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=2, label='total solver time')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        pres = self.get_value_from_log('PRES ', column=-4)
        self.assertAlmostEqualDelayed(pres, target_val=0., delta=76., label='PRES')

        self.assertDelayedFailures()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        u_press = self.get_value_from_log('U-Press ', column=-5)
        self.assertAlmostEqualDelayed(u_press, target_val=0., delta=104., label='U-Press')

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=4, label='total solver time')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        u_press = self.get_value_from_log('U-Press ', column=-5)
        self.assertAlmostEqualDelayed(u_press, target_val=0., delta=104., label='U-Press')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()
#
# ####################################################################
# #  benard: ray_9.rea, ray_dd.rea, ray_dn.rea, ray_nn.rea
# ####################################################################

class Benard_Ray9(NekTestCase):
    example_subdir = 'benard'
    case_name = 'ray_9'

    def setUp(self):
        self.size_params = dict (
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lelg      = '500',
        )

        self.build_tools(['genmap'])
        self.run_genmap(rea_file='ray_9')

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9', step_limit=1000)

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=30., label='total solver time')

        gmres = self.get_value_from_log('gmres ', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=23., label='gmres')

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2']='lx1'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9', step_limit=1000)

        gmres = self.get_value_from_log('gmres ', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=23., label='gmres')

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9', step_limit=1000)

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=40., label='total solver time')

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=11., label='gmres')

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2']='lx1-2'
        self.config_size()
        self.build_nek(usr_file='ray_9')
        self.run_nek(rea_file='ray_9', step_limit=1000)

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=11., label='gmres')

    def tearDown(self):
        self.move_logs()

# ####################################################################
# #  eddy; eddy_uv.rea, amg_eddy.rea, htps_ed.rea
# ####################################################################

# TODO: implement eddy for amg_eddy.rea, htps_ed.rea

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
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=34., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=80, label='total solver time')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-7,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=34., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.007702E-07, delta=1E-08, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=6.489061E-07, delta=1E-08, label='Y err')

        self.assertDelayedFailures()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=22., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.759103E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=7.842019E-05, delta=1E-06, label='Y err')

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, 0.1, delta=80, label='total solver time')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6,)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=22., label='gmres')

        xerr = self.get_value_from_log('X err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(xerr, target_val=6.759103E-05, delta=1E-06, label='X err')

        yerr = self.get_value_from_log('Y err', column=-6, row=-1)
        self.assertAlmostEqualDelayed(yerr, target_val=7.842019E-05, delta=1E-06, label='Y err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()


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

    @pn_pn_serial
    def test_PnPn_Serial(self):
        # Update SIZE parameters for PnPn
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

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=80, label='total solver time')

        self.assertDelayedFailures()

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

        self.assertDelayedFailures()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
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

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, 0.1, delta=80, label='total solver time')

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

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()
####################################################################
#  kov_st_state; kov_st_stokes.rea
####################################################################

class KovStState(NekTestCase):
    # Note: Legacy Analysis.py script only checked Pn-Pn-2 test cases
    example_subdir = 'kov_st_state'
    case_name = 'kov_st_stokes'

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

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        # solver_time = self.get_value_from_log(label='total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=5, label='total solver time')

        err = self.get_value_from_log(label='err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err, target_val=8.55641E-10, delta=1e-11, label='err')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        err = self.get_value_from_log(label='err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err, target_val=8.55641E-10, delta=1e-11, label='err')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()

####################################################################
#  lowMach_test; lowMach_test.rea
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

        # Tweak the .rea file and run genmap
        from re import sub
        cls = self.__class__
        rea_path = os.path.join(self.examples_root, cls.example_subdir, cls.case_name + '.rea')
        with open(rea_path, 'r') as f:
            lines = [sub(r'^.*IFNAV.*$', '  T T IFNAV & IFADVC', l) for l in f]
        with open(rea_path, 'w') as f:
            f.writelines(lines)
        self.run_genmap()

    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        # solver_time = self.get_value_from_log(label='total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=40, label='total solver time')

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=100, label='gmres')

        vx = self.get_value_from_log(label='ERROR VX', column=-5, row=-1)
        self.assertAlmostEqualDelayed(vx, target_val=2.4635E-09, delta=1e-10, label='VX')

        errt = self.get_value_from_log(label='ERROR T', column=-5, row=-1)
        self.assertAlmostEqualDelayed(errt, target_val=4.5408E-12, delta=1e-13, label='T')

        qtl = self.get_value_from_log(label='ERROR QTL', column=-5, row=-1)
        self.assertAlmostEqualDelayed(qtl, target_val=2.6557E-06, delta=1e-07, label='QTL')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=200)

        gmres = self.get_value_from_log(label='gmres', column=-7)
        self.assertAlmostEqualDelayed(gmres, target_val=0, delta=100, label='gmres')

        vx = self.get_value_from_log(label='ERROR VX', column=-5, row=-1)
        self.assertAlmostEqualDelayed(vx, target_val=2.4635E-09, delta=1e-10, label='VX')

        errt = self.get_value_from_log(label='ERROR T', column=-5, row=-1)
        self.assertAlmostEqualDelayed(errt, target_val=4.5408E-12, delta=1e-13, label='T')

        qtl = self.get_value_from_log(label='ERROR QTL', column=-5, row=-1)
        self.assertAlmostEqualDelayed(qtl, target_val=2.6557E-06, delta=1e-07, label='QTL')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()


####################################################################
#  mv_cyl with CVODE
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
    def test_PnPn_Parallel_Steps1e3(self):
        if not "CVODE" in self.pplist:
            self.fail("\"CVODE\" is not listed in $PPLIST. This test cannot be run.".format(self.id()))

        self.log_suffix += '.steps_1e3'
        self.config_parfile({'GENERAL' : {'numSteps' : '1e3', 'dt' : '1e-3'}})
        self.build_nek()
        self.run_nek()

        err3 = self.get_value_from_log('err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err3, target_val=0.1743079E-03, delta=1e-6, label='err (column -3)')

        err2 = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=0.6348537E-06, delta=1e-9, label='err (column -2)')

        self.assertDelayedFailures()

    @pn_pn_parallel
    def test_PnPn_Parallel_Steps1e4(self):
        if not "CVODE" in self.pplist:
            self.fail("\"CVODE\" is not listed in $PPLIST. This test cannot be run.".format(self.id()))

        self.log_suffix += '.steps_1e4'
        self.config_parfile({'GENERAL' : {'numSteps' : '1e4', 'dt' : '1e-4'}})
        self.build_nek()
        self.run_nek()

        err3 = self.get_value_from_log('err', column=-3, row=-1)
        self.assertAlmostEqualDelayed(err3, target_val=0.1693853E-05, delta=1e-8, label='err (column -3)')

        err2 = self.get_value_from_log('err', column=-2, row=-1)
        self.assertAlmostEqualDelayed(err2, target_val=0.6344692E-09, delta=1e-12, label='err (column -2)')

        self.assertDelayedFailures()

    def tearDown(self):
        self.move_logs()


####################################################################
#  var_vis; var_vis.rea
####################################################################

class VarVis(NekTestCase):
    example_subdir = 'var_vis'
    case_name = 'st2'

    def setUp(self):
        self.size_params = dict(
            ldim      = '2',
            lx1       = '8',
            lxd       = '12',
            lx2       = 'lx1-2',
            lx1m      = 'lx1',
            lelg      = '500',
        )
        self.build_tools(['genmap'])
        self.run_genmap()

    @unittest.expectedFailure
    @pn_pn_serial
    def test_PnPn_Serial(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        # TODO: This fails here and in legacy tests
        phrase = self.get_phrase_from_log('ABORT: ')
        self.assertIsNotNullDelayed(phrase, label='ABORT')
        self.assertDelayedFailures()

    @unittest.expectedFailure
    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        self.size_params['lx2'] = 'lx1'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        # TODO: This fails here and in legacy tests
        phrase = self.get_phrase_from_log('ABORT: ')
        self.assertIsNotNullDelayed(phrase, label='ABORT')
        self.assertDelayedFailures()

    @pn_pn_2_serial
    def test_PnPn2_Serial(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        # solver_time = self.get_value_from_log('total solver time', column=-2)
        # self.assertAlmostEqualDelayed(solver_time, target_val=0.1, delta=30, label='total solver time')

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=19, label='gmres')

        self.assertDelayedFailures()

    @pn_pn_2_parallel
    def test_PnPn2_Parallel(self):
        self.size_params['lx2'] = 'lx1-2'
        self.config_size()
        self.build_nek()
        self.run_nek(step_limit=None)

        gmres = self.get_value_from_log('gmres ', column=-6)
        self.assertAlmostEqualDelayed(gmres, target_val=0., delta=19, label='gmres')

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
        cls = self.__class__
        try:
            os.remove(os.path.join(self.examples_root, cls.example_subdir, 'l2norms.dat'))
        except OSError:
            pass

    @pn_pn_serial
    def test_PnPn_Serial(self):
        if "CMTNEK" not in self.pplist:
            self.fail("\"CMTNEK\" is not listed in $PPLIST. This test cannot be run.".format(self.id()))
        self.build_nek()
        self.run_nek(step_limit=1000)
        self.diff_l2norms()

    @pn_pn_parallel
    def test_PnPn_Parallel(self):
        if "CMTNEK" not in self.pplist:
            self.fail("\"CMTNEK\" is not listed in $PPLIST. This test cannot be run.".format(self.id()))
        self.build_nek()
        self.run_nek(step_limit=1000)
        self.diff_l2norms()


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
    os.environ['F77'] = args.f77
    os.environ['IFMPI'] = args.ifmpi
    os.environ['PARALLEL_PROCS'] = args.nprocs
    if args.verbose:
        os.environ['VERBOSE_TESTS'] = 'true'
        ut_verbose = 2
    else:
        os.environ['VERBOSE_TESTS'] = 'false'
        ut_verbose = 1

    testList = (
               Axi, 
               Eddy_EddyUv, 
               Benard_Ray9, 
               Eddy_LegacySize, 
               KovStState, 
               LowMachTest, 
               MvCylCvode, 
               VarVis, 
               CmtInviscidVortex
               ) 

    suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
