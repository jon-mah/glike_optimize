""" Simulate Ancestral Recombination Graph (ARG) using a given set of parameters.

JCM 20260801
"""


import sys
import os
import logging
import time
import argparse
import warnings

sys.path.append('../glike/')
import estimate

import glike
import msprime
import json
import tskit
import numpy
import scipy.stats.distributions
import scipy.integrate
import scipy.optimize
import pandas as pd
import math

class ArgumentParserNoArgHelp(argparse.ArgumentParser):
    """Like *argparse.ArgumentParser*, but prints help when no arguments."""

    def error(self, message):
        """Print error message, then help."""
        sys.stderr.write('error: %s\n\n' % message)
        self.print_help()
        sys.exit(2)


class SimulateARG():
    """Wrapper class to allow functions to reference each other."""

    def ExistingFile(self, fname):
        """Return *fname* if existing file, otherwise raise ValueError."""
        if os.path.isfile(fname):
            return fname
        else:
            raise ValueError("%s must specify a valid file name" % fname)

    ########## Three popopulation out of Africa
    def three_pop_out_of_africa_demo_no_m(self,
        t1 = 848, t2 = 5600, t3 = 8800,
        N_anc = 7300, N_yri = 12300, N_ooa = 2100,
        N_ceu = 1000, N_chb = 510,
        gr_ceu = 0.004, gr_chb = 0.0055,
    ):
        # Migration matrices for each demographic epoch
        Q0 = numpy.array([
            [0,  0,                 0],
            [ 0,            0,      0],
            [ 0,             0,                0]
        ]) if max(0, 0, 0) > 0 else None

        Q1 = numpy.array([
            [0,  0],
            [ 0, 0]
        ]) if 0 > 0 else None

        demo = glike.Demo()

        # ---------------------------------------------------------
        # Epoch 0 -> t1
        #
        # Three populations:
        #   yri
        #   ceu
        #   chb
        #
        # CEU and CHB are growing populations.
        # ---------------------------------------------------------

        demo.add_phase(
            glike.Phase(
                0,
                1e-6,
                [
                    1/N_yri,
                    1/(N_ceu * math.exp(gr_ceu * t1)),
                    1/(N_chb * math.exp(gr_chb * t1))
                ],
                grs=[0, gr_ceu, gr_chb],
                populations=["yri", "ceu", "chb"]
            )
        )

        demo.add_phase(
            glike.Phase(
                1e-6,
                t1,
                [
                    1/N_yri,
                    1/(N_ceu * math.exp(gr_ceu * t1)),
                    1/(N_chb * math.exp(gr_chb * t1))
                ],
                grs=[0, gr_ceu, gr_chb],
                Q=Q0,
                populations=["yri", "ceu", "chb"]
            ),
            discretize=200
        )

        # ---------------------------------------------------------
        # t1 -> t2
        #
        # CEU and CHB merge into the OOA population.
        #
        # Going backwards in time:
        #
        #     YRI       -> YRI
        #     CEU + CHB -> OOA
        #
        # After t1, CEU/CHB migration is replaced by YRI-OOA
        # migration.
        # ---------------------------------------------------------

        P_ceu_chb_split = numpy.array([
            [1, 0],
            [0, 1],
            [0, 1]
        ])

        demo.add_phase(
            glike.Phase(
                t1,
                t2,
                [
                    1/N_yri,
                    1/N_ooa
                ],
                [0, 0],
                P=P_ceu_chb_split,
                Q=Q1,
                populations=["yri", "ooa"]
            ),
            discretize=500
        )

        # ---------------------------------------------------------
        # t2 -> t3
        #
        # OOA merges into the ancestral YRI population.
        #
        # Going backwards:
        #
        #     YRI + OOA -> ancestral population
        # ---------------------------------------------------------

        P_ooa_split = numpy.array([
            [1],
            [1]
        ])

        demo.add_phase(
            glike.Phase(
                t2,
                t3,
                [1/N_ooa],
                P=P_ooa_split,
                populations=["anc"]
            )
        )

        # ---------------------------------------------------------
        # t3 -> infinity
        #
        # The ancestral population changes size to N_anc.
        # ---------------------------------------------------------

        demo.add_phase(
            glike.Phase(
                t3,
                math.inf,
                [1/N_anc],
                populations=["anc"]
            )
        )

        return demo

    def three_pop_out_of_africa_demography_no_m(self,
        t1=848,
        t2=5600,
        t3=8800,
        N_anc=7300,
        N_yri=12300,
        N_ooa=2100,
        N_ceu=1000,
        N_chb=510,
        gr_ceu=0.004,
        gr_chb=0.0055,
    ):
        demography = msprime.Demography()

        # Present-day populations
        demography.add_population(
            name="yri",
            initial_size=N_yri,
        )

        demography.add_population(
            name="ceu",
            initial_size=N_ceu,
            growth_rate=gr_ceu,
        )

        demography.add_population(
            name="chb",
            initial_size=N_chb,
            growth_rate=gr_chb,
        )

        # t1: CEU and CHB split.
        # Going backwards in time, CHB merges into CEU.
        demography.add_mass_migration(
            time=t1,
            source="chb",
            dest="ceu",
            proportion=1,
        )

        # After the CEU/CHB split, CEU and CHB no longer exchange migrants.
        demography.set_symmetric_migration_rate(
            populations=["ceu", "chb"],
            rate=0,
        )

        # OOA ancestral population
        demography.add_population(
            name="ooa",
            initial_size=N_ooa,
        )

        # t2: OOA population joins the YRI population.
        # Going backwards in time, OOA merges into YRI.
        demography.add_mass_migration(
            time=t2,
            source="ceu",
            dest="yri",
            proportion=1,
        )

        # At t2, CEU/CHB ancestry is represented by the OOA population.
        # Set the appropriate ancestral population size.
        demography.add_population_parameters_change(
            time=t2,
            initial_size=N_ooa,
            growth_rate=0,
            population="yri",
        )

        # t3: ancestral African population size changes.
        demography.add_population_parameters_change(
            time=t3,
            initial_size=N_anc,
            growth_rate=0,
            population="yri",
        )

        return demography


    ########## Three popopulation out of Africa
    def three_pop_out_of_africa_demo(self,
        t1 = 848, t2 = 5600, t3 = 8800,
        N_anc = 7300, N_yri = 12300, N_ooa = 2100,
        N_ceu = 1000, N_chb = 510,
        gr_ceu = 0.004, gr_chb = 0.0055,
        m_yri_ooa = 25e-5, m_yri_ceu = 3e-5,
        m_yri_chb = 1.9e-5, m_ceu_chb = 9.6e-5
    ):
        # Migration matrices for each demographic epoch
        Q0 = numpy.array([
            [-m_yri_ceu-m_yri_chb,  m_yri_ceu,                 m_yri_chb],
            [ m_yri_ceu,            -m_yri_ceu-m_ceu_chb,      m_ceu_chb],
            [ m_yri_chb,             m_ceu_chb,                -m_yri_chb-m_ceu_chb]
        ]) if max(m_yri_ceu, m_yri_chb, m_ceu_chb) > 0 else None

        Q1 = numpy.array([
            [-m_yri_ooa,  m_yri_ooa],
            [ m_yri_ooa, -m_yri_ooa]
        ]) if m_yri_ooa > 0 else None

        demo = glike.Demo()

        # ---------------------------------------------------------
        # Epoch 0 -> t1
        #
        # Three populations:
        #   yri
        #   ceu
        #   chb
        #
        # CEU and CHB are growing populations.
        # ---------------------------------------------------------

        demo.add_phase(
            glike.Phase(
                0,
                1e-6,
                [
                    1/N_yri,
                    1/(N_ceu * math.exp(gr_ceu * t1)),
                    1/(N_chb * math.exp(gr_chb * t1))
                ],
                grs=[0, gr_ceu, gr_chb],
                populations=["yri", "ceu", "chb"]
            )
        )

        demo.add_phase(
            glike.Phase(
                1e-6,
                t1,
                [
                    1/N_yri,
                    1/(N_ceu * math.exp(gr_ceu * t1)),
                    1/(N_chb * math.exp(gr_chb * t1))
                ],
                grs=[0, gr_ceu, gr_chb],
                Q=Q0,
                populations=["yri", "ceu", "chb"]
            ),
            discretize=200
        )

        # ---------------------------------------------------------
        # t1 -> t2
        #
        # CEU and CHB merge into the OOA population.
        #
        # Going backwards in time:
        #
        #     YRI       -> YRI
        #     CEU + CHB -> OOA
        #
        # After t1, CEU/CHB migration is replaced by YRI-OOA
        # migration.
        # ---------------------------------------------------------

        P_ceu_chb_split = numpy.array([
            [1, 0],
            [0, 1],
            [0, 1]
        ])

        demo.add_phase(
            glike.Phase(
                t1,
                t2,
                [
                    1/N_yri,
                    1/N_ooa
                ],
                [0, 0],
                P=P_ceu_chb_split,
                Q=Q1,
                populations=["yri", "ooa"]
            ),
            discretize=500
        )

        # ---------------------------------------------------------
        # t2 -> t3
        #
        # OOA merges into the ancestral YRI population.
        #
        # Going backwards:
        #
        #     YRI + OOA -> ancestral population
        # ---------------------------------------------------------

        P_ooa_split = numpy.array([
            [1],
            [1]
        ])

        demo.add_phase(
            glike.Phase(
                t2,
                t3,
                [1/N_ooa],
                P=P_ooa_split,
                populations=["anc"]
            )
        )

        # ---------------------------------------------------------
        # t3 -> infinity
        #
        # The ancestral population changes size to N_anc.
        # ---------------------------------------------------------

        demo.add_phase(
            glike.Phase(
                t3,
                math.inf,
                [1/N_anc],
                populations=["anc"]
            )
        )

        return demo

    def three_pop_out_of_africa_demography(self,
        t1=848,
        t2=5600,
        t3=8800,
        N_anc=7300,
        N_yri=12300,
        N_ooa=2100,
        N_ceu=1000,
        N_chb=510,
        gr_ceu=0.004,
        gr_chb=0.0055,
        m_yri_ooa=25e-5,
        m_yri_ceu=3e-5,
        m_yri_chb=1.9e-5,
        m_ceu_chb=9.6e-5,
    ):
        demography = msprime.Demography()

        # Present-day populations
        demography.add_population(
            name="yri",
            initial_size=N_yri,
        )

        demography.add_population(
            name="ceu",
            initial_size=N_ceu,
            growth_rate=gr_ceu,
        )

        demography.add_population(
            name="chb",
            initial_size=N_chb,
            growth_rate=gr_chb,
        )

        # Migration between populations
        demography.set_symmetric_migration_rate(
            populations=["yri", "ceu"],
            rate=m_yri_ceu,
        )

        demography.set_symmetric_migration_rate(
            populations=["yri", "chb"],
            rate=m_yri_chb,
        )

        demography.set_symmetric_migration_rate(
            populations=["ceu", "chb"],
            rate=m_ceu_chb,
        )

        # t1: CEU and CHB split.
        # Going backwards in time, CHB merges into CEU.
        demography.add_mass_migration(
            time=t1,
            source="chb",
            dest="ceu",
            proportion=1,
        )

        # After the CEU/CHB split, CEU and CHB no longer exchange migrants.
        demography.set_symmetric_migration_rate(
            populations=["ceu", "chb"],
            rate=0,
        )

        # OOA ancestral population
        demography.add_population(
            name="ooa",
            initial_size=N_ooa,
        )

        # t2: OOA population joins the YRI population.
        # Going backwards in time, OOA merges into YRI.
        demography.add_mass_migration(
            time=t2,
            source="ceu",
            dest="yri",
            proportion=1,
        )

        # At t2, CEU/CHB ancestry is represented by the OOA population.
        # Set the appropriate ancestral population size.
        demography.add_population_parameters_change(
            time=t2,
            initial_size=N_ooa,
            growth_rate=0,
            population="yri",
        )

        # t3: ancestral African population size changes.
        demography.add_population_parameters_change(
            time=t3,
            initial_size=N_anc,
            growth_rate=0,
            population="yri",
        )

        return demography


    def simulateARGParser(self):
        """Return *argparse.ArgumentParser* for ``simulate_ARG.py``."""
        parser = ArgumentParserNoArgHelp(
            description=(
                '''
                Simulate a given number of equally distant trees using 
                the true demographic parameters of a Native Hawaiian 
                Ancestral Recombination Graph (ARG). We then fit those 
                trees to the demographic model of the NH scenario.
                '''),
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        parser.add_argument(
            '--N_SAMPLES', type=float,
            dest='N_SAMPLES',
            help=('Number of samples to be simulated.'),
            default=1000.0)
        parser.add_argument(
            '--SEQUENCE_LENGTH', type=float,
            dest='SEQUENCE_LENGTH',
            help=('Length of sequence to be simulated.'),
            default=3e7)
        parser.add_argument(
            '--RECOMBINATION_RATE', type=float,
            dest='RECOMBINATION_RATE',
            help=('Rate of recombination.'),
            default=1e-8)
        parser.add_argument(
            '--SEED', type=int,
            dest='SEED',
            help=('Random seed.'),
            default=1)
        parser.add_argument(
            '--MODEL', type=str,
            dest='MODEL',
            help=('Assumed demographic scenario.'),
            default='3G09_no_m')
        parser.add_argument(
            '--KAPPA', type=int,
            dest='KAPPA',
            help=('Kappa parameter for the simulation.'),
            default=10000)
        parser.add_argument(
            '--TRUE_START', action='store_true',
            default=False,
            help=('Initial starting conditions based on simulation.'))
        parser.add_argument(
            '--CMA_params', action='store_true',
            default=False,
            help=('Evaluating parameters based on CMA-ES output.')
        )
        parser.add_argument(
            '--gLike_params', action='store_true',
            default=False,
            help=('Evaluating parameters based on gLike output.')
        )
        parser.add_argument(
            'outprefix', type=str,
            help='The file prefix for the output files.')
        return parser

    def main(self):
        """Execute main function."""
        # Parse command line arguments
        parser = self.simulateARGParser()
        args = vars(parser.parse_args())
        prog = parser.prog

        # Assign arguments
        N_SAMPLES = args['N_SAMPLES']
        SEQUENCE_LENGTH = args['SEQUENCE_LENGTH']
        RECOMBINATION_RATE = args['RECOMBINATION_RATE']
        SEED = args['SEED']
        MODEL = args['MODEL']
        KAPPA = args['KAPPA']
        TRUE_START = args['TRUE_START']
        CMA_params = args['CMA_params']
        gLike_params = args['gLike_params']
        outprefix = args['outprefix']

        # Numpy options
        numpy.set_printoptions(linewidth=numpy.inf)

        # create output directory if needed
        outdir = os.path.dirname(args['outprefix'])
        if outdir:
            if not os.path.isdir(outdir):
                if os.path.isfile(outdir):
                    os.remove(outdir)
                os.mkdir(outdir)

        # Output files: logfile
        # Remove output files if they already exist
        underscore = '' if args['outprefix'][-1] == '/' else '_'
        logfile = '{0}{1}benchmark_ARG_gLike.log'.format(args['outprefix'], underscore)
        if TRUE_START:
            logfile = '{0}{1}benchmark_ARG_true.log'.format(args['outprefix'], underscore)
            output_csv = \
                '{0}{1}num_trees_benchmark_true.csv'.format(
                    args['outprefix'], underscore)
        if CMA_params:
            logfile = '{0}{1}benchmark_ARG_CMA.log'.format(args['outprefix'], underscore)
            output_csv = \
                '{0}{1}num_trees_benchmark_CMA.csv'.format(
                    args['outprefix'], underscore)
        if gLike_params:
            logfile = '{0}{1}benchmark_ARG_gLike.log'.format(args['outprefix'], underscore)
            output_csv = \
                '{0}{1}num_trees_benchmark_gLike.csv'.format(
                    args['outprefix'], underscore)
        to_remove = [logfile, output_csv]
        for f in to_remove:
            if os.path.isfile(f):
                os.remove(f)

        # Set up to log everything to logfile.
        logging.shutdown()
        logging.captureWarnings(True)
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.INFO)
        logger = logging.getLogger(prog)
        warning_logger = logging.getLogger("py.warnings")
        logfile_handler = logging.FileHandler(logfile)
        logger.addHandler(logfile_handler)
        warning_logger.addHandler(logfile_handler)
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
        logfile_handler.setFormatter(formatter)
        logger.setLevel(logging.INFO)

        # print some basic information
        logger.info('Beginning execution of {0} in directory {1}\n'.format(
            prog, os.getcwd()))
        logger.info('Progress is being logged to {0}\n'.format(logfile))
        logger.info('Parsed the following arguments:\n{0}\n'.format(
            '\n'.join(['\t{0} = {1}'.format(*tup) for tup in args.items()])))

        # inherit true params from glike
        if MODEL == '3G09' or MODEL == '3G09_no_m':
            demography = self.three_pop_out_of_africa_demography( 
                t1=848, t2=5600, t3=8800,
                N_anc=7300, N_yri=12300, N_ooa=2100, N_ceu=1000, N_chb=510,
                gr_ceu=0.004, gr_chb=0.0055,
                m_yri_ooa=0, m_yri_ceu=0, m_yri_chb=0, m_ceu_chb=0
            )
            arg = msprime.sim_ancestry(
                {"ceu": N_SAMPLES},
                sequence_length=SEQUENCE_LENGTH,
                recombination_rate=RECOMBINATION_RATE,
                demography=demography,
                ploidy=1,
                random_seed=SEED
            )
        elif MODEL == '3I21':
            demography = glike.neandertal_admixture_demography(
                t1=30, t2=50, t3=73.95, t4=290, 
                N_yri=10000, N_ceu=10000, N_nea = 10000,
                m1=0.029
            )
            arg = msprime.sim_ancestry(
                {"ceu": N_SAMPLES},
                sequence_length=SEQUENCE_LENGTH,
                recombination_rate=RECOMBINATION_RATE,
                demography=demography,
                ploidy=1,
                random_seed=SEED
            )
        elif MODEL == '4A21':
            demography = glike.ancient_europe_demography()
            arg = msprime.sim_ancestry(
                {"bronze": N_SAMPLES},
                sequence_length=SEQUENCE_LENGTH,
                recombination_rate=RECOMBINATION_RATE,
                demography=demography,
                ploidy=1,
                random_seed=SEED
            )
        else:
            demography = glike.native_hawaiians_demography()
            arg = msprime.sim_ancestry(
                {"admix": N_SAMPLES},
                sequence_length=SEQUENCE_LENGTH,
                recombination_rate=RECOMBINATION_RATE,
                demography=demography,
                ploidy=1,
                random_seed=SEED,
            )
        logger.info('Finished simulating ARG.')

        # Demographic model
        if MODEL == '3G09_no_m':
            x_true = {
                "t1": 848, "t2": 5600, "t3": 8800,
                "N_anc": 7300, "N_yri": 12300, "N_ooa": 2100, "N_ceu": 1000, "N_chb": 510,
                "gr_ceu": 0.004, "gr_chb": 0.0055                
            }
            true_demo = self.three_pop_out_of_africa_demo_no_m(**x_true)
        elif MODEL == '3G09':
            x_true = {
                "t1": 848, "t2": 5600, "t3": 8800,
                "N_anc": 7300, "N_yri": 12300, "N_ooa": 2100, "N_ceu": 1000, "N_chb": 510,
                "gr_ceu": 0.004, "gr_chb": 0.0055,
                "m_yri_ooa": 0, "m_yri_ceu": 0, "m_yri_chb": 0, "m_ceu_chb": 0,
            }
            # m1 temporarily set to 0
            true_demo = self.three_pop_out_of_africa_demo(**x_true)
        elif MODEL == '3I21':
            x_true = {
                't1': 30, 't2': 50, 't3': 73.95, 't4': 290,
                'N_yri': 10000, 'N_ceu': 10000, 'N_nea': 10000, 
                'm1': 0.029
            }
            # m1 temporarily set to 0, originally set to 0.029
            true_demo = glike.neandertal_admixture_demo(**x_true)
        elif MODEL == '4A21':
            x_true = {
                't1': 140, 't2': 180, 't3': 200, 't4': 600, 
                't5': 800, 't6': 1500,'r1': 0.5, 'r2': 0.5, 'r3': 0.75,
                'N_ana': 50000, 'N_neo': 50000, 'N_whg': 10000, 
                'N_bronze': 50000, 'N_yam': 5000, 'N_ehg': 10000,
                'N_ne': 5000, 'N_wa': 5000, 'N_ooa': 5000, 'gr': 0.067
            }
            true_demo = glike.ancient_europe_demo(**x_true)
        else: 
            x_true = {
                't1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0, 
                'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000, 
                'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702, 
                'N_ooa':2470, 'N_anc':2665, 'gr':0.078
            }
            true_demo = glike.native_hawaiians_demo(**x_true)

        if TRUE_START:
            if MODEL == '3G09_no_m':
                x_true = {
                    "t1": 848, "t2": 5600, "t3": 8800,
                    "N_anc": 7300, "N_yri": 12300, "N_ooa": 2100, "N_ceu": 1000, "N_chb": 510,
                    "gr_ceu": 0.004, "gr_chb": 0.0055
                }
            elif MODEL == '3G09':
                x_true = {
                    "t1": 848, "t2": 5600, "t3": 8800,
                    "N_anc": 7300, "N_yri": 12300, "N_ooa": 2100, "N_ceu": 1000, "N_chb": 510,
                    "gr_ceu": 0.004, "gr_chb": 0.0055,
                    "m_yri_ooa": 0, "m_yri_ceu": 0, "m_yri_chb": 0, "m_ceu_chb": 0,
                }
            elif MODEL == '3I21':
                x_true = {
                    't1': 30, 't2': 50, 't3': 73.95, 't4': 290,
                    'N_yri': 10000, 'N_ceu': 10000, 'N_nea': 10000,
                    'm1': 0.029
                }
            elif MODEL == '4A21':
                x_true = {
                    't1': 140, 't2': 180, 't3': 200, 't4': 600,
                    't5': 800, 't6': 1500,'r1': 0.5, 'r2': 0.5, 'r3': 0.75,
                    'N_ana': 50000, 'N_neo': 50000, 'N_whg': 10000,
                    'N_bronze': 50000, 'N_yam': 5000, 'N_ehg': 10000,
                    'N_ne': 5000, 'N_wa': 5000, 'N_ooa': 5000, 'gr': 0.067
                }
            else:
                x_true = {
                    't1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0,
                    'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000,
                    'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702,
                    'N_ooa':2470, 'N_anc':2665, 'gr':0.078
                }
        if CMA_params:
            x_true = {'t1': 55.011621001610386, 't2': 133.83139566053288, 
                't3': 38846.11092152734, 'N_anc': 22453.430698161006, 
                'N_yri': 383.05284687097054, 'N_ooa': 525.26487741864, 
                'N_ceu': 293.54608899316605, 'N_chb': 188.844355602063, 
                'gr_ceu': 8.477556224701735e-05, 'gr_chb': 0.01086077550025692}

        if gLike_params:
            x_true = {'t1': 33.66786, 't2': 97.30886, 't3': 97.81987, 
                'N_anc': 34966.222, 'N_yri': 19341.99342, 'N_ooa': 30878.45804, 
                'N_ceu': 13583.28307, 'N_chb': 11372.06889, 'gr_ceu': 0.01875, 
                'gr_chb': 0.03938}

        num_trees = range(10, 10001, 10)
        likelihood = []
        t_start = time.time()
        for num_tree in num_trees:
            logger.info('Evaluating log likelihood of:')
            logger.info(num_tree)
            trees = [
                arg.at((i + 0.5) * SEQUENCE_LENGTH // num_tree).copy()
                for i in range(num_tree)
            ]
            i = 0
            likelihood.append(glike.glike_trees(trees, true_demo))

        logger.info(num_trees)
        logger.info(likelihood)

        data = [num_trees, likelihood]
        df = pd.DataFrame(data)
        df.to_csv(output_csv, index=False)
        
        elapsed = time.time() - t_start
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)

        logger.info('Finished benchmarking ARG.')
        logger.info(f"Finished in {hours}h {minutes}m {elapsed % 60:.1f}s")
        logger.info('Pipeline executed succesfully.')



if __name__ == '__main__':
    SimulateARG().main()

