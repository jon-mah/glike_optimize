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
            help=('Lenght of sequence to be simulated.'),
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
            '--NUM_TREES', type=int,
            dest='NUM_TREES',
            help=('Number of equally distant trees to simulate.'),
            default=10)
        parser.add_argument(
            '--OPTIMIZER', type=str,
            dest='OPTIMIZER',
            help=('Optimization strategy for fitting ARGs.'),
            default='maximize',
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
        NUM_TREES = args['NUM_TREES']
        OPTIMIZER = args['OPTIMIZER']
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
        logfile = '{0}{1}simulate_ARG.log'.format(args['outprefix'], underscore)
        output_posfile = \
            '{0}{1}NH_positions.json'.format(
                args['outprefix'], underscore)
        to_remove = [logfile, output_posfile]
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
        demography = glike.native_hawaiians_demography()
        arg = msprime.sim_ancestry(
            {"admix": N_SAMPLES},
            sequence_length=SEQUENCE_LENGTH,
            recombination_rate=RECOMBINATION_RATE,
            demography=demography,
            ploidy=1,
            random_seed=SEED,
        )

        trees = [
            arg.at((i + 0.5) * SEQUENCE_LENGTH // NUM_TREES).copy()
            for i in range(NUM_TREES)
        ]
        i = 0
        for tree in trees:
            # Save ARG
            output_treefile = \
                '{0}{1}NH_ARG_{2}.trees'.format(
                    args['outprefix'], underscore, i)
            arg_path = output_treefile
            arg.dump(arg_path)
            i = i + 1

        # Compute and save tree positions
        N_TREES = 20
        step = int(SEQUENCE_LENGTH) // (N_TREES + 1)
        positions = list(range(step, int(SEQUENCE_LENGTH), step))[:N_TREES]
        positions_path = output_posfile
        with open(positions_path, "w") as f:
            json.dump(positions, f)

        logger.info('Finished simulating ARG.')

        # Demographic model
        x_true = {'t1':19, 't2':411, 't3':1040, 't4':2004, 'r1':0.0, 
                  'r2':0.198, 'r3':0.334, 'N_admix':35682, 'N_afr':10000, 
                  'N_eur':13388, 'N_asia':25234, 'N_pol':15695, 'N_aa':2702, 
                  'N_ooa':2470, 'N_anc':2665, 'gr':0.078}
        true_demo = glike.native_hawaiians_demo(**x_true)
        true_demo.print()

        logp_true = glike.glike_trees(trees, true_demo)

        def fun(t1, t2, t3, t4, r1, r2, r3, N_admix, N_afr, N_eur, N_asia, 
                N_pol, N_aa, N_ooa, N_anc, gr):
            demo = glike.native_hawaiians_demo(t1, t2, t3, t4, r1, r2, r3, 
                                        N_admix, N_afr, N_eur, N_asia, 
                                        N_pol, N_aa, N_ooa, N_anc, gr) 
            return glike.glike_trees(trees, demo)

        # x0 = {'t1':10, 't2': 100, 't3': 1000, 't4': 2000,
        #       'r1':0.25, 'r2':0.25, 'r3':0.25,
        #       'N_admix': 10000, 'N_afr': 10000, 'N_eur': 10000,
        #       'N_asia': 10000, 'N_pol': 10000, 'N_aa': 10000,
        #       'N_ooa': 10000, 'N_anc': 10000, 'gr': 0.1}
        x0 = x_true
        bounds = [(1, "t2"), ("t1", "t3"), ("t2", "t4"), ("t3", 1e4), 
                  (0.0, 1.0), (0.0, 1.0), (0.0, 1.0), 
                  (100, 100000), (100, 100000), (100, 100000), (100, 100000), 
                  (100, 100000), (100, 100000), (100, 100000), (100, 100000),
                  (0, 0.5)]

        logger.info('Starting glike optimization.')
        t_start = time.time()
        if OPTIMIZER == 'CMA_ES':
            x, logp = estimate.maximize_CMA_ES(fun, x0, bounds = bounds)
        else: 
            x, logp = glike.maximize(fun, x0, bounds = bounds)
        elapsed = time.time() - t_start
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)

        logger.info('Finished fitting ARG.')
        logger.info(f"Finished in {hours}h {minutes}m {elapsed % 60:.1f}s")
        logger.info(f"Estimated: {x}")
        logger.info(f"logp = {logp}  (true = {logp_true})")
        logger.info('Pipeline executed succesfully.')



if __name__ == '__main__':
    SimulateARG().main()

