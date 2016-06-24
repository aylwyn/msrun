#!/usr/bin/python
# Aylwyn Scally 2016

import sys
#import getopt
#import subprocess
import os
import os.path
import logging
from logging import error, warning, info, debug, critical
import gzip
from itertools import combinations
#import locale
#import string
#import decimal
import argparse
from numpy import log

#sys.path.insert(0, '/home/aos21/msprime')
import msprime

import aosutils

#p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p = argparse.ArgumentParser()
p.add_argument('-n', '--nreps', type = int, default = 1, help = 'number of repetitions (default = 1)')
p.add_argument('-s', '--nsamps', type = int, default = 2, help = 'number of samples per rep (default = 4)')
p.add_argument('--npops', type = int, default = 1, help = 'number of populations (default = 1)')
p.add_argument('-p', action = 'append', default = [], nargs = '*', help = 'add present-day population: -p samples [Ne [growth_rate]]; (overrides -s, --npops)', metavar = ('samples', 'Ne'))#, 'growth_rate'))
p.add_argument('-l', '--seqlen', type = int, default = 10000, help = 'sequence length simulated (default = 10000)')
p.add_argument('-u', '--mugen', type = float, default = 1.25e-8, help = 'per-generation mutation rate (default = 1.25e-8)')
p.add_argument('-N', '--N0', type = int, default = 20000, help = 'base effective population size (default = 20000)')
p.add_argument('-r', '--recgen', type = float, default = 0.0, help = 'per-generation recombination rate (default = 0.0)')
p.add_argument('-g', '--tgen', type = float, default = 29.0, help = 'generation time (y) (default = 29.0)')
#p.add_argument('--eN', action='append', default = [], nargs = 2, help='global Ne history', metavar = ('time', 'Ne'))
#p.add_argument('--eG', action='append', default = [], nargs = 2, help='global growth history', metavar = ('time', 'alpha'))
#p.add_argument('--eNG', action='append', default = [], nargs = 4, help='global growth history in terms of time intervals and Ne change (Note: t2 is older than t1)', metavar = ('time1', 'time2', 'Ne1', 'Ne2'))
#p.add_argument('--en', action='append', default = [], nargs = 3, help='population Ne history', metavar = ('time', 'pop', 'Ne'))
#p.add_argument('--eg', action='append', default = [], nargs = 3, help='population growth history', metavar = ('time', 'pop', 'alpha'))
#p.add_argument('--eng', action='append', default = [], nargs = 5, help='population growth history in terms of time intervals and Ne change (Note: t2 is older than t1)', metavar = ('pop', 'time1', 'time2', 'Ne1', 'Ne2'))
p.add_argument('--ej', action='append', default = [], nargs = 3, help='population merge_history', metavar = ('time', 'from_pop', 'to_pop'))
p.add_argument('--em', action='append', default = [], nargs = 4, help='migration_history', metavar = ('time', 'to_pop', 'from_pop', 'to_pop_mig_frac'))

p.add_argument('-R', '--recfile', help='HapMap-format recombination rate file (ignores -r)' )
#p.add_argument('--mst', type = float, help = 'ms theta value; overrides --mugen if set, otherwise ms_theta = 4 * MUGEN * N0 * SEQLEN (e.g. = 4.8e-4 * SEQLEN for human)')
#p.add_argument('--msr', type = float, help = 'ms rho value; overrides --recgen if set, otherwise ms_rho = 4 * RECGEN * N0 * SEQLEN (e.g. = 4.0e-4 * SEQLEN for human)')
p.add_argument('--cmd', action='store_true', default = False, help = 'include command in output')
#p.add_argument('-T', '--trees', action='store_true', default = False, help = 'include trees in output')
p.add_argument('--mrca', action='store_true', default = False, help = 'include TMRCA in ms output')
p.add_argument('--diversity', action='store_true', default = False, help = 'output expected and measured pairwise diversity')
p.add_argument('--haplotypes', action='store_true', default = False, help = 'output haplotypes for each sample')
p.add_argument('--segsites', action='store_true', default = False, help = 'output alleles at segregating sites')
p.add_argument('--outfile', action='store_true', default = False, help = 'redirect output to file')
p.add_argument('--prefix', default = '', help='output file prefix')
p.add_argument('--suffix', default= 'ms', help='output file suffix')
p.add_argument('--encode', action='store_true', default = False, help = 'encode demographic parameters in output file name (sets --outfile)')
#p.add_argument('--encode_all', action='store_true', default = False, help = 'encode all arguments in output file name (sets --outfile)')
#p.add_argument('--batch', action='store_true', default = False, help = 'output submit.py nohup command')
p.add_argument('--sf', type=int, default = 4, help = 'sig figs to use for numbers')
p.add_argument('--demdebug', action='store_true', default = False, help= 'run DemographyDebugger')
p.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
p.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)

args = p.parse_args()

def fnum(x): return aosutils.fnum(x, args.sf)

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

cmdname = 'msp'

encmd = []
#encmd_t = [] # time-ordered commands
popconfig = []
demevents = []

if not args.p:
	popsize = args.nsamps/args.npops
	if not popsize * args.npops == args.nsamps:
		warning('nsamps %d not a multiple of npops %d; setting nsamps = %d' % (args.nsamps, args.npops, popsize * args.npops))
		args.nsamps = popsize * args.npops
	if args.npops > 1:
		info('distributing %d samples among %d populations' %  (args.nsamps, args.npops))
	for pop in range(args.npops):
		popconfig.append(msprime.PopulationConfiguration(sample_size=popsize, initial_size=args.N0))
		encmd.append('-p %d' % popsize)
else:
	for x in args.p:
		if len(x) == 1:
			popconfig.append(msprime.PopulationConfiguration(sample_size=int(x[0]), initial_size=args.N0))
		elif len(x) == 2:
			popconfig.append(msprime.PopulationConfiguration(sample_size=int(x[0]), initial_size=eval(x[1])))
		elif len(x) == 3:
			popconfig.append(msprime.PopulationConfiguration(sample_size=int(x[0]), initial_size=eval(x[1]), growth_rate=eval(x[2])))
		else:
			error('max 3 arguments for -p: samples [initial_Ne [growth_rate]]')
			sys.exit(1)
		encmd.append('-p ' + '_'.join(x))
	args.npops = len(popconfig)

migmatrix = [[0] * args.npops] * args.npops

if args.recfile:
	recmap = msprime.RecombinationMap.read_hapmap(args.recfile)
	args.seqlen = int(recmap.get_positions()[-1])
	encmd.append('-R ' + args.recfile)
#else:
#	encmd.append('-r %s %d' % (fnum(args.msr), args.seqlen))

#if args.trees:
#	encmd.append('-T')
#if args.mrca:
#	encmd.append('-L')

for x in args.ej:
	demevents.append(msprime.MassMigration(time=eval(x[0]), source = int(x[1]), destination = int(x[2]), proportion=1.0))
	encmd.append('--ej ' + '_'.join(x))

#for x in args.eN:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	Nev = float(eval(x[1])) / args.N0
#	encmd_t.append((tev, '-eN %s %s' % (fnum(tev), fnum(Nev))))
#for x in args.eG:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	alpha = float(eval(x[1]))
#	encmd_t.append((tev, '-eG %s %s' % (fnum(tev), fnum(alpha))))
#for x in args.eNG:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	tev2 = eval(x[1]) / (4 * args.N0 * args.tgen)
#	Nev1 = float(eval(x[2])) / args.N0
#	Nev2 = float(eval(x[3])) / args.N0
#	alpha = -(1 / (tev2 - tev)) * log(Nev2 / Nev1)
#	encmd_t.append((tev, '-eN %s %s -eG %s %s' % (fnum(tev), fnum(Nev1), fnum(tev), fnum(alpha))))
#	encmd_t.append((tev2, '-eG %s 0' % (fnum(tev2))))
#for x in args.en:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	pnum = int(x[1])
#	Nev = float(eval(x[2])) / args.N0
#	encmd_t.append((tev, '-en %s %d %s' % (fnum(tev), pnum, fnum(Nev))))
#for x in args.eg:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	pnum = int(x[1])
#	alpha = float(eval(x[2]))
#	encmd_t.append((tev, '-eg %s %d %s' % (fnum(tev), pnum, fnum(alpha))))
#for x in args.eng:
#	pnum = int(x[0])
#	tev = eval(x[1]) / (4 * args.N0 * args.tgen)
#	tev2 = eval(x[2]) / (4 * args.N0 * args.tgen)
#	Nev1 = float(eval(x[3])) / args.N0
#	Nev2 = float(eval(x[4])) / args.N0
#	alpha = -(1 / (tev2 - tev)) * log(Nev2 / Nev1)
#	encmd_t.append((tev, '-en %s %d %s -eg %s %d %s' % (fnum(tev), pnum, fnum(Nev1), fnum(tev), pnum, fnum(alpha))))
#	encmd_t.append((tev2, '-eg %s %d 0' % (fnum(tev2), pnum)))
#for x in args.em:
#	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
#	p1, p2 = x[1:3]
#	mig = eval(x[3]) * 4 * args.N0
#	encmd_t.append((tev, '-em %s %s %s %s' % (fnum(tev), p1, p2, fnum(mig))))

#encmd_t.sort()
#encmd.append(' '.join([x[1] for x in encmd_t]))

#outargs = ' '.join([args.msargs, ' '.join(encmd)])
outargs = ' '.join(encmd)
outname = args.prefix
if args.encode:
	outname += '_'.join(outargs.split()).replace('_-', '-')
if args.suffix:
	outname += '.' + args.suffix

if args.outfile or args.encode:
	boutname = outname + '.bout'

cmd = ' '.join(['msprun.py'] + sys.argv[1:])

if args.demdebug:
	dp = msprime.DemographyDebugger(Ne=args.N0, population_configurations=popconfig, migration_matrix=migmatrix, demographic_events=demevents)
	dp.print_history()
else:
#	info('running \'%s\'' % (cmd))
	if not args.sim:
		if args.seqlen:
			if args.recfile:
				ts = msprime.simulate(Ne=args.N0, population_configurations=popconfig, migration_matrix=migmatrix, demographic_events=demevents, mutation_rate=args.mugen, recombination_map=recmap)
			else:
				ts = msprime.simulate(Ne=args.N0, population_configurations=popconfig, migration_matrix=migmatrix, demographic_events=demevents, mutation_rate=args.mugen, recombination_rate=args.recgen, length=args.seqlen)
		else:
			ts = msprime.simulate(Ne=args.N0, population_configurations=popconfig, migration_matrix=migmatrix, demographic_events=demevents)

		if args.cmd:
			sys.stdout.write('CMD\t%s\n' % (cmd))

		if args.mrca:
			for tree in ts.trees():
				sys.stdout.write('MRCA\t%f\n' % tree.get_time(tree.get_root()) / (4 * args.N0))

		if args.diversity:
			sys.stdout.write('SEQLEN\t%d\n' % args.seqlen)
			sys.stdout.write('PI\t%f\n' % (float(ts.get_pairwise_diversity())/args.seqlen))
	#		sys.stdout.write('THETA %f\n' % (4 * args.N0 * args.mugen))
			if args.npops > 1:
				pi = [0.0] * args.npops
				for pop in range(args.npops):
					pi[pop] = float(ts.get_pairwise_diversity(ts.get_samples(pop))) / args.seqlen
					sys.stdout.write('PI_%d\t%f\n' % (pop, pi[pop]))
			sys.stdout.write('SEGSITES\t%d\n' % ts.get_num_mutations())


		if args.haplotypes:
			for ih, hap in enumerate(ts.haplotypes()):
				sys.stdout.write('SAMPLE\t%d\t%s\n' % (ih, hap))

		if args.segsites:
			alsep = ''
			for position, variant in ts.variants():
				sys.stdout.write('SITE\t%d\t%s\n' % (position, variant))
	#		hap = list(ts.haplotypes())
	#		for im, mutn in enumerate(ts.mutations()):
	#			sys.stdout.write('SITE\t%d\t%s\n' % (mutn[0], alsep.join([x[im] for x in hap])))
