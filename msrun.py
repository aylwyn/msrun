#!/usr/bin/python
# Aylwyn Scally 2014

import sys
import getopt
import subprocess
#from glob import glob
import os
import os.path
#from time import strftime
import logging
from logging import error, warning, info, debug, critical
import gzip
import locale
import string
import decimal
import argparse

from aosutils import *

p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('-s', '--nsamps', type = int, default = 4, help = 'number of samples per rep')
p.add_argument('-n', '--nreps', type = int, default = 1, help = 'number of repetitions')
p.add_argument('-l', '--seqlen', type = int, default = 10000, help = 'sequence length simulated')
p.add_argument('-p', '--npops', type = int, default = 1, help = 'number of populations')
p.add_argument('-u', '--mugen', type = float, default = 1.25e-8, help = 'per-generation mutation rate')
p.add_argument('-N', '--N0', type = int, default = 1e4, help = 'base effective population size')
p.add_argument('-r', '--recgen', type = float, default = 0.0, help = 'per-generation recombination rate')
p.add_argument('-g', '--tgen', type = float, default = 30.0, help = 'generation time (y)')
p.add_argument('-T', '--trees', action='store_true', default = False, help = 'inlude trees in ms output')
p.add_argument('--mrca', action='store_true', default = False, help = 'inlude TMRCA in ms output')
p.add_argument('--eN', help='global Ne history: "time,Ne ..."' )
p.add_argument('--en', help='population Ne history: "time,pop,Ne ..."' )
p.add_argument('--ej', help='population merge_history: "time,from_pop,to_pop ..."')
p.add_argument('--em', help='migration_history: "time,to_pop,from_pop,to_pop_mig_frac ..."')
p.add_argument('--macs', action='store_true', default = False, help = 'output MaCS command (otherwise use ms)')
p.add_argument('--recfile', help='MaCS-formatted recombination rate file (enforces --macs, ignores -r)' )
p.add_argument('--chrmap', help='chromosomal recombination rate file. Sets SEQLEN equal to length of map and writes recfile to PREFIX.macsrec. (Enforces --macs, ignores -r)' )
p.add_argument('--mst', type = float, help = 'ms theta value; overrides --mugen if set, otherwise ms_theta = 4 * MUGEN * N0 * SEQLEN (e.g. = 4.8e-4 * SEQLEN for human)')
p.add_argument('--msr', type = float, help = 'ms rho value; overrides --recgen if set, otherwise ms_rho = 4 * RECGEN * N0 * SEQLEN (e.g. = 4.0e-4 * SEQLEN for human)')
p.add_argument('--msargs', help = 'ms arguments: "nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <ms_options>" (see MS documentation for more options)')
p.add_argument('--macsargs', help = 'MaCS arguments: "nsamps seqlen -i nreps -t macst [-r macsr] [-I npops pop1_nsamps [pop2_nsamps ...]] <macs_options>" (see MaCS documentation for more options)')
p.add_argument('--outfile', action='store_true', default = False, help = 'redirect output to OUTFILE')
p.add_argument('--prefix', default = 'sim', help='output prefix')
p.add_argument('--suffix', default='ms', help='output suffix')
p.add_argument('--encode_pars', action='store_true', default = False, help = 'encode parameters in output file name')
p.add_argument('--unscale_string', help='print unscaled parameters for simulation output name UNSCALE_STRING')
p.add_argument('--bsub', action='store_true', default = False, help = 'output submit.py bsub command')
p.add_argument('--batch', action='store_true', default = False, help = 'output submit.py nohup command')
p.add_argument('--zipout', action='store_true', default = False, help = 'add zipout flag to submit.py call')
p.add_argument('-P', '--pipecmds', help = 'pipe commands')
p.add_argument('-a', '--allargs', action='store_true', default = False, help = 'encode all arguments in output name')
p.add_argument('--run', action='store_true', default = False, help = 'run ms command')
p.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
p.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
#p.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)
#pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if args.macs:
	if (args.mst or args.msr or args.msargs or args.mrca):
		error('cannot combine ms and MaCS arguments')
	if args.suffix == 'ms':
		args.suffix = 'macs'

if args.recfile or args.chrmap:
	args.macs = True

	if args.chrmap:
		cstart = -1
		info('reading recombination map file %s' % args.chrmap)
		if args.chrmap.endswith('.gz'):
			tok = [line.split() for line in gzip.open(args.chrmap)]
		else:
			tok = [line.split() for line in open(args.chrmap)]
		chrpos = [int(float(x[0])) for x in tok[1:]]
		recrate = [x[1] for x in tok[1:]]
		args.seqlen = chrpos[-1] - chrpos[0] + 1
		args.recfile = '.'.join([args.prefix, args.suffix, 'recfile'])
		info('map length %d' % args.seqlen)
		info('writing recombination rate file %s' % args.recfile)
		recfile = open(args.recfile, 'w')
		for i in range(len(chrpos) - 1):
			recfile.write('%.9f\t%.9f\t%s\n' % (float(chrpos[i] - chrpos[0]) / args.seqlen, float(chrpos[i + 1] - chrpos[0]) / args.seqlen, recrate[i]))
		recfile.close()

encmd = []
if args.macs and args.nreps:
	encmd.append('-i %d' % (args.nreps))
if not args.mst:
	args.mst = args.mugen * 4 * args.N0 * args.seqlen
	args.macst = args.mugen * 4 * args.N0
	if args.macs:
		encmd.append('-t %s' % (fnum(args.macst)))
if args.npops > 1:
	popsize = args.nsamps/args.npops
	if not popsize * args.npops == args.nsamps:
		error('only equal pop sizes supported; nsamps %d not a multiple of npops %d' % (args.nsamps, args.npops))
	encmd.append('-I %d %s' % (args.npops, ' '.join([str(popsize)] * args.npops)))
if args.recfile:
	encmd.append('-R %s' % (args.recfile))
else:
	if args.msr:
		encmd.append('-r %s %d' % (fnum(args.msr), args.seqlen))
	elif args.recgen > 0.0:
		args.msr = args.recgen * 4 * args.N0 * args.seqlen
		args.macsr = args.recgen * 4 * args.N0
		if args.macs:
			encmd.append('-r %s' % (fnum(args.macsr)))
		else:
			encmd.append('-r %s %d' % (fnum(args.msr), args.seqlen))
if args.trees:
	encmd.append('-T')
if args.mrca:
	encmd.append('-L')
if args.eN:
	for x in args.eN.split():
	#	print(x)
		tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
		Nev = float(eval(x.split(',')[1])) / args.N0
		encmd.append('-eN %s %s' % (fnum(tev), fnum(Nev)))
if args.en:
	for x in args.en.split():
		tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
		pnum = int(x.split(',')[1])
		Nev = float(eval(x.split(',')[2])) / args.N0
		encmd.append('-en %s %d %s' % (fnum(tev), pnum, fnum(Nev)))
if args.ej:
	for x in args.ej.split():
		tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
		p1, p2 = x.split(',')[1:]
		encmd.append('-ej %s %s %s' % (fnum(tev), p1, p2))
if args.em:
	for x in args.em.split():
		tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
		p1, p2 = x.split(',')[1:3]
		mig = eval(x.split(',')[3]) * 4 * args.N0
		encmd.append('-em %s %s %s %s' % (fnum(tev), p1, p2, fnum(mig)))

if args.unscale_string:
	if not args.unscale_string.endswith('.mssim'):
		error('--unscale_string currently only works with .mssim files')
	print('Using mugen = %e per bp per generation, tgen = %.1f y:' % (args.mugen, args.tgen))
	msargs = {}
#	for tok in (x.split('_') for x in simname.split('-')[1:]):
	for i, x in enumerate(args.unscale_string.replace('.mssim', '').split('-')):
		if i == 0 and x.isdigit():
			print('%s samples' % x)
			continue
		if i == 1 and x.isdigit():
			print('%s reps; each with:' % x)
			continue
		argl = x.split('_')[0]
		if argl in msargs:
			msargs[argl].append(x)
		else:
			msargs[argl] = [x]
#	print(msargs)
	if 'r' in msargs:
		seqlen = int(msargs['r'][0].split('_')[2])
#	print('\tseqlen = %d kbp' % (seqlen/1e3))
	print('\tseqlen = %s bp' % "{:,}".format(seqlen))
	if 't' in msargs:
		args.mst = float(msargs['t'][0].split('_')[1])
		args.N0 = args.mst / (4 * args.mugen * seqlen)
#		print('\tN0 = %d' % N0)
		print('\tN0 = %s' % "{:,}".format(int(args.N0)))
	if 'r' in msargs:
		rho = float(msargs['r'][0].split('_')[1]) / (4 * args.N0 * seqlen)
		print('\trec rate = %e per bp per generation' % rho)
	if 'I' in msargs:
		popn = msargs['I'][0].split('_')
		print('\t%s populations; numbers of samples: %s' % (popn[1], ', '.join(popn[2:])))
	if 'en' in msargs:
		for arg in msargs['en']:
			tj = float(arg.split('_')[1]) * args.N0 * args.tgen
			pop = arg.split('_')[2]
			siz = float(arg.split('_')[3]) * args.N0
			print('\tAt %d yrs: population %s size %d' % (tj, pop, siz))
	if 'ej' in msargs:
		for arg in msargs['ej']:
			tj = float(arg.split('_')[1]) * args.N0 * args.tgen
			pop = arg.split('_')[2:]
			print('\tAt %d yrs: population %s joins population %s' % (tj, pop[0], pop[1]))
	if 'em' in msargs:
		for arg in msargs['em']:
			tj = float(arg.split('_')[1]) * args.N0 * args.tgen
			pop = arg.split('_')[2:4]
			mig = float(arg.split('_')[4]) / (4 * args.N0)
			print('\tAt %d yrs: migration fraction of population %s from population %s set to %e' % (tj, pop[0], pop[1], mig))
	sys.exit(0)

if args.macs:
	if not args.macsargs:
		args.macsargs = '%d %d' % (args.nsamps, args.seqlen)

	outargs = ' '.join([args.macsargs, ' '.join(encmd)])
	cmd = ' '.join(['macs', outargs])
else:
	if not args.msargs:
		args.msargs = '%d %d -t %s' % (args.nsamps, args.nreps, fnum(args.mst))

	outargs = ' '.join([args.msargs, ' '.join(encmd)])
	cmd = ' '.join(['ms', outargs])

outname = args.prefix
if args.encode_pars:
	if outname:
		outname += '.'
	if args.allargs or args.macs:
		outname += '_'.join(outargs.split()).replace('_-', '-')
	else:
		outname += '_'.join(outargs.split()[2:]).replace('_-', '-')
if args.suffix:
	outname += '.' + args.suffix

if args.pipecmds:
	cmd += ' | ' + args.pipecmds

if args.bsub:
	cmd = 'submit.py bsub \'' + cmd + '\' -M 2 -o %s -v' % outname
	if args.zipout:
		cmd += ' --zipout'
elif args.batch:
	cmd = 'submit.py nohup \'' + cmd + '\' -o %s -v' % outname
	if args.zipout:
		cmd += ' --zipout'
elif args.outfile:
	boutname = outname + '.bout'
	cmd += ' > %s 2> %s' % (outname, boutname)
#redirect = '>'
#if pipecmds:
#	redirect = ' '.join(['|', pipecmds, redirect])

#cmd = ' '.join(['ms', msargs, redirect, outname])
#print('ms %s %s %s' % (msargs, redirect, outname))

if args.run:
#	info('args \'%s\'' % (' '.join(sys.argv[1:])))
	info('running \'%s\'' % (cmd))
	subcall(cmd, args.sim, wait = True)
else:
	print(cmd)

