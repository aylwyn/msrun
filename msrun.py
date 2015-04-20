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
from numpy import log

import aosutils

#p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p = argparse.ArgumentParser()
p.add_argument('-s', '--nsamps', type = int, default = 4, help = 'number of samples per rep (default = 4)')
p.add_argument('-n', '--nreps', type = int, default = 1, help = 'number of repetitions (default = 1)')
p.add_argument('-l', '--seqlen', type = int, default = 10000, help = 'sequence length simulated (default = 10000)')
p.add_argument('-p', '--npops', type = int, default = 1, help = 'number of populations (default = 1)')
p.add_argument('-u', '--mugen', type = float, default = 1.25e-8, help = 'per-generation mutation rate (default = 1.25e-8)')
p.add_argument('-N', '--N0', type = int, default = 10000, help = 'base effective population size (default = 10000)')
p.add_argument('-r', '--recgen', type = float, default = 0.0, help = 'per-generation recombination rate')
p.add_argument('-g', '--tgen', type = float, default = 30.0, help = 'generation time (y) (default = 30.0)')
p.add_argument('-T', '--trees', action='store_true', default = False, help = 'include trees in ms output')
p.add_argument('--mrca', action='store_true', default = False, help = 'include TMRCA in ms output')

p.add_argument('--eN', action='append', default = [], nargs = 2, help='global Ne history', metavar = ('time', 'Ne'))
p.add_argument('--eG', action='append', default = [], nargs = 2, help='global growth history', metavar = ('time', 'alpha'))
p.add_argument('--eNG', action='append', default = [], nargs = 4, help='global growth history in terms of time intervals and Ne change (Note: t2 is older than t1)', metavar = ('time1', 'time2', 'Ne1', 'Ne2'))
p.add_argument('--en', action='append', default = [], nargs = 3, help='population Ne history', metavar = ('time', 'pop', 'Ne'))
p.add_argument('--eg', action='append', default = [], nargs = 3, help='population growth history', metavar = ('time', 'pop', 'alpha'))
p.add_argument('--eng', action='append', default = [], nargs = 5, help='population growth history in terms of time intervals and Ne change (Note: t2 is older than t1)', metavar = ('pop', 'time1', 'time2', 'Ne1', 'Ne2'))
p.add_argument('--ej', action='append', default = [], nargs = 3, help='population merge_history', metavar = ('time', 'from_pop', 'to_pop'))
p.add_argument('--em', action='append', default = [], nargs = 4, help='migration_history', metavar = ('time', 'to_pop', 'from_pop', 'to_pop_mig_frac'))

p.add_argument('--macs', action='store_true', default = False, help = 'output MaCS command (otherwise use ms)')
p.add_argument('--scrm', action='store_true', default = False, help = 'output scrm command (otherwise use ms)')
p.add_argument('--recfile', help='MaCS-formatted recombination rate file (enforces --macs, ignores -r)' )
p.add_argument('--chrmap', help='chromosomal recombination rate file. Sets SEQLEN equal to length of map and writes recfile to PREFIX.macsrec. (Enforces --macs, ignores -r)' )
p.add_argument('--mst', type = float, help = 'ms theta value; overrides --mugen if set, otherwise ms_theta = 4 * MUGEN * N0 * SEQLEN (e.g. = 4.8e-4 * SEQLEN for human)')
p.add_argument('--msr', type = float, help = 'ms rho value; overrides --recgen if set, otherwise ms_rho = 4 * RECGEN * N0 * SEQLEN (e.g. = 4.0e-4 * SEQLEN for human)')
p.add_argument('--msargs', help = 'ms arguments: "nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <ms_options>" (see MS documentation for more options)')
p.add_argument('--macsargs', help = 'MaCS arguments: "nsamps seqlen -i nreps -t macst [-r macsr] [-I npops pop1_nsamps [pop2_nsamps ...]] <macs_options>" (see MaCS documentation for more options)')
p.add_argument('--scrmargs', help = 'scrm arguments: "nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <scrm_options>" (see scrm documentation for more options)')
p.add_argument('--outfile', action='store_true', default = False, help = 'redirect output to file')
p.add_argument('--prefix', default = 'mssim', help='output file prefix')
p.add_argument('--suffix', default='', help='output file suffix')
p.add_argument('--encode_pars', action='store_true', default = False, help = 'encode demographic parameters in output file name (sets --outfile)')
p.add_argument('--encode_all_args', action='store_true', default = False, help = 'encode all arguments in output file name (sets --outfile)')
p.add_argument('--unscale_string', help='print unscaled parameters for simulation output name UNSCALE_STRING')
p.add_argument('--bsub', action='store_true', default = False, help = 'output submit.py bsub command')
p.add_argument('--batch', action='store_true', default = False, help = 'output submit.py nohup command')
p.add_argument('--zipout', action='store_true', default = False, help = 'add zipout flag to submit.py call')
p.add_argument('-P', '--pipecmds', help = 'pipe commands')
p.add_argument('--sigfigs', type=int, default = 4, help = 'sig figs to use for numbers')
p.add_argument('--run', action='store_true', default = False, help = 'run ms command')
p.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
p.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
#p.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)
#pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')

args = p.parse_args()

def fnum(x): return aosutils.fnum(x, args.sigfigs)

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if args.scrm:
	cmdname = 'scrm'
	if args.suffix == 'ms':
		args.suffix = 'scrm'
elif args.macs:
	cmdname = 'macs'
	if (args.mst or args.msr or args.msargs or args.mrca):
		error('cannot combine ms and MaCS arguments')
	if args.suffix == 'ms':
		args.suffix = 'macs'
else:
	cmdname = 'ms'

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
encmd_t = [] # time-ordered commands
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

for x in args.eN:
##	print(x)
#	tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
#	Nev = float(eval(x.split(',')[1])) / args.N0
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	Nev = float(eval(x[1])) / args.N0
	encmd_t.append((tev, '-eN %s %s' % (fnum(tev), fnum(Nev))))
for x in args.eG:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	alpha = float(eval(x[1]))
	encmd_t.append((tev, '-eG %s %s' % (fnum(tev), fnum(alpha))))
for x in args.eNG:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	tev2 = eval(x[1]) / (4 * args.N0 * args.tgen)
	Nev1 = float(eval(x[2])) / args.N0
	Nev2 = float(eval(x[3])) / args.N0
	alpha = -(1 / (tev2 - tev)) * log(Nev2 / Nev1)
	encmd_t.append((tev, '-eN %s %s -eG %s %s' % (fnum(tev), fnum(Nev1), fnum(tev), fnum(alpha))))
	encmd_t.append((tev2, '-eG %s 0' % (fnum(tev2))))
for x in args.en:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	pnum = int(x[1])
	Nev = float(eval(x[2])) / args.N0
	encmd_t.append((tev, '-en %s %d %s' % (fnum(tev), pnum, fnum(Nev))))
for x in args.eg:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	pnum = int(x[1])
	alpha = float(eval(x[2]))
	encmd_t.append((tev, '-eg %s %d %s' % (fnum(tev), pnum, fnum(alpha))))
for x in args.eng:
	pnum = int(x[0])
	tev = eval(x[1]) / (4 * args.N0 * args.tgen)
	tev2 = eval(x[2]) / (4 * args.N0 * args.tgen)
	Nev1 = float(eval(x[3])) / args.N0
	Nev2 = float(eval(x[4])) / args.N0
	alpha = -(1 / (tev2 - tev)) * log(Nev2 / Nev1)
	encmd_t.append((tev, '-en %s %d %s -eg %s %d %s' % (fnum(tev), pnum, fnum(Nev1), fnum(tev), pnum, fnum(alpha))))
	encmd_t.append((tev2, '-eg %s %d 0' % (fnum(tev2), pnum)))
for x in args.ej:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	p1, p2 = x[1:]
	encmd_t.append((tev, '-ej %s %s %s' % (fnum(tev), p1, p2)))
for x in args.em:
	tev = eval(x[0]) / (4 * args.N0 * args.tgen)
	p1, p2 = x[1:3]
	mig = eval(x[3]) * 4 * args.N0
	encmd_t.append((tev, '-em %s %s %s %s' % (fnum(tev), p1, p2, fnum(mig))))

encmd_t.sort()
encmd.append(' '.join([x[1] for x in encmd_t]))

if args.unscale_string:
	if not args.unscale_string.endswith('.ms'):
		error('--unscale_string currently only works with .ms files')
	print('Using mugen = %e per bp per generation, tgen = %.1f y:' % (args.mugen, args.tgen))
	msargs = {}
#	for tok in (x.split('_') for x in simname.split('-')[1:]):
	for i, x in enumerate(args.unscale_string.replace('.ms', '').split('-')):
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
	cmd = ' '.join([cmdname, outargs])
else:
	if not args.msargs:
		args.msargs = '%d %d -t %s' % (args.nsamps, args.nreps, fnum(args.mst))

	outargs = ' '.join([args.msargs, ' '.join(encmd)])
	cmd = ' '.join([cmdname, outargs])

outname = args.prefix
if args.encode_pars:
#	if outname:
#		outname += '.'
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
elif args.outfile or args.encode_pars:
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
	aosutils.subcall(cmd, args.sim, wait = True)
else:
	print(cmd)

