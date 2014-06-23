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
#import re
import locale
import string
import decimal
import argparse

from aosutils import *

# global defaults
loglevel = logging.WARNING

def fnum(num, sf = 3): # round to 3sf and format compactly
	s = []
	nf = 0
	ppos = -1
	for x in str(num):
#		print((x, s))
		if x == '.':
			ppos = len(s)
			continue
		if nf == 0 and ppos < 0 and x == '0':
			continue
		s.append(x)
		if x != '-':
			nf += 1
		if ppos >= 0 and nf > sf:
			if int(s[-1]) >= 5:
				s[-2] = str(int(s[-2]) + 1)
			s = s[:-1]
			break
	if len(s) == 0:
		s = ['0']
	if ppos >= 0:
		s.insert(ppos, '.')
		if s[0] == '.':
			s.insert(0, '0')
		return(''.join(s).rstrip('0').rstrip('.'))
	else:
		return(''.join(s))


#def usage():
#	print('usage: %s [msargs] OPTIONS' % (os.path.basename(sys.argv[0])))
#	print('usage: %s [-s nsamps] [-n nreps] [-p npops] [-N N0] [-l seqlen] [-u mu_gen [-r rho_gen] | --ms_theta=ms_theta [--ms_rho=ms_rho]] [-g tgen] [--eN=Ne_all_history] [--en=Ne_pop_history] [--ej=merge_history] [--em=migration_history] [--trees] [--mrca] OPTIONS' % (os.path.basename(sys.argv[0])))
#	print('usage: %s --unscale MSCMD [-l seqlen] [-u mu_gen] [-g tgen]' % (os.path.basename(sys.argv[0])))
#	print('''
#	msargs: 'nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <MSOPTIONS>'
#		(see msdoc for more options)
#	ms scaling:
#		ms_theta = 4 * mu_gen * N0 * seqlen (e.g. = 4.8e-4 * seqlen for human)
#		ms_rho = 4 * rec_gen * N0 * seqlen (e.g. = 4e-4 * seqlen for human)
#	Ne_all_history: 'time,Ne ...'
#	Ne_pop_history: 'time,pop,Ne ...'
#	merge_history: 'time,from_pop,to_pop ...'
#	migration_history: 'time,to_pop,from_pop,to_pop_mig_frac ...'
#	''')

os.umask(0002)

p = argparse.ArgumentParser()
p.addargument('MSARGS', nargs = '?', help = 'ms arguments: "nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <ms_options>" (see msdoc for more options)')
p.addargument('-n', '--nreps', type = int, default = 1, 'number of repetitions')
p.addargument('-u', '--mugen', type = float, default = 1.25e-8, 'per-generation mutation rate')
p.addargument('-r', '--recgen', type = float, default = 0.0, 'per-generation recombination rate'))
p.addargument('-l', '--seqlen', type = int, default = 10e3, 'sequence length simulated'))
p.addargument('-P', '--pipecmds', default = '', help = 'pipe commands')
p.addargument('-a', '--allargs', action='store_true', default = False, help = 'encode all arguments in output name')
p.addargument('-p', '--npops', type = int, default = 0, 'number of populations'))
p.addargument('-N', '--N0', type = int, default = 1e4, 'base effective population size')
p.addargument('-g', '--tgen', type = float, default = 30.0, 'generation time (y)'))
p.addargument('-s', '--nsamps', type = int, default = 4, 'number of samples per rep')
p.addargument('--trees', action='store_true', default = False, help = 'inlude trees in ms output')
p.addargument('--mrca', action='store_true', default = False, help = 'inlude TMRCA in ms output')
p.addargument('--run', action='store_true', default = False, help = 'run ms command')
p.addargument('--eN', default='', help='global Ne history: "time,Ne ..."' )
p.addargument('--en', default='', help='population Ne history: "time,pop,Ne ..."' )
p.addargument('--ej', default='', help='population merge_history: "time,from_pop,to_pop ..."')
p.addargument('--em', default='', help='migration_history: "time,to_pop,from_pop,to_pop_mig_frac ..."')
p.addargument('--outfile', action='store_true', default = False, help = 'redirect output to OUTFILE')
p.addargument('--bsub', action='store_true', default = False, help = 'output bsub.py command')
p.addargument('--prefix', default='', help='output prefix')
p.addargument('--suffix', default='mssim', help='output suffix')
p.addargument('--mst', type = float, default = 0.0, 'ms theta value; otherwise ms_theta = 4 * MUGEN * N0 * SEQLEN (e.g. = 4.8e-4 * SEQLEN for human)')
p.addargument('--msr', type = float, default = 0.0, 'ms rho value; otherwise ms_rho = 4 * RECGEN * N0 * SEQLEN (e.g. = 4.0e-4 * SEQLEN for human)')
p.addargument('--unscale_string', default='', help='print unscaled parameters for simulation output name UNSCALE_STRING')
p.add_argument('--sim', action='store_true', default = False, help = 'dry run')
p.add_argument('-v', '--verbose', action='store_true', default = False)#, help = 'dry run')
p.add_argument('--debug', action='store_true', default = False, help=argparse.SUPPRESS)
pp.add_argument('--bsim', action='store_true', default = False, help=argparse.SUPPRESS)
#pp.add_argument('--replace', action='store_true', default = False, help = 'replace existing files')

args = p.parse_args()

loglevel = logging.WARNING
if args.verbose:
	loglevel = logging.INFO
if args.debug:
	loglevel = logging.DEBUG
logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

encmd = []
if args.mst > 0.0:
	args.mst = args.mugen * 4 * args.N0 * args.seqlen
if args.msr > 0.0:
	encmd.append('-r %s %d' % (fnum(args.msr), seqlen))
if args.npops:
	popsize = args.nsamps/args.npops
	if not popsize * args.npops == args.nsamps:
		error('only equal pop sizes supported; nsamps %d not a multiple of npops %d' % (args.nsamps, args.npops))
	encmd.append('-I %d %s' % (args.npops, ' '.join([str(popsize)] * args.npops)))
if args.recgen > 0.0:
	args.msr = args.recgen * 4 * args.N0 * seqlen
	encmd.append('-r %s %d' % (fnum(args.msr), seqlen))
if args.trees:
	encmd.append('-T')
if args.mrca:
	encmd.append('-L')
for x in args.eN.split():
#	print(x)
	tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
	Nev = float(eval(x.split(',')[1])) / args.N0
	encmd.append('-eN %s %s' % (fnum(tev), fnum(Nev)))
for x in args.en.split():
	tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
	pnum = int(x.split(',')[1])
	Nev = float(eval(x.split(',')[2])) / args.N0
	encmd.append('-en %s %d %s' % (fnum(tev), pnum, fnum(Nev)))
for x in args.ej.split():
	tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
	p1, p2 = x.split(',')[1:]
	encmd.append('-ej %s %s %s' % (fnum(tev), p1, p2))
for x in args.em.split():
	tev = eval(x.split(',')[0]) / (4 * args.N0 * args.tgen)
	p1, p2 = x.split(',')[1:3]
	mig = eval(x.split(',')[3]) * 4 * args.N0
	encmd.append('-em %s %s %s %s' % (fnum(tev), p1, p2, fnum(mig)))

if args.unscale_string:
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


if not args.msargs:
	args.msargs = '%d %d -t %s' % (args.nsamps, nreps, fnum(args.mst))

msargs = ' '.join([args.msargs, ' '.join(encmd)])

if args.allargs:
	outname = '_'.join([args.pref] + msargs.split()).replace('_-', '-')
else:
	outname = '_'.join([args.pref] + msargs.split()[2:]).replace('_-', '-')
if args.suffix:
	outname += '.' + args.suffix

cmd = ' '.join(['ms', msargs])

if args.pipecmds:
	cmd += ' | ' + args.pipecmds

if args.bsub:
	cmd = 'bsub.py \'' + cmd + '\' -M 2 -o %s' % outname
elif args.outfile:
	cmd += ' > ' + outname
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

