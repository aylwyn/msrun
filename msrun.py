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

from bputil import *

# global defaults
loglevel = logging.WARNING

sim = 0
allargs = True
run = False
pref = ''
pipecmds = ''
suf = 'mssim'
outfile = False
bsubout = False
unscale = False
simname = ''

os.umask(0002)
encmd = []
N0 = 1e4
tgen = 30.0
nsamps = 4
nreps = 1
seqlen = 10e3
mugen = 1.25e-8

mst = mugen * 4 * N0 * seqlen

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


def usage():
	print('usage: %s [msargs] OPTIONS' % (os.path.basename(sys.argv[0])))
	print('usage: %s [-s nsamps] [-n nreps] [-p npops] [-N N0] [-l seqlen] [-u mu_gen [-r rho_gen] | --ms_theta=ms_theta [--ms_rho=ms_rho]] [-g tgen] [--eN=Ne_all_history] [--en=Ne_pop_history] [--ej=merge_history] [--em=migration_history] [--trees] [--mrca] OPTIONS' % (os.path.basename(sys.argv[0])))
	print('usage: %s --unscale MSCMD [-l seqlen] [-u mu_gen] [-g tgen]' % (os.path.basename(sys.argv[0])))
	print('''
	msargs: 'nsamps nreps -t mst [-r msr seqlen] [-I npops pop1_nsamps [pop2_nsamps ...]] <MSOPTIONS>'
		(see msdoc for more options)
	ms scaling:
		ms_theta = 4 * mu_gen * N0 * seqlen (e.g. = 4.8e-4 * seqlen for human)
		ms_rho = 4 * rec_gen * N0 * seqlen (e.g. = 4e-4 * seqlen for human)
	Ne_all_history: 'time,Ne ...'
	Ne_pop_history: 'time,pop,Ne ...'
	merge_history: 'time,from_pop,to_pop ...'
	migration_history: 'time,to_pop,from_pop,to_pop_mig_frac ...'
	OPTIONS: [-P pipe_cmds] [--prefix=prefix] [--run] [--suffix=suffix] [--sim] [--outf | --bsub]
	''')
	sys.exit(2)

try:
	opts, args = getopt.gnu_getopt(sys.argv[1:], 'n:u:r:l:P:ap:N:vg:h:s:', ['trees', 'mrca', 'run', 'eN=', 'en=', 'ej=', 'em=', 'sim', 'debug', 'outf', 'bsub', 'prefix=', 'suffix=', 'mst=', 'msr=', 'unscale'])
except getopt.GetoptError:
	usage()
	sys.exit(2)
for (oflag, oarg) in opts:
	if oflag == '--sim':
		sim = 2
	if oflag == '--simv':
		sim = 1
		loglevel = logging.INFO
	elif oflag == '-a':
		allargs = True
	elif oflag == '--suffix':
		suf = oarg
	elif oflag == '-P':
		pipecmds = oarg
	elif oflag == '--prefix':
		pref = oarg
	elif oflag == '-s':
		nsamps = int(eval(oarg))
	elif oflag == '-p':
		npops = int(eval(oarg))
		popsize = nsamps/npops
		if not popsize * npops == nsamps:
			error('only equal pop sizes supported; nsamps %d not a multiple of npops %d' % (nsamps, npops))
		encmd.append('-I %d %s' % (npops, ' '.join([str(popsize)] * npops)))
	elif oflag == '-n':
		nreps = int(eval(oarg))
	elif oflag == '-N':
		N0 = int(eval(oarg))
	elif oflag == '-l':
		seqlen = int(eval(oarg))
	elif oflag == '-u':
		mst = float(eval(oarg)) * 4 * N0 * seqlen
	elif oflag == '-r':
		msr = float(eval(oarg)) * 4 * N0 * seqlen
		encmd.append('-r %s %d' % (fnum(msr), seqlen))
	elif oflag == '-g':
		tgen = float(oarg)
	elif oflag == '--trees':
		encmd.append('-T')
	elif oflag == '--mrca':
		encmd.append('-L')
	elif oflag == '--mst':
		mst = float(oarg)
	elif oflag == '--msr':
		msr = float(oarg)
		encmd.append('-r %s %d' % (fnum(msr), seqlen))
	elif oflag == '--eN':
		for x in oarg.split():
#			print(x)
			tev = eval(x.split(',')[0]) / (4 * N0 * tgen)
			Nev = float(eval(x.split(',')[1])) / N0
			encmd.append('-eN %s %s' % (fnum(tev), fnum(Nev)))
	elif oflag == '--en':
		for x in oarg.split():
			tev = eval(x.split(',')[0]) / (4 * N0 * tgen)
			pnum = int(x.split(',')[1])
			Nev = float(eval(x.split(',')[2])) / N0
			encmd.append('-en %s %d %s' % (fnum(tev), pnum, fnum(Nev)))
	elif oflag == '--ej':
		for x in oarg.split():
			tev = eval(x.split(',')[0]) / (4 * N0 * tgen)
			p1, p2 = x.split(',')[1:]
			encmd.append('-ej %s %s %s' % (fnum(tev), p1, p2))
	elif oflag == '--em':
		for x in oarg.split():
			tev = eval(x.split(',')[0]) / (4 * N0 * tgen)
			p1, p2 = x.split(',')[1:3]
			mig = eval(x.split(',')[3]) * 4 * N0
			encmd.append('-em %s %s %s %s' % (fnum(tev), p1, p2, fnum(mig)))
	elif oflag == '-v':
		loglevel = logging.INFO
	elif oflag == '--debug':
		loglevel = logging.DEBUG
	elif oflag == '--run':
		run = True
	elif oflag == '--outf':
		outfile = True
	elif oflag == '--bsub':
		bsubout = True
	elif oflag == '--unscale':
		unscale = True

logging.basicConfig(format = '%(module)s:%(lineno)d:%(levelname)s: %(message)s', level = loglevel)

if len(args) >= 1:
	if unscale:
		simname = args[0]
	else:
		cargs = args[0]
else:
	cargs = '%d %d -t %s' % (nsamps, nreps, fnum(mst))

if unscale:
	print('Using mugen = %e per bp per generation, tgen = %.1f y:' % (mugen, tgen))
	msargs = {}
#	for tok in (x.split('_') for x in simname.split('-')[1:]):
	for i, x in enumerate(simname.replace('.mssim', '').split('-')):
		if i == 0 and x.isdigit():
			print('%s samples' % x)
			continue
		if i == 1 and x.isdigit():
			print('%s reps' % x)
			continue
		argl = x.split('_')[0]
		if argl in msargs:
			msargs[argl].append(x)
		else:
			msargs[argl] = [x]
#	print(msargs)
	if 'r' in msargs:
		seqlen = int(msargs['r'][0].split('_')[2])
	print('seqlen = %d bp' % (seqlen))
	if 't' in msargs:
		mst = float(msargs['t'][0].split('_')[1])
		N0 = mst / (4 * mugen * seqlen)
		print('N0 = %d' % N0)
	if 'r' in msargs:
		rho = float(msargs['r'][0].split('_')[1]) / (4 * N0 * seqlen)
		print('rec rate = %e per bp per generation' % rho)
	if 'I' in msargs:
		popn = msargs['I'][0].split('_')
		print('%s populations; numbers of samples: %s' % (popn[1], ', '.join(popn[2:])))
	if 'en' in msargs:
		for arg in msargs['en']:
			tj = float(arg.split('_')[1]) * N0 * tgen
			pop = arg.split('_')[2]
			siz = float(arg.split('_')[3]) * N0
			print('At %d yrs: population %s size %d' % (tj, pop, siz))
	if 'ej' in msargs:
		for arg in msargs['ej']:
			tj = float(arg.split('_')[1]) * N0 * tgen
			pop = arg.split('_')[2:]
			print('At %d yrs: population %s joins population %s' % (tj, pop[0], pop[1]))
	if 'em' in msargs:
		for arg in msargs['em']:
			tj = float(arg.split('_')[1]) * N0 * tgen
			pop = arg.split('_')[2:4]
			mig = float(arg.split('_')[4]) / (4 * N0)
			print('At %d yrs: migration fraction of population %s from population %s set to %e' % (tj, pop[0], pop[1], mig))
	sys.exit(0)


msargs = ' '.join([cargs, ' '.join(encmd)])

if allargs:
	outname = '_'.join([pref] + msargs.split()).replace('_-', '-')
else:
	outname = '_'.join([pref] + msargs.split()[2:]).replace('_-', '-')
if suf:
	outname += '.' + suf

cmd = ' '.join(['ms', msargs])

if pipecmds:
	cmd += ' | ' + pipecmds

if bsubout:
	cmd = 'bsub.py \'' + cmd + '\' -M 2 -o %s' % outname
elif outfile:
	cmd += ' > ' + outname
#redirect = '>'
#if pipecmds:
#	redirect = ' '.join(['|', pipecmds, redirect])

#cmd = ' '.join(['ms', msargs, redirect, outname])
#print('ms %s %s %s' % (msargs, redirect, outname))

if run:
#	info('args \'%s\'' % (' '.join(sys.argv[1:])))
	info('running \'%s\'' % (cmd))
	subcall(cmd, sim, wait = True)
else:
	print(cmd)

