# -*- coding: utf-8 -*-
'''
Implements the used casutools tasks in python scripts

Functions require the associated tools to be available on the
system.
'''


import subprocess as sp
from threading import Lock

def find_imstack():
    '''
    Function to find imstack, as it has been renamed on ngtshead
    '''
    names = ['casu_imstack', 'imstack']
    for name in names:
        try:
            sp.check_call(['which', name], stderr=sp.PIPE, stdout=sp.PIPE)
        except sp.CalledProcessError:
            pass
        else:
            return name

def construct_filelist_argument(filelist):
    '''
    Wrapper around constructing a filelist
    '''
    return '@{0}'.format(filelist)


lock = Lock()
def run_command(cmd, verbose=False):
    '''
    Wraps subprocess to run the command
    '''
    str_cmd = map(str, cmd)

    if verbose:
        with lock:
            print ' '.join(str_cmd)

    sp.check_call(str_cmd)

def imstack(filelist, confidence_map, outstack='outstack.fits', outconf='outconf.fits',
        verbose=False, catalogues=''):
    '''
    Runs the casu task `imstack`
    '''
    cmd = [find_imstack(),
            construct_filelist_argument(filelist),
            confidence_map,
            catalogues,
            outstack,
            outconf]

    run_command(cmd, verbose=verbose)


def imcore(input_file, output_table, ipix=2, threshold=2.0, confidence_map='noconf', rcore=2,
        filtfwhm=1, ellfile=False, casu_verbose=False, verbose=False):
    '''
    Runs the casu task `imcore`
    '''
    cmd = ['imcore', input_file, confidence_map, output_table, ipix, threshold,
            '--filtfwhm', filtfwhm,
            '--rcore', rcore]

    if casu_verbose:
        cmd.append('--verbose')

    if not ellfile:
        cmd.append('--noell')

    run_command(cmd, verbose=verbose)

def imcore_list(input_file, listfile, output_file, threshold=2.0, confidence_map='noconf',
        rcore=5, casu_verbose=False, noell=True, verbose=False):
    '''
    Runs the casu task `imcore_list`
    '''
    cmd = ['imcore_list', input_file, confidence_map, listfile, output_file,
            threshold, '--rcore', rcore]

    if noell:
        cmd.append('--noell')

    if casu_verbose:
        cmd.append('--verbose')

    run_command(cmd, verbose=verbose)

def wcsfit(infile, incat, catsrc='viz2mass', site='cds', catpath=False, verbose=False):
    '''
    Runs the casu task `wcsfit`
    '''
    cmd = ['wcsfit', infile, incat, '--catsrc', catsrc, '--site', site]

    if 'local' in catsrc:
      print 'using local catalog!'
      cmd = cmd[:-2]
      if catpath:
        cmd.append('--catpath')
        cmd.append(catpath)
      else:
        print 'no path supplied!'
        quit()

    run_command(cmd, verbose=verbose)
