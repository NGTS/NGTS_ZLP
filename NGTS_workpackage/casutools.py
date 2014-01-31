import subprocess as sp

def find_imstack():
    names = ['casu_imstack', 'imstack']
    for name in names:
        try:
            sp.check_call(['which', name], stderr=sp.PIPE, stdout=sp.PIPE)
        except sp.CalledProcessError:
            pass
        else:
            return name

def construct_filelist_argument(filelist):
    return '@{}'.format(filelist)

def run_command(args, verbose=False):
    str_cmd = map(str, cmd)

    if verbose:
        print str_cmd

    sp.check_call(str_cmd)

def imstack(filelist, confidence_map, outstack='outstack.fits', outconf='outconf.fits',
        verbose=False, casu_verbose=False):
    """
    Runs the casu task `imstack`
    """
    catalogues = ""

    cmd = [find_imstack(),
            construct_filelist_argument(filelist),
            confidence_map,
            catalogues,
            outstack,
            outconf]

    if casu_verbose:
        cmd.append('--verbose')

    run_command(cmd, verbose=verbose)

