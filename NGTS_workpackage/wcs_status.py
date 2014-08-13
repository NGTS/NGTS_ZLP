import fitsio

HEADER_KEY = 'wcscompl'

def wcs_succeeded(fname):
    '''
    Return true if the wcs has succeeded for a file
    '''
    header = fitsio.read_header(fname)
    return HEADER_KEY not in header or header[HEADER_KEY]

def set_wcs_status(fname, succeeded):
    '''
    Set the corresponding header key to `succeeded`
    '''
    with fitsio.FITS(fname, 'rw') as outfile:
        outfile[0].write_key(HEADER_KEY, succeeded, comment='WCS succeeded?')
