import fitsio

__all__ = ['wcs_succeeded', 'set_wcs_status', 'filter_wcs_successes']

HEADER_KEY = 'wcscompl'

def wcs_succeeded(fname, hdu=0):
    '''
    Return true if the wcs has succeeded for a file
    '''
    header = fitsio.read_header(fname, ext=hdu)
    return HEADER_KEY not in header or header[HEADER_KEY]

def set_wcs_status(fname, succeeded, hdu=0):
    '''
    Set the corresponding header key to `succeeded`
    '''
    with fitsio.FITS(fname, 'rw') as outfile:
        outfile[hdu].write_key(HEADER_KEY, succeeded, comment='WCS succeeded?')
