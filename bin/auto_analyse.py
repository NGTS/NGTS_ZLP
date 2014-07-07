#!/usr/bin/env python

import os, time
import glob

path_to_watch = "."

label = '*.fits'

filelist = 'temp_filelist'

confmap = '/ngts/pipedev/InputCatalogue/srw_confidence.fits'
catfile = '/ngts/pipedev/InputCatalogue/output/SimonTest2/SimonTest2_dither_NG190335+491133/catfile.fits'

before = dict([(f,None) for f in glob.glob(label)])

while 1:
	time.sleep(1)
	after = dict([(f,None) for f in glob.glob(label)])
	added = [os.path.abspath(f) for f in after if not f in before]
	if len(added) > 0:
		with open(filelist,'w') as out:
			for f in added:
				out.write(f + '\n')
		command = 'ZLP_app_photom.py --confmap '+confmap+' --catfile '+catfile+' --filelist '+os.path.abspath(filelist)
		print command
		os.system(command)
	before = after
