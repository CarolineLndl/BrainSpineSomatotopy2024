#!/usr/bin/env python
# -*- coding: utf-8	-*-

import os, sys , argparse
import numpy as	np
import pydicom
from shutil import copyfile, rmtree, move

def	get_arguments():
	parser	= argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="",
		epilog="""
		Sort DICOM Files	
		Input: Folder
		""")
	
	parser.add_argument(
		"-d", "--idir",
		required=True, nargs="+",
		help="Folder to be sorted",
		)
	
	parser.add_argument(
		"-r", "--rmRaw",
		required=False, nargs="+",
		help="Remove	raw	folder",
		)	 
	
	parser.add_argument(
		"-k", "--keepName",
		required=False, nargs="+",
		help="Keep old name",
		)
	
	parser.add_argument(
		"-o", "--odir",
		required=True, nargs="+",
		help="Output	folder - if	doesnt exist it	will be	created",
		)
	
	parser.add_argument(
		"-v", "--verbose",
		required=False, nargs="+",
		help="Verbose to	get	more information about whats going on",
		)

	args =	parser.parse_args()
	if	len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	else:
		return args


class sortDCM(object):
	"""
	"""
	def __init__(
		self, idir, odir, rmRaw=False, keepName=False,
		verbose=False, log_level="INFO"):
	
		self.inputFolder =	idir[0]
		if not os.path.exists(self.inputFolder):
			print('Input dir does not exit: {}'.format(self.inputFolder))
			sys.exit()

		if not os.path.exists(odir[0]):
			os.mkdir(odir[0])

		self.outputFolder = odir[0]
		self.keepName = keepName
		self.verbose = verbose
		self.rmRaw = rmRaw

	def run(self):
		onlyFiles = [f for f in os.listdir(self.inputFolder) if os.path.isfile(os.path.join(self.inputFolder, f))]
		onlyFiles.sort()
		meSeq = []

		for nFile in onlyFiles:

			iFile = os.path.join(self.inputFolder, nFile)
			ds = pydicom.dcmread(iFile)	# Read File
			
			seriesFolder = os.path.join(self.outputFolder, '{:02d}'.format(ds.SeriesNumber) + '-' + str(ds.SeriesDescription.replace(' ','_')))
			
			if not os.path.exists(seriesFolder): # Create Serie Directory
				os.mkdir(seriesFolder)
				if self.verbose:
					print('Create new series of dicoms: {}'.format(seriesFolder))

			if not self.keepName: # Change Name
				newName	= os.path.join(seriesFolder, str(ds.PatientName) + '-' + '{:03d}'.format(ds.InstanceNumber) + '.dcm')
			else:
				newName = os.path.join(seriesFolder, nFile)

			if not os.path.exists(newName):
				copyfile(iFile,	newName)
				if self.verbose:
					print('Copy file {} to {}'.format(nFile, newName))
			else:
				meSeq.append(seriesFolder)
				copyfile(iFile, newName.replace('.dcm','_' + str(ds.EchoNumbers) + '.dcm'))
				print('ERROR: {} already exists with this new name {}'.format(iFile, newName))

				meSeq = list(set(meSeq))

				for nSeries in meSeq:
					onlyFiles = [f for f in os.listdir(nSeries) if os.path.isfile(os.path.join(nSeries, f))]
					for nFile in onlyFiles:
						iFile = os.path.join(nSeries, nFile)
						ds = pydicom.dcmread(iFile)
						seriesFolder = os.path.join(nSeries, 'echo_' + str(ds.EchoNumbers))

						if not os.path.exists(seriesFolder):
							os.mkdir(seriesFolder)

						newName = os.path.join(seriesFolder, str(ds.PatientName) + '-' + '{:03d}'.format(ds.InstanceNumber) + '.dcm')
						move(iFile, newName)

		# Remove input Folder
		if self.rmRaw:
			rmtree(self.inputFolder)
			if self.verbose:
				print('Remove RAW folder: {}'.format(self.inputFolder))


def	main():
	"""Let's go"""
	args =	get_arguments()
	app = sortDCM(**vars(args))
	return	app.run()

if __name__	== '__main__':
	sys.exit(main())
