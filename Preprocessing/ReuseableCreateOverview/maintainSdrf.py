#!/usr/bin/env python

# Title			: 	maintainSdrf.py
# Created by	: 	Marissa Dubbelaar
# Created on	: 	24-02-2016
# Purpose		:	1) Makes sure that the SDRF files are downloaded.
#					2) Reads the content of the SDRF files.

import urllib
import os

class MaintainSdrf():
	def __init__(self, pathway):
		'''
			A directory is made where the SDRF files are stored (temporarily).
			The rest of the function can be used to be called within the main.
		'''
		self.pathway = pathway
		self.createSdrfDir()

		self.allInfo = []
		self.material = []
		self.tissue = []
		self.celltype = []
		self.age = []
		self.strain = []
		self.sampleInfo = []

	def createSdrfDir(self):
		'''
			Check whether the given pathway exsists, if not it will be made.	
		'''
		if not os.path.exists(self.pathway):
			os.makedirs(self.pathway)
		os.chdir(self.pathway)

	def downloadSdrf(self, EGEODNR):
		'''
			The sdrf files are obtained for the given EGEOD number.
		'''
		urllib.urlretrieve("https://www.ebi.ac.uk/arrayexpress/files/"+EGEODNR+"/"+EGEODNR+".sdrf.txt", EGEODNR+".sdrf.txt")
		with open(EGEODNR+".sdrf.txt") as f:
			lines = f.readlines()
			# Some studies do not contain a .sdrf.txt, they have an .hyb.sdrf.txt file instead.
			# So when the data of the .sdrf.txt is not available, the .hyb.sdrf.txt file is downloaded and the .sdrf.txt file is removed.
			if "not publicly available" in lines[0]:
				os.remove(EGEODNR+".sdrf.txt")
				urllib.urlretrieve("https://www.ebi.ac.uk/arrayexpress/files/"+EGEODNR+"/"+EGEODNR+".hyb.sdrf.txt", EGEODNR+".sdrf.txt")
		# Information from this file is retreived and returned.
		self.obtainSdrfInfo(EGEODNR+".sdrf.txt")
		return self.sampleInfo

	def obtainSdrfInfo(self, EGEODfile):
		'''
			The function obtainSdrfInfo obtains information about the samples such as:
				- Material
				- Tissue type
				- Cell type
				- Age
				- Strain (mutation)
		'''

		self.sampleInfo = []
		# File is opened and the content is read.
		f = open(EGEODfile)
		lines = f.readlines()
		# Files with less than 100 lines (100 samples) are used further.
		if len(lines)-1 < 100:
			for line in lines:
				# The header is used to determine the index of the information.
				if line.startswith("Source "):
					# The header is tab separated.
					line = line.lower().split("\t")
					# The information defined by "library_strategy" is used as information about the material.
					self.material = [i for i, x in enumerate(line) if "library_strategy" in x]
					# The information defined by "comment Sample_source_name" or "characteristics organism part" is used as information about the tissue.
					self.tissue = [i for i, x in enumerate(line) if x.startswith("comment [Sample_source_name") or x.startswith("characteristics [organism part")]
					# The information defined by "cell type" or "celltype" is used as information about the cell type.
					self.celltype = [i for i, x in enumerate(line) if "cell type" in x or "celltype" in x]
					# The information defined by "characteristics age" or "age" is used as information about the age.
					self.age = [i for i, x in enumerate(line) if x.startswith("characteristics [age") or "age" in x]
					# The information defined by "strain", "genetic background" or "disease" is used as information about the strain.
					self.strain = [i for i, x in enumerate(line) if "strain" in x or "genetic background" in x or "disease" in x]
				else:
					try:
						self.allInfo = []
						line = line.split("\t")
						# It is determined if the samples are RNA sequencing data.
						if line[self.material[0]].lower() == "rna-seq":
							# The identifiers are obtained (SRA, SRX and SRS numbers).
							self.identifier = [i for i, x in enumerate(line) if x.startswith("SRR") or x.startswith("SRX") or x.startswith("SRS")]
							# The information is obtained from each rows with the use of the indexes that were obtained.
							for spec in [self.celltype, self.age, self.strain, self.tissue]:
								if spec == []:
									spec = [0]
									line[spec[0]] = " "	
								self.allInfo.insert(len(self.allInfo), line[spec[0]])
							self.allInfo.insert(len(self.allInfo), line[self.identifier[0]])
							# Information is saved into the self.sampleInfo list so it can be returned by the function downloadSdrf.
							self.sampleInfo.append(self.allInfo)
					except:
						continue

	def removeSdrf(self):
		'''Removes the folder where all of the sdrf files are stored.'''
		shutil.rmtree(self.pathway)