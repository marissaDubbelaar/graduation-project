#!/usr/bin/env python

# Title			: 	samplesSheetCreater.py
# Created by	: 	Marissa Dubbelaar
# Created on	: 	26-05-2016
# Purpose		:	Create a samples sheet which can be used within the RNA sequencing pipeline of Calculon.

import os 

class CreateSampleSheet():
	def __init__(self, pathwayFastQ, outputFile):
		'''
			As input the pathway to the FastQ files and the output file are given.
			The init calls the necessary functions to obtain the names of the files and 
			to create and close the file where the information is written.
		'''
		self.pathwayFastQ = pathwayFastQ
		self.outputFile = outputFile
		# A count is made to write the index of the file.
		self.count = 1
		# The file is made
		self.createCSV()
		# The names of the files within the pathwayFastQ are obtained
		self.readFiles()
		# CSV files is closed again.
		self.closeCSV()

	def readFiles(self):
		'''
			This function obtains all of the known fastq files within the given pathway.
		'''
		# The root and file names are obtained that are available in the self.pathwayFastQ
		for root, directories, files in os.walk(self.pathwayFastQ):
			# The root is splitted to make sure that the ensemble id is obtained for further use.
			root = root.split("/")
			# For every filename in the pathway.
			for f in files:
				# Make sure that the files end with .fastq or .fastq.gz
				if f.endswith(".fastq") or f.endswith(".fastq.gz"):
					# The file name is splitted on the dot, resulting in the SRA number and the extention.
					f= f.split(".")
					# The SRA number and the ensemble id is used with the function "writeCSV"
					self.writeCSV(f[0], root[-2])

	def createCSV(self):
		'''
			Creates a csv where all of the necessary information for the RNA sequencing pipeline needs to be stored.
		'''
		self.csv = open(self.outputFile, 'w')
		self.csv.write("internalId,project,sampleName,reads1FqGz,reads2FqGz,sortedBamFile\n")

	def closeCSV(self):
		'''
			The CSV file is closed.
		'''
		self.csv.close()

	def writeCSV(self, SRANumber, EGEODNumber):
		'''
			All of the necessary information will be written into the file with the use of the function writeCSV.
		'''
		# If the fastq file is part of a pairwise sequencing.
		if "_" in SRANumber and SRANumber.endswith("_1"):
			# The SRNA number splitted and used to be written into the file and the column "reads2FqGz" is filled in as well.
			SRANumber, addition = SRANumber.split("_")
			self.csv.write(str(self.count)+",EGEOD_"+EGEODNumber+","+SRANumber+",/groups/umcg-gcc/tmp04/umcg-mdubbelaar/GOAD/"
				+EGEODNumber+"/"+SRANumber+"_1.fastq,/groups/umcg-gcc/tmp04/umcg-mdubbelaar/GOAD/"
				+EGEODNumber+"/"+SRANumber+"_2.fastq,/groups/umcg-gcc/tmp04/umcg-mdubbelaar/tryKallisto/miceData/EGEOD_"
				+EGEODNumber+"/run01/results/sortedBam/"+SRANumber+"_"+str(self.count)+".bam\n")
			self.count += 1
		# If the fastq file is sequenced singlewise
		elif "_" not in SRANumber:
			# The SRA number is used and the column "reads2FqGz" remains empty.
			self.csv.write(str(self.count)+",EGEOD_"+EGEODNumber+","+SRANumber+",/groups/umcg-gcc/tmp04/umcg-mdubbelaar/GOAD/"
				+EGEODNumber+"/"+SRANumber+".fastq,,/groups/umcg-gcc/tmp04/umcg-mdubbelaar/tryKallisto/miceData/EGEOD_"
				+EGEODNumber+"/run01/results/sortedBam/"+SRANumber+"_"+str(self.count)+".bam\n")
			self.count += 1

css = CreateSampleSheet("/Volumes/Elements/School/Afstuderen/GOAD3/Data/53737/", "/Users/marissa/Desktop/samplesSheet53737.csv")
