# -*- coding: utf-8 -*-
#!/usr/bin/env python

# Title			: 	collectInfoPerStudy.py
# Created by	: 	Marissa Dubbelaar
# Created on	: 	03-03-2016
# Purpose		:	1) Creates the overview.
#					2) Obtains the global information about the study.
#					3) Obtains information about the samples of the studies.

from obtainAccessions import AccessionObtainer
from connectWebpage import HtmlConnector
from maintainSdrf import MaintainSdrf
from createOverviewXlsx import CreateOverviewXlsx

import re

class CreateOverview():
	def __init__(self):
		'''
			CreateOverview acts as the main for all of the functions.
			It obtaines the accession numbers from the arrayexpress website.
			Using this to obtain the global information of the study.
			The SRDF files of these studies are downloaded to obtain information such as the SRA, SRX or SRS numbers.
			Information about the studies that are defined online are obtained with the SRA, SRX or SRS numbers.
		'''

		self.accessions = []
		# This list contains all of the search terms that need to be checked on arrayexpress.
		self.cellTypes = ["Astrocyte", "Astrocytes", "Astrocytic", "Schwann", "Schwann+Cells", "Microglia", "Microglial", "Oligodendrocyte", "Oligodendrocytes", "Glia",
		"OPC", "Oligodendrocyte+progenitor+cell", "Oligo+progenitor+cell", "Oligo"]

		# The xlsx file is created.
		self.xlsxWriter = CreateOverviewXlsx("/Users/marissa/Desktop/GOAD_GliaOnly_Overview")
		# The knownAccessions are obtained with the use of the function self.xlsxWriter.returnKnownSra().
		self.knownAccessions = self.xlsxWriter.returnKnownSra()

		self.genericInfo = []
		self.abstract = []
		self.sample = []
		self.year = []
		self.author = []

		# Obtains the accession numbers from the arrayexpress page
		self.obtainAcc = AccessionObtainer()
		# Connects to the html page to retreive the content.
		self.connectHtml = HtmlConnector()
		# Obtains the SDRF files.
		self.sdrfObtainer = MaintainSdrf("/Users/marissa/Desktop/TestScripts/ReuseableCreateOverview/SDRF/")

		# Get all EGEOD numbers from the arrayexpress searchpage
		self.getEGEODnr() 
		# Information such as the EGEOD number, the author, year, abstract, title, author and so on are obtained.
		self.obtainNecessaryInfo()

	def getEGEODnr(self):
		'''
			Each searchterm is used to obtain a list of unique accession numbers.
		'''
		for celltype in self.cellTypes:
			# The html content is saved.
			htmlContent = self.connectHtml.obtainHtmlContent("https://www.ebi.ac.uk/arrayexpress/search.html?query="+celltype+"+AND+rna+NOT+mirna&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22sequencing+assay%22&page=1&pagesize=1000")
			# Obtains the E-GEOD number and creates a list of all of the unique found accession numbers on the homepage.
			self.accessions = self.obtainAcc.getUniqueAccessionName(htmlContent)

	def obtainNecessaryInfo(self):
		'''
			Each E-GEOD number is used to connect to their corresponding page.
			Obtaining information as the author and year, 
			And information from the sample page (with the use of the SRA, SRX or SRS number).
		'''
		for EGEOD in self.accessions:
			if EGEOD not in self.knownAccessions:
				print(EGEOD)
				# HTML content is obtained for each study.
				htmlContent = self.connectHtml.obtainHtmlContent("http://www.ebi.ac.uk/arrayexpress/experiments/"+EGEOD)
				# Obtaining the sample information from the SDRF file.
				sampleInfo = self.sdrfObtainer.downloadSdrf(EGEOD)
				if sampleInfo:
					sra =""
					# The author and year are obtained by the function self.obtainSampleInfo().
					authorYear = self.obtainSampleInfo(htmlContent)
					for item in sampleInfo:
						if sra != item[-1]:
							sra = item[-1]
							# Whitespaces are removed from the column containing the SRA, SRX or SRS number.
							item[-1] = item[-1].split()
							# The html content of the sample NCBI page is obtained
							htmlContent = self.connectHtml.obtainHtmlContent("http://www.ncbi.nlm.nih.gov/sra?term="+item[-1][0])
							# Necessary information is obtained from this page.
							genericInfo = self.obtainGlobalInfo(htmlContent)					
							# Information about the study and sample is written into the Xlsx file.
							self.xlsxWriter.writeToXlsx([EGEOD, authorYear[1], authorYear[0], genericInfo[0], genericInfo[4], " ", genericInfo[1], genericInfo[2], genericInfo[3], item[3], item[0], item[2], item[1], genericInfo[5]])

	def obtainSampleInfo(self,htmlContent):
		'''
			The function obtainSampleInfo is used to obtain the information abou the publication year and the first author from the study.
			The HTML content is used to obtain this information.
			Since arrayexpress defines the author in several different ways.
		'''
		author = []
		year = re.search("<div>Status</div></td><td class=\"value\"><div>[A-z ]+[0-9 ]+[A-z ]+([0-9]+)", htmlContent)
		author = re.search('Citation.{41}="http:\/\/.+">.{1,900}<\/a>[.] ([A-Za-z ]+),', htmlContent)
		if not author:
			re.search('Citation.?</div></td><td class="value"><div>[A-z ]+. ([A-z ]+),', htmlContent)
			if not author:
				author = re.search('Contact.{1,130}>([A-Za-z .-\u00C0-\u00ff]+)?[,|&]', repr(htmlContent))
				if not author:
					# The last re pattern is used when the name is not found with the pattern above. This is possible in some rare cases.
					author = re.search('Contact.{1,130}>([A-Za-z \.\-]+)[,|&]?(</div></td></tr>)?', htmlContent)
					if not author:					
						author = re.search('Contact.{1,130}\>[A-Za-z \.\-\u00C0-\u00ff]+?[,|&]lt\;[A-z\.\@\&\;]+\<\/a\>, ([A-z ]+),', repr(htmlContent))
		# Data is saved and changed into a capitalized string.
		author = str(author.group(1)).strip().split(" ")
		return year.group(1), str(author[0]).capitalize() +" et al."


	def obtainGlobalInfo(self, htmlContent):
		'''
			Generic information such as the title, organism, intrument, paired/single end and abstract information is saved.
		'''
		self.genericInfo = re.findall('Study: <span>([A-z0-9 ,\-\(\)\'\:\;\.ü]+).+Organism: <span><a href=.+>([A-z0-9 \.\-]+)</a>.+Instrument: <span>([A-z0-9 ]+)</span>.+Layout: <span>([A-z]+)</span>', htmlContent)
		if not self.genericInfo:
			# The following regex pattern is used when the organism is unknown.
			self.genericInfo = re.findall('Study: <span>([A-z0-9 ,\-\(\)\'\:\;\.ü]+).+Instrument: <span>([A-z0-9 ]+)</span>.+Layout: <span>([A-z]+)</span>', htmlContent)
			self.genericInfo.insert(1, "Unknown")
		# The abstract of the study is obtained.
		abstract = re.findall('Abstract</span><span class="less">hide Abstract</span></a><div class="expand-body">([A-z0-9 \.\,\:\(\)\-\+\/\'ü\’=;\?]+)</div></div></span></div>', htmlContent)
		# The SRA number is retreived.
		sra = re.search('run=(SRR[0-9]+)', htmlContent)
		if abstract:
			self.abstract = abstract[0]
			self.sample = " "
		else:
			self.abstract = " "
			sample = re.findall("Sample: <span>([A-z0-9 ]+.+)<div class=\"expand-body\"><a href=", htmlContent)
			self.sample = " "
		# Information is saved into the self.genericInfo list.
		self.genericInfo[0] = list(self.genericInfo[0])
		self.genericInfo[0].insert(len(self.genericInfo[0]), self.abstract)
		self.genericInfo[0].insert(len(self.genericInfo[0]), sra.group(1))
		self.genericInfo = self.genericInfo[0]
		return self.genericInfo

if __name__ == "__main__":
	co = CreateOverview()

