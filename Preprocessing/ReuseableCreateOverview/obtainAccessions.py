#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# Title			: 	obtainAccessions.py
# Created by	: 	Marissa Dubbelaar
# Created on	: 	10-02-2016
# Purpose		:	1) Obtain the unique Accession numbers

import re

class AccessionObtainer():
	def __init__(self):
		'''Makes sure that the class HtmlConnector can be used to retreive information from a given HTML page.'''
		self.uniqueAccession = set()

	def getUniqueAccessionName(self, htmlContent):
		'''
		The function getUniqueAccessionName is used to obtain the E-GEOD numbers that are found on arrayexpress with the RNA-seq specification.
		The HTML of the given page is saved into a list with multiple strings, with the use of the join function the html is joined to make sure that
		interesting information can be obtained more easely with the use of a regex pattern.
		The found E-GEOD numbers are saved into the list: self.uniqueAccession and can be used further in other defs.  
		'''
		# Opens the url with all of the files that are connected to the right cell type.
		EGEODnumbers = re.findall('a href="/arrayexpress/experiments/(E-\w{3,4}-[\d]+)/\?[\w\=\&\;\%\+\)]+"', htmlContent)
		for EGEODNr in EGEODnumbers:
			if EGEODNr:
				# The EGEODNr is added into the list when it is not known yet in the self.knownAccessions
				self.uniqueAccession.add(EGEODNr)
			else:
				print("No new studies found")
		return self.uniqueAccession
