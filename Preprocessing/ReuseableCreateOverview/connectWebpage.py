# Title			: 	connectWebpage.py
# Created by	: 	Marissa Dubbelaar
# Created on	: 	12-02-2016
# Purpose		:	1) Make connection to the given website.
#					2) Return the HTML content of the site.

import urllib2
import time
import socket

class HtmlConnector():
	def __init__(self):
		self.html_contents = []

	def obtainHtmlContent(self, linkToSite):
		'''
		The function obtainHtmlContent makes sure that the connection is make a connection to the given
		HTML link. To make sure that this process is not seen as an attack on the website, a second of sleep
		is used. Headers are added to secure the accessebility to the website.
		The html content is read in the end and is returned.
		'''
		count = 0
		try:
			opener = urllib2.build_opener()
			# Change the user agent.
			opener.addheaders = [('User-Agent', 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_9_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/35.0.1916.47 Safari/537.36')]
			response = opener.open(linkToSite)
			# Obtain the HTML content
			html_contents = response.read()
			# Closes the connection with the HTML page
			response.close()
			return html_contents
		except socket.error:
			# Socket.error is raised when the pages are opened to soon after each other.
			# A second of sleep is included and the function is called again.
			time.sleep(1)
			obtainHtmlContent()
		except:
			# This exception is raised when there is no socket.error
			print("Something went wrong while connecting to the html page: " + linkToSite+"\nPlease try again!")
			exit()

# Test
# Prints the htmlContent when this code is called as __main__
if __name__ == "__main__":
	htmlCon = HtmlConnector()
	htmlContent = htmlCon.obtainHtmlContent("https://www.ebi.ac.uk/arrayexpress/search.html?query=microglia&exptype%5B%5D=%22rna+assay%22&exptype%5B%5D=%22sequencing+assay%22&page=1&pagesize=500")
	print(htmlContent)	
