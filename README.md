# graduation-project
Scripts provided in this repository were used during the development of GOAD

#Preprocessing steps
Files found within this folder were used to preprocessing the data.
This folder contains the following content:
  - createOverview        (Folder)
    Where the scripts can be found that can be used to create an overview for GOAD.
    All of the functions were called by the script "collectInfoPerStudy.py".
  - createSimpleOverview  (Python Script) 
    Used to create a overview without all sample information, this overview can be used to validate the studies.
  - sampleSheetCreate     (Python Script)
    This script can be used to create an sampleSheet that can be used within the Kallisto pipeline of MOLGENIS compute.
  - mergeCounts           (Bash Script)
    Counts of samples are saved within one file per sample, this script enables the merge of these files into one count file.
    
#htmlContent
The folder consists of html, css and javascript files that are used to create the GOAD web application.
  - QE_Bargraph           (CSS)
    Contains the css for the bargraph shown with QE data.
  - QE_Dashboard          (CSS)
    Contains the css for the dashboard shown with QE data.
  - goad                  (CSS)
    Contains the css for GOAD.
  
  - goad-contactPage      (HTML)
    Contains the contact information part of the web application.
  - goad-informationPage  (HTML)
    Contains information about the transformation of the raw data into data for the web application.
  - goad-publicationPage  (HTML)
    Contains the publication part of the web application.
  - goad-tutorialPage     (HTML)
    Contains the tutorial information for the GOAD web application.
  - view-goadmanager      (HTML)
    Contains all of the HTML content.
  
  - goadBarGraph          (JavaScript)
    JavaScript content to use the QE Bargraph.
  - goadDashboard         (JavaScript)
    JavaScript content to use the QE Dashboard.
  - goadFunctions         (JavaScript)
    Contains all of the functions that where used on the GOAD web application.
  - goadPublicationPart   (JavaScript)
    Contains the code needed for the publication part of the web application.
  - goadmanager           (JavaScript)
    Contains the rest of the code for the GOAD web application

#rScripts
This folder consists of R scripts that were used with the GOAD web application
  - CalculateTPMbyGTF     (R Script)
    Makes sure that the TPM values are calculated with the use of the GTF files (and the raw counts within the MOLGENIS database).
  - datatableDE           (R Script)
    This R script is used by the R API of MOLGENIS, making sure that raw counts are transformed into logFC and FDR values.
    These values can be used to write into a table.
  - scatterplotDE         (R Script)
    This R Script is used by the R API of MOLGENIS, creating an interactive D3js scatterplot in the end.
  
