# graduation-project
Scripts provided in this repository were used during the development of GOAD

#Preprocessing steps
Files found within this folder were used to preprocessing the data.<br />
This folder contains the following content:<br />
  - createOverview        (Folder)<br />
    Where the scripts can be found that can be used to create an overview for GOAD.<br />
    All of the functions were called by the script "collectInfoPerStudy.py".<br />
  - createSimpleOverview  (Python Script)<br /> 
    Used to create a overview without all sample information, this overview can be used to validate the studies.<br />
  - sampleSheetCreate     (Python Script)<br />
    This script can be used to create an sampleSheet that can be used within the Kallisto pipeline of MOLGENIS compute.<br />
  - mergeCounts           (Bash Script)<br />
    Counts of samples are saved within one file per sample, this script enables the merge of these files into one count file.<br />
    <br />
#htmlContent
The folder consists of html, css and javascript files that are used to create the GOAD web application.<br />
  - QE_Bargraph           (CSS)<br />
    Contains the css for the bargraph shown with QE data.<br />
  - QE_Dashboard          (CSS)<br />
    Contains the css for the dashboard shown with QE data.<br />
  - goad                  (CSS)<br />
    Contains the css for GOAD.<br />
<br />  
  - goad-contactPage      (HTML)<br />
    Contains the contact information part of the web application.<br />
  - goad-informationPage  (HTML)<br />
    Contains information about the transformation of the raw data into data for the web application.<br />
  - goad-publicationPage  (HTML)<br />
    Contains the publication part of the web application.<br />
  - goad-tutorialPage     (HTML)<br />
    Contains the tutorial information for the GOAD web application.<br />
  - view-goadmanager      (HTML)<br />
    Contains all of the HTML content.<br />
  <br />
  - goadBarGraph          (JavaScript)<br />
    JavaScript content to use the QE Bargraph.<br />
  - goadDashboard         (JavaScript)<br />
    JavaScript content to use the QE Dashboard.<br />
  - goadFunctions         (JavaScript)<br />
    Contains all of the functions that where used on the GOAD web application.<br />
  - goadPublicationPart   (JavaScript)<br />
    Contains the code needed for the publication part of the web application.<br />
  - goadmanager           (JavaScript)<br />
    Contains the rest of the code for the GOAD web application.<br />

#rScripts
This folder consists of R scripts that were used with the GOAD web application.<br />
  - CalculateTPMbyGTF     (R Script)<br />
    Makes sure that the TPM values are calculated with the use of the GTF files (and the raw counts within the MOLGENIS database).<br />
  - datatableDE           (R Script)<br />
    This R script is used by the R API of MOLGENIS, making sure that raw counts are transformed into logFC and FDR values.<br />
    These values can be used to write into a table.<br />
  - scatterplotDE         (R Script)<br />
    This R Script is used by the R API of MOLGENIS, creating an interactive D3js scatterplot in the end.<br />
  
