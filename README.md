# GliomaDataBase
This Repo is a Tool set to analyse the IvyGap Giloma Database 

Review of Repo Components

User Interface tools and Explaination 

The User Interface tools (found in the AppScripts folder in the GalaxyToolBranch of this Repo) were created to make statistical analysis of IvyGap data more efficient and User Friendly
Each Script runs statistical anyalysis of the Gene Expression of a provided Gene and the name sake varible of the script, Tumour Region, Molecular Subtype, or Survival Time.

If you are a first time user follow there steps to use this interface locally 
1. Make a Copy of this Repo into the code editor of your choice which supports python
2. Ensure the the IvyGapInformation File provided in this script is downloaded in an accessible location; this will be refered to as file 1
3. Go to the Ivy Gap Glioma Database website and use RNA Seq to find the gene you would like to analyse; a zip folder will be downloaded, save the file titled Expression; this will be refered to as file 2
4. Open the Copy of the Repo in your code editor and click on the folder titled AppScripts
5. Find the Varible you would like to analysis and open the script
6. Check that all nessary libraries are installed to run the script; the first lines outlin all the imports
7. Ensure Streamlit in downloaded on you computer
8. After all libraries are downloaded (which is only needed when using the tool for the first time) open an instance in the terminal and run the following lines of code
    cd AppScripts (This may changes based on how the script is saved just ensure you are working in the directory where the python file is saved)
    streamlit run TumorRegionStatsSummaryUI.py (or the script of your choosing)
9. This should open a window for the Streamlit app in your web browser
10. Then use the drag and drop feature to put file 1 in IVYGap file location and File in the Gene Expression drop menu

Note: using this tool the first time will take the longest as everything must be downloaded locally, after the first use you will only need to open the script in a streamlit instance and then it can be used repeated in the web browser

User Interface Users may stop here 


Mannual Python Script Information 

Base Python Scripts- These are the orginal Python Scripts that the User interface is built based off of; these files require mannual editing

RegionExpressionBoxPlots

Files-> comparing Tumor Regions and Gene Expression are included and are example of comparison of IvyGap Data which can be done by this script 

WorkflowScripts

Files 
IvyTumourInformation-This File contains the Patient information of each Tumor Sample provided by the IvyGap Database 
new_file.csv- This File is contains that Patient Information combine with selected Gene expression; It is edited from the ExpressionInformation script
ANOVAResults.csv- This File is assigned the ANOVA results from StatsAnalysis Script
DunnResults.csv- This File is assigned the Dunnes results from StatsAnalysis Script
SumStat.csv-This File is assigned the Summary Statistics results from StatsAnalysis Script

Scripts
ExpressionInformation- This File is used to Format the expression information of the gene or genes desired to be analysised and combine it with the Patient and tumor Sample Information; the user is require to download and edit path infomation to suit local computer 
StatAnalysis- This File Runs summary Statistics, ANOVA, Shapiro-Wilk, Kruskal-Wallis, and Dunns Post Hoc Tests on the values from file composed from new_file.csv
BoxplotTrial- This File Runs Using ExpressionInformation and StatAnalysis scripts to create summary Box Plots Representative of expression in various tumour regions; Labels may need to be edited 


To use this tool, you must download this file and edit the path according to the computer being used 


