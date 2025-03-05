# GliomaDataBase
This Repo is a Tool set to analyse the IvyGap Giloma Database 

Review of Repo Components

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
