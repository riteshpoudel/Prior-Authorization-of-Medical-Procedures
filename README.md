# Prior-Authorization-of-Medical-Procedures

This study used a cross-sectional statistical analysis of health plan administrative data for two months data of PAs submitted by physicians. 
These data consisted of 9393 PA records for the months July and August and January 
2013, the two months data elements had been recorded for the health plan. 
Although this was two months of recording, the PA program had been in existence for many years. Thus, these data do not represent new or extraordinary circumstances for the program, although this is one of the limitations of the study. 
Data analysis consisted of calculating
•Univariate frequencies,  Bivariate Analysis
And the calculation of a binary logistic regression model of the likelihood that a request for PA would be denied, given drug types and patient demographics. 
The binary multiple logistic regressions and other machine learning classification algorithms are chosen for the main data analysis and selected the best for prediction, because it helps to identify which factors independently predict , if the prescribed medication requires PA. 
Administrative data were used in this project because they contain information about every request for PA. These data are limited since the diagnosis and reasons for the request are not recorded. 
The steps involved in the process of building the system for PA prediction: 
Data Exploration 
Pre-Processing 
Apply classification models 
Validate the results 
And Select the best model and Make Prediction 
The Algorithms used to train and test the data and compare the results between the models and select the best performing algorithms: 
• GLM 
• Decision Tree 
• SVM 
• Random Forest 

•H2o , AutoEncoder 
• H20 Deeplearning 
• Stacking (GLM, rpart, SVM) 

Bottlenecks: 
The data is fully categorical, and having too many levels (as shown in the fig1) that a classification algorithm cannot handle so easily, as the attributes will grow 10000+ by dummying the levels, we will have to face the Curse of Dimensionality. 
Now the challenges are to reduce the levels and make the data set to run algorithms. 
Solution for the bottlenecks: 
Methodologies adopted to reduce the levels in each individual independent variable are: 
1) Take the target frequency for each level and cluster it using K-Means clustering technique and assign the cluster id for the levels accordingly. 

2) Take the true count and divide it by the total frequency of a particular level and assign the true response rate build models on it. 

Method 1: As the levels of each individual variables is very high, reducing the levels is the 
Taking the frequency of the Dependent variable against each Independent variable, and clustering it using k-means clustering algorithm, so that the levels will be dramatically reduced from high to low (for e.g. from 1609 to 4-7) 
Create frequency table for each individual attribute against the Target, standardize the frequency, and apply k-means to find the clusters, check the elbow curve and identify the optimized K value for the K-Means algorithm, once clustered the individual attribute, then applying the cluster id to the levels accordingly will reduces the huge levels in the data. 
Observations: The observation at this stage is that, table() function is used to get the target frequency against the level of particular attribute. 
e.g., 
table(PriorAuth_Data$Drug , PriorAuth_Data$Target) 
OutPut: FALSE 	TRUE 
Drug001 	6 	10 
Drug002 	2 	1 
Drug003 	8 	28 
Drug004 	2 	1 
Drug005 	3 	10 
Drug006 	2 	4 
Drug007 	32 	100 

After converting it to data frame it looks below 
Data.frame(table(PriorAuth_Data$Drug,PriorAuth_Data$Target)) 

Out Put: 
Var1 Var2 Freq 
Drug001 FALSE 6 
Drug002 FALSE 2 
Drug003 FALSE 8 
Drug004 FALSE 2 
Drug005 FALSE 3 
Drug006 FALSE 2 
Drug007 FALSE 32 

As I need to normalize the data before clustering it, we need to exclude the level data which is categorical so converting it to a data frame is required, and the data frame conversion results as above. To overcome the above issue, I have used sqldf() , so that the SQL queries can be executed and the data frame can be created with the data in required structure. 
After getting the frequency table for each independent variable, normalizing the data and applying k-means to get the cluster id, which will be replaced with the level data in the data set to reduce the levels drastically, from 2000 to 5-7. 
Method2: 
Take the true count and divide it by the total frequency of a particular level and assign the true response rate build models on it. 
As just clustering the frequency was not very intuitive, and the correlation between the target and the independent variables of the data is not reflecting the actual relationship between the attributes, taking the true frequency and calculating the probability against the total frequency of that particular level speaks the relationship well between the variables as it is intended to as shown below. 


Data: 
Data is gathered from the 
Attribute Analysis and insights to introduce new variables, which help in prediction prior authorization for the prescription. 
By doing the grouped data analysis, and also as per the definition of the variables such as: 
GPI : Generic Product Identifier 
NDC : National Drug Code 
These two (GPI , NDC ) are the identifiers and are for similar purpose to identify the nature of the drug, its group, its name, its vendor, its product code, the drug class and sub class..!

And if we see the data, we have the individual components of the combinations mentioned for GPI and NDC, which will be having correlation between the grouped data and as well as the GPI and NDC, by which the dimensions that needs to be passed to build the predictive model is even reduced from 14 to 3 or 4, which makes the model’s job easier and the computation will also be significantly low if the data grows periodically. 

