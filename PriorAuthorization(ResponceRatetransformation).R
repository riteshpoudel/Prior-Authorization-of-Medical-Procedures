##using all data, without removing row level duplicates

##Remove all environment variables, to relinguis the RAM space and also to avoid overlaps...

rm(list=ls(all=TRUE))

#Uncomment the installer part if the packages are not available
#pkgs <- c("readr",
#          "spam","LICORS",
#          "cluster","vegan","sqldf","C50","ROCR","rpart","e1071","ggplot2","ggthemes","scales","h2o")

#install.packages(pkgs)


library(readr)
# library(LICORS)
# library(spam)
library(vegan)
library(cluster)
#Try usig this later all in the end
library(sqldf)
library(RSQLite)
library(tcltk)
library(ROCR)
library(rpart) # classification algorithm
#library(C50) # classification algorithm
library(e1071)# classification algorithm
library(ggplot2) # visualization
library(ggthemes) # visualization
library(scales) # visualization
library(randomForest) # classification algorithm
library(h2o) # Neural Network
setwd("E:/Workspace/Project")

getwd()

#Loading data into a DataFrame using read_csv, as to apply format to the date variable while loading itself.
#as I have abserved that the date is having diffrent format in the data set  provided.
PriorAuth_Data <- read_csv("PriorAuth_Data.csv",
                           col_types = cols(TransDate = col_date(format = "%m/%d/%Y")))

sum(is.na(PriorAuth_Data)) # no NA 's found

PriorAuth_Data <- PriorAuth_Data[setdiff(names(PriorAuth_Data),'UserID')]

head(PriorAuth_Data)
str(PriorAuth_Data)
# By Definition GPI Generic Product Identifier :
# Combination of Drug group + Drug class + Drug sub-class + Drug name + Drug name extension
PriorAuth_Data$comb1 <- paste(PriorAuth_Data$Drug
                              ,PriorAuth_Data$DrugClass
                              ,PriorAuth_Data$DrugSubClass
                              ,PriorAuth_Data$Drug_Chemical_Name
                              ,PriorAuth_Data$DrugGroup,sep = '')

# Seeing into the data, another Combination of
#Drug group + Drug class + Drug sub-class + Drug name + Drug name extension + GPI gives the same NDC Code.
#So forming another combination attribute...!
PriorAuth_Data$comb2 <- paste(PriorAuth_Data$Drug
                              ,PriorAuth_Data$DrugClass
                              ,PriorAuth_Data$DrugSubClass
                              ,PriorAuth_Data$Drug_Chemical_Name
                              ,PriorAuth_Data$DrugGroup
                              ,PriorAuth_Data$GPI,sep = '')

PriorAuth_Data$comb3 <- paste(PriorAuth_Data$RxGroupId
                              ,PriorAuth_Data$Bin
                              ,PriorAuth_Data$PCN
                              ,sep = '')

##Create frequency table for each individual attribute against the Target, and standardising the freque
##to find the clusters, to reduce the huge levels in the data.

########################################################################################
# Whys SQLDF, why not table() for getting frequency ?
# Table is the easy way to get frequency, but then while converting it to dartaframe,
# the data is transformed undesired manner.
#eg: samle table(PriorAuth_Data$Drug,PriorAuth_Data$Target) out put
#           FALSE TRUE
# Drug001     6   10
# Drug002     2    1
# Drug003     8   28
# Drug004     2    1
# Drug005     3   10
# Drug006     2    4
# Drug007    32  100
# after converting it to dataframe it looks below
# Data.frame(table(PriorAuth_Data$Drug,PriorAuth_Data$Target))
#    Var1  Var2 Freq
# Drug001 FALSE    6
# Drug002 FALSE    2
# Drug003 FALSE    8
# Drug004 FALSE    2
# Drug005 FALSE    3
# Drug006 FALSE    2
# Drug007 FALSE   32
########################################################################################

attributes <- setdiff(names(PriorAuth_Data),c('Target','TransDate'))

rm(list= ls()[!(ls() %in% c('PriorAuth_Data','attributes'))])

# Taking the count of True and total count for each attribute and its levels and further taking 
# "True Responce Rate" to handle the high leveled categorical data.

select <- c('Select ')
case <- ', sum(case when Target = 1 then 1 else 0 end) as True, sum(case when Target = 1 then 1 else 0 end)+(sum(case when Target = 0 then 1 else 0 end)) as Total from PriorAuth_Data group by '

for (i in 1:length(attributes)){
  
  ##Creating Frequency tables for all Individual attributes against the Class variable dynamically
  assign(attributes[i],sqldf(paste(select,attributes[i],case,attributes[i],sep = '')))
  
}

## Calculating True responce rate for each attribute.
Bin                 <-cbind(Bin[1],Bin[2]/Bin[3])
Drug                <-cbind(Drug[1],Drug[2]/Drug[3])
DrugSubClass        <-cbind(DrugSubClass[1],DrugSubClass[2]/DrugSubClass[3])
DrugClass           <-cbind(DrugClass[1],DrugClass[2]/DrugClass[3])
Drug_Chemical_Name  <-cbind(Drug_Chemical_Name[1],Drug_Chemical_Name[2]/Drug_Chemical_Name[3])
GPI                 <-cbind(GPI[1],GPI[2]/GPI[3])
NDC                 <-cbind(NDC[1],NDC[2]/NDC[3])
DrugGroup           <-cbind(DrugGroup[1],DrugGroup[2]/DrugGroup[3])
DoctorID            <-cbind(DoctorID[1],DoctorID[2]/DoctorID[3])
RxGroupId           <-cbind(RxGroupId[1],RxGroupId[2]/RxGroupId[3])
PCN                 <-cbind(PCN[1],PCN[2]/PCN[3])
State               <-cbind(State[1],State[2]/State[3])
comb1               <-cbind(comb1[1],comb1[2]/comb1[3])
comb2               <-cbind(comb2[1],comb2[2]/comb2[3])
comb3               <-cbind(comb3[1],comb3[2]/comb3[3])

rm(select,case,i)

PriorAuth_num_Data <- sqldf('select d.True DrugId,sc.True DrugSubClassId , dc.True DrugClassId ,c.True Drug_ChemicalId, g.True GpiId , n.True NdcId , dg.True DrugGroupId,
                            Dr.True DoctorID , gi.True RxGroupId , b.True BinId, p.True PcnId, s.True StateId, c1.True com1Id, c2.True com2Id ,c3.True com3Id , PA.Target,PA.TransDate
                            from PriorAuth_Data PA
                            inner join Drug d on PA.Drug = d.Drug
                            inner join DrugSubClass sc on PA.DrugSubClass = sc.DrugSubClass
                            inner join DrugClass dc on PA.DrugClass = dc.DrugClass
                            inner join Drug_Chemical_Name c on PA.Drug_Chemical_Name = c.Drug_Chemical_Name
                            inner join GPI g on PA.GPI = g.GPI
                            inner join NDC n on PA.NDC = n.NDC
                            inner join DrugGroup dg on PA.DrugGroup = dg.DrugGroup
                            inner join DoctorID Dr on PA.DoctorID = Dr.DoctorID
                            inner join RxGroupId gi on PA.RxGroupId = gi.RxGroupId
                            inner join Bin b on PA.Bin = b.Bin
                            inner join PCN p on PA.PCN = p.PCN
                            inner join State s on PA.State = s.State
                            inner join comb1 c1 on PA.comb1 = c1.comb1
                            inner join comb2 c2 on PA.comb2 = c2.comb2
                            inner join comb3 c3 on PA.comb3 = c3.comb3
                            ')

#Removing unwanted/unused variables and relinquishing memory.
rm(list= ls()[!(ls() %in% c('PriorAuth_num_Data'))])

# as all other attributes are now integer for except Trandate, converting Target to Integer
PriorAuth_num_Data$Target <- ifelse(PriorAuth_num_Data$Target == TRUE, 1, 0)
str(PriorAuth_num_Data)
#for some reasons the below two attributes are chr, so converting it to num
PriorAuth_num_Data$DoctorID<- as.numeric(PriorAuth_num_Data$DoctorID)
PriorAuth_num_Data$RxGroupId <- as.numeric(PriorAuth_num_Data$RxGroupId)

str(PriorAuth_num_Data)

cor(PriorAuth_num_Data[setdiff(names(PriorAuth_num_Data),c('TransDate'))])

#########################################################################################
#Data is Ready for Model building
#########################################################################################
#Split Train and test using Time-based sampling approach as suggested
# Total 18 days data
# data grouped by date info
data.frame(table(PriorAuth_num_Data$TransDate))

traindata <- subset(PriorAuth_num_Data, TransDate<=as.Date("2013-08-04"))
nrow(traindata) # 7255 rows for training i.e 77%
testdata <- subset(PriorAuth_num_Data, TransDate>as.Date("2013-08-04"))
nrow(testdata)#2138 rows for testing teh model i.e 23%
str(traindata)

rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

#########################################################################################
### Building Modelz
#########################################################################################
### Logistic Regression
#########################################################################################

# LogReg <- glm(Target ~ ., data=traindata, family=binomial)
# summary(LogReg)
#
step(glm(Target ~ .,data=traindata),direction = 'backward')
# step gave me NdcId + DoctorID + com3Id +Transdate, But removing Transdate gave a slight improvement.
LogReg <- glm(formula = Target ~ NdcId + DoctorID + com3Id ,
              data = traindata)

# train results
prob<-predict(LogReg, type="response")
pred_class <- ifelse(prob> 0.5, 1, 0)
table(traindata$Target,pred_class)

# Error Metric

conf.mat = table(traindata$Target,pred_class)
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

# Test results
fitted.results <- predict(LogReg,testdata,type='response')
fitted.class <- ifelse(fitted.results > 0.5,1,0)
table(testdata$Target,fitted.class)

# Error Metric
conf.mat = table(testdata$Target,fitted.class)
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

#Ploting the ROC curve and calculate the AUC
#(area under the curve) which are typical performance measurements
#for a binary classifier.
#The ROC (Receiver Operating Characteristic curve) is a curve generated by plotting the true positive rate (TPR = sensitivity) against
# the false positive rate (FPR= specificity) at various threshold settings while the AUC is
# the area under the ROC curve. As a rule of thumb, a model with good
#predictive ability should have an AUC closer to 1 (1 is ideal) than to 0.5.


p <- predict(LogReg,testdata, type="response")
pr <- prediction(p, testdata$Target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf,colorize = TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7))

abline(a=0, b= 1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc 

#########################################################################################
### Decision Tree
#########################################################################################

rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

PA_rpart <- rpart(Target~ ., data=traindata, method="class")
#By Variable Importance got below as top  
PA_rpart <- rpart(Target~DoctorID + NdcId + com3Id, data=traindata, method="class")
plot(PA_rpart,main="Classification Tree for loan Class",margin=0.15,uniform=TRUE)
text(PA_rpart,use.n=T)
summary(PA_rpart)

conf.mat = table(traindata$Target, predict(PA_rpart, newdata=traindata, type="class"))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

conf.mat = table(testdata$Target, predict(PA_rpart, newdata=testdata, type="class"))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

# Compared the results of glm and rpart,  it is almost the same
#write.csv(cbind(testdata$Target, testresult),'desition_testresult.csv')
#write.csv(cbind(testdata$Target,fitted.class),'glm_testresult.csv')

#########################################################################################
### SVM
#########################################################################################
rm(list= ls()[!(ls() %in% c('traindata','testdata'))])
# Build best SVM model
PA_SVM <- svm(Target~ NdcId + DoctorID + com3Id, data=traindata,  kernel = "linear")

# Look at the model summary
summary(PA_SVM)

plot(PA_SVM$index)

# Predict on train data
pred_Train  =  predict(PA_SVM, traindata)

plot(pred_Train) # Plot shows more than 0.5 

conf.mat = table(traindata$Target, ifelse(pred_Train> 0.5, 1, 0))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

# Predict on test data
pred_Test  =  predict(PA_SVM, testdata[setdiff(names(testdata),c('Target'))])
conf.mat = table(testdata$Target, ifelse(pred_Test> 0.5, 1, 0))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

#write.csv(cbind(testdata$Target, ifelse(pred_Test> 0.5, 1, 0)),'SVM_results.csv')

######################################################################################################
# Random Forest
######################################################################################################

rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

pa_rf <- randomForest(factor(Target) ~ ., data = traindata)

str(traindata)

pa_rf$importance  
round(importance(pa_rf), 2)   

# Extract and store important variables obtained from the random forest model
rf_Imp_Attr = data.frame(pa_rf$importance)
rf_Imp_Attr = data.frame(row.names(rf_Imp_Attr),rf_Imp_Attr[,1])
colnames(rf_Imp_Attr) = c('Attributes', 'Importance')
rf_Imp_Attr = rf_Imp_Attr[order(rf_Imp_Attr$Importance, decreasing = TRUE),]

# plot (directly prints the important attributes) 
varImpPlot(pa_rf)

# Importance is teh same as my other models, just that the Trans date is having higher importance, still 
# I am not including it in my model, because it is giving very poor performance in test as espected...
pa_rf <- randomForest(factor(Target) ~ NdcId + DoctorID + com3Id, data = traindata)

# Predict on Train data
pred_Train = predict(pa_rf, traindata[,setdiff(names(traindata),"Target")],
                     type="response", norm.votes=TRUE)

conf.mat = table(traindata$Target, pred_Train)

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

# Predict on test data
pred_Test  =  predict(pa_rf, testdata[,setdiff(names(testdata),"Target")],
                      type="response", norm.votes=TRUE)

conf.mat = table(testdata$Target, pred_Test)

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

######################################################################################################
### H2o Autoencoder
######################################################################################################
rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

h2o.init(ip='localhost', port = 54321, max_mem_size = '2g')

train.hex <- as.h2o(x = traindata, destination_frame = "train.hex")
test.hex <- as.h2o(x = testdata, destination_frame = "test.hex")

y = "Target"
x = setdiff(colnames(train.hex), y)

aec <- h2o.deeplearning(x = x, autoencoder = T, 
                        training_frame=train.hex,
                        activation = "Tanh",
                        hidden = c(20),
                        epochs = 100)

# Extract features from train data
features_train <- as.data.frame(h2o.deepfeatures(data = train.hex[,x], object = aec))
# Extract features from test data
features_test <- as.data.frame(h2o.deepfeatures(data = test.hex[,x], object = aec))

train_new <-data.frame(traindata,features_train)
test_new <-data.frame(testdata,features_test)

rf_DL <- randomForest(Target ~ ., data=train_new, keep.forest=TRUE, ntree=50)

# importance of attributes
round(importance(rf_DL), 2)
importanceValues = data.frame(attribute=rownames(round(importance(rf_DL), 2)),MeanDecreaseGini = round(importance(rf_DL), 2))
row.names(importanceValues)=NULL
importanceValues = importanceValues[order(-importanceValues$IncNodePurity),]
# Top 10 Important attributes
Top10ImpAttrs = as.character(importanceValues$attribute[1:16])

Top10ImpAttrs

train_Imp = subset(train_new,select = c(Top10ImpAttrs,"Target"))
test_Imp = subset(test_new,select = c(Top10ImpAttrs,"Target"))

rm(train_new,test_new)
# Build the classification model
model_Imp = glm(Target ~ ., data = train_Imp,family = 'binomial')


# train results
prob<-predict(model_Imp, type="response")
pred_class <- ifelse(prob> 0.5, 1, 0)
table(train_Imp$Target,pred_class)

# Error Metric

conf.mat = table(train_Imp$Target,pred_class)
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

# Test results
fitted.results <- predict(model_Imp,test_Imp,type='response')
fitted.class <- ifelse(fitted.results > 0.5,1,0)
table(test_Imp$Target,fitted.class)

# Error Metric
conf.mat = table(test_Imp$Target,fitted.class)
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

######################################################################################
## Deep learning with two hidden nodes
######################################################################################

rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

h2o.init(ip='localhost', port = 54321, max_mem_size = '2g')

y = "Target"
x = c('NdcId' , 'DoctorID' , 'com3Id')

train.hex <- as.h2o(x = traindata, destination_frame = "train.hex")
test.hex <- as.h2o(x = testdata, destination_frame = "test.hex")

dlModel <- h2o.deeplearning(x = x, y = y,
                            training_frame=train.hex,
                            activation = "TanhWithDropout",
                            hidden = c(150,150),
                            input_dropout_ratio = 0.001,
                            l1 = 1e-10,
                            epochs = 100)

# View specified parameters of the deep learning model
dlModel@parameters

# Examine the performance of the trained model
dlModel # display all performance metrics

# Metrics
h2o.performance(dlModel)

# Get MSE only
h2o.mse(dlModel)

pred = h2o.predict(dlModel, train.hex[x])

preddf <- as.data.frame(pred)

conf.mat <- table(traindata$Target,ifelse(preddf$predict > 0.5,1,0))
conf.mat
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))


#----Prediction on test data

pred1 = h2o.predict(dlModel, test.hex[x])

pred_test_df <- as.data.frame(pred1)

conf.mat <- table(testdata$Target,ifelse(pred_test_df$predict > 0.5,1,0))
conf.mat
cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

write.csv(x = ifelse(pred_test_df$predict > 0.5,1,0),file = "DL.csv")


#######################################################################
#Stagging 
#######################################################################

#Glm
LogReg <- glm(formula = Target ~ NdcId + DoctorID + com3Id ,
              data = traindata)
#SVM
PA_SVM <- svm(Target~ NdcId + DoctorID + com3Id, data=traindata,  kernel = "linear")

#Decision Tree
PA_rpart <- rpart(Target~ DoctorID + NdcId + com3Id, data=traindata, method="class")

# Prediction on Train

# Train results
glm_train<-predict(LogReg, type="response")
glm_Train <- ifelse(glm_train> 0.5, 1, 0)
table(traindata$Target,glm_Train)

rpart_train <- predict(PA_rpart, newdata=traindata, type="class")
as.vector(rpart_train)


SVM_Train  =  predict(PA_SVM, traindata)
SVM_Train <- ifelse(SVM_Train> 0.5, 1, 0)

#combine 

train_Pred_All_Models = data.frame(GLM = glm_Train, 
                                   rpart = rpart_train,
                                   SVM = SVM_Train)
train_Pred_All_Models = data.frame(sapply(train_Pred_All_Models, as.factor))

str(train_Pred_All_Models)
summary(train_Pred_All_Models)
rm(rpart_Train, glm_Train, SVM_Train)

table(train_Pred_All_Models$GLM) #GLM
table(train_Pred_All_Models$rpart)  #rpart
table(train_Pred_All_Models$SVM)  #SVm
table(traindata$Target) #Original Dataset DV

train_Pred_All_Models = cbind(train_Pred_All_Models, Target = traindata$Target)

# Ensemble Model with GLM as Meta Model

ensemble_Model = glm(Target ~ ., train_Pred_All_Models, family = binomial)
summary(ensemble_Model)

# Check the "ensemble_Model model" on the train data
ensemble_Train = predict(ensemble_Model, train_Pred_All_Models, 
                         type = "response")
ensemble_Train = ifelse(ensemble_Train > 0.5, 1, 0)
table(ensemble_Train)

cm_Ensemble = table(ensemble_Train, train_Pred_All_Models$Target)
sum(diag(cm_Ensemble))/sum(cm_Ensemble)

cat("Accuracy : ",sum(diag(cm_Ensemble))/sum(cm_Ensemble))
cat("Recall : ",cm_Ensemble[2,2]/sum(cm_Ensemble[2,]))
cat("precision : ", cm_Ensemble[2,2]/sum(cm_Ensemble[,2]))
cat("F1 Score : ", 2*(cm_Ensemble[2,2]/sum(cm_Ensemble[,2])*cm_Ensemble[2,2]/sum(cm_Ensemble[2,]))/((cm_Ensemble[2,2]/sum(cm_Ensemble[,2])+cm_Ensemble[2,2]/sum(cm_Ensemble[2,]))))