##using all data, without removing row level duplicates

##Remove all environment variables, to relinguis the RAM space and also to avoid overlaps...

rm(list=ls(all=TRUE))

pkgs <- c("readr",
#          "spam","LICORS",
          "cluster","vegan","sqldf","C50")
install.packages(pkgs)

library(readr)
# library(spam)
# library(LICORS)
library(vegan)
library(cluster)
#Try usig this later all in the end
library(sqldf)
library(RSQLite)
library(tcltk)
library(ROCR)
library(caret)
library(C50)
setwd("E:/Workspace/Project")

getwd()

#Loading data into a DataFrame using read_csv, as to apply format to the date variable while loading itself.
#as I have abserved that the date is having diffrent format in the data set  provided.
PriorAuth_Data <- read_csv("PriorAuth_Data.csv",
                           col_types = cols(TransDate = col_date(format = "%m/%d/%Y")))


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

nrow(PriorAuth_Data)

attributes <- setdiff(names(PriorAuth_Data),c('Target','TransDate'))

rm(list= ls()[!(ls() %in% c('PriorAuth_Data','attributes'))])

##Create frequency table for each individual attribute against the Target, and standardising the freque
##to find the clusters, to reduce the huge levels in the data.

select <- c('Select ')
case <- ', sum(case when Target = 1 then 1 else 0 end) as True, sum(case when Target = 0 then 1 else 0 end) as False from PriorAuth_Data group by '
decoSelect <- 'select True,false from '

for (i in 1:length(attributes)){

  ##Creating Frequency tables for all Individual attributes against the Class variable dynamically
  assign(attributes[i],sqldf(paste(select,attributes[i],case,attributes[i],sep = '')))

  ##Standardizing the True False Frequency of all newly created dataframes, and make it ready for
  assign(attributes[i],cbind(sqldf(paste(select,attributes[i],'from',attributes[i])),decostand(sqldf(paste(decoSelect,attributes[i])),method = 'range')))

}

rm(select,case,decoSelect,i,Bin1)

## Running kmeans++ for the k-values ranging 1:10 to Identifying the optimal k for all the Indipendent variable
Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(Drug[-1],centers = j)$withinss) } # best K = 4
plot(Clusterz)
clusplot(Drug[-1],kmeans(Drug[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
Drug <- cbind(Drug[1],kmeans(Drug [-1],centers = 4)$cluster)
colnames(Drug) <-c('Drug','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j] <- sum(kmeans(DrugSubClass[-1],centers = j)$withinss) } # best K = 3
plot(Clusterz)
clusplot(DrugSubClass[-1],kmeans(DrugSubClass[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
DrugSubClass <- cbind(DrugSubClass[1],kmeans(DrugSubClass [-1],centers = 3)$cluster)
colnames(DrugSubClass) <-c('DrugSubClass','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(DrugClass[-1],centers = j)$withinss)  } # best K = 3
plot(Clusterz)
clusplot(DrugClass[-1],kmeans(DrugClass[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
DrugClass <- cbind(DrugClass[1],kmeans(DrugClass[-1],centers = 4)$cluster)
colnames(DrugClass) <-c('DrugClass','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(Drug_Chemical_Name[-1],centers = j)$withinss)  }# best K = 4
plot(Clusterz)
clusplot(Drug_Chemical_Name[-1],kmeans(Drug_Chemical_Name[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
Drug_Chemical_Name <- cbind(Drug_Chemical_Name[1],kmeans(Drug_Chemical_Name[-1],centers = 4)$cluster)
colnames(Drug_Chemical_Name) <-c('Drug_Chemical_Name','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(GPI[-1],centers = j)$withinss)  }# best K = 4
plot(Clusterz)
clusplot(GPI[-1],kmeans(GPI[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
GPI <- cbind(GPI[1],kmeans(GPI[-1],centers = 4)$cluster)
colnames(GPI) <-c('GPI','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(NDC[-1],centers = j)$withinss)  }# best K = 5
plot(Clusterz)
clusplot(NDC[-1],kmeans(NDC[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
NDC <- cbind(NDC[1],kmeans(NDC[-1],centers = 5)$cluster)
colnames(NDC) <-c('NDC','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(DrugGroup[-1],centers = j)$withinss)  }# best K = 2
plot(Clusterz)
clusplot(DrugGroup[-1],kmeans(DrugGroup[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
DrugGroup <- cbind(DrugGroup[1],kmeans(DrugGroup[-1],centers = 3)$cluster)
colnames(DrugGroup) <-c('DrugGroup','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(DoctorID[-1],centers = j)$withinss)  }# best K = 5
plot(Clusterz)
clusplot(DoctorID[-1],kmeans(DoctorID[-1],centers = 6)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
DoctorID <- cbind(DoctorID[1],kmeans(DoctorID[-1],centers = 6)$cluster)
colnames(DoctorID) <-c('DoctorID','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(RxGroupId[-1],centers = j)$withinss)  }# best K = 4
plot(Clusterz)
clusplot(RxGroupId[-1],kmeans(RxGroupId[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
RxGroupId <- cbind(RxGroupId[1],kmeans(RxGroupId[-1],centers = 4)$cluster)
colnames(RxGroupId) <-c('RxGroupId','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(Bin[-1],centers = j)$withinss)  }# best K = 2
plot(Clusterz)
clusplot(Bin[-1],kmeans(Bin[-1],centers = 2)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
Bin <- cbind(Bin[1],kmeans(Bin[-1],centers = 2)$cluster)
colnames(Bin) <-c('Bin','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(PCN[-1],centers = j)$withinss)  }# best K = 3
plot(Clusterz)
clusplot(PCN[-1],kmeans(PCN[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
PCN <- cbind(PCN[1],kmeans(PCN[-1],centers = 3)$cluster)
colnames(PCN) <-c('PCN','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz [j]  <-sum(kmeans(State [-1],centers = j)$withinss)  }# best K = 3
plot(Clusterz)
clusplot(State[-1],kmeans(State[-1],centers = 2)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
State <- cbind(State[1],kmeans(State [-1],centers = 2)$cluster)
colnames(State) <-c('State','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(comb1[-1],centers = j)$withinss)  }# best K = 5
plot(Clusterz)
clusplot(comb1[-1],kmeans(comb1[-1],centers = 4)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
comb1 <- cbind(comb1[1],kmeans(comb1[-1],centers = 5)$cluster)
colnames(comb1) <-c('comb1','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(comb2[-1],centers = j)$withinss)  }# best K = 3
plot(Clusterz)
clusplot(comb2[-1],kmeans(comb2[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=3,lines=0)
comb2 <- cbind(comb2[1],kmeans(comb2[-1],centers = 3)$cluster)
colnames(comb2) <-c('comb2','cluster')

Clusterz <- 0
for (j in 1:10) {  Clusterz[j]  <-sum(kmeans(comb3[-1],centers = j)$withinss)  }# best K = 3
plot(Clusterz)
clusplot(comb3[-1],kmeans(comb3[-1],centers = 3)$cluster,color=TRUE,shade=TRUE,lables=2,lines=0)

comb3 <- cbind(comb3[1],kmeans(comb3[-1],centers = 3)$cluster)
colnames(comb3) <-c('comb3','cluster')

rm(j,UserID,Clusterz)

PriorAuth_num_Data <- sqldf('select d.cluster DrugId,sc.cluster DrugSubClassId , dc.cluster DrugClassId ,c.cluster Drug_ChemicalId, g.cluster GpiId , n.cluster NdcId , dg.cluster DrugGroupId,
       Dr.cluster DoctorID , gi.cluster RxGroupId , b.cluster BinId, p.cluster PcnId, s.cluster StateId, c1.cluster com1Id, c2.cluster com2Id ,c3.cluster com3Id , PA.Target,PA.TransDate
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

# Below is not used for factorised data
PriorAuth_num_Data$DoctorID<- as.integer(PriorAuth_num_Data$DoctorID)
PriorAuth_num_Data$RxGroupId <- as.integer(PriorAuth_num_Data$RxGroupId)
# till here

attributes <- setdiff(names(PriorAuth_num_Data),c('TransDate'))
PriorAuth_fact_Data <- PriorAuth_num_Data
PriorAuth_fact_Data[attributes] <- lapply(PriorAuth_fact_Data[attributes], factor)

str(PriorAuth_fact_Data)

# Below is not used for factorised data
cor(PriorAuth_num_Data[setdiff(names(PriorAuth_num_Data),c('TransDate'))])

#Split Train and test using Time-based sampling approach
# Total 18 days data
# data grouped by date info
data.frame(table(PriorAuth_num_Data$TransDate))

traindata <- subset(PriorAuth_num_Data, TransDate<=as.Date("2013-08-04"))
nrow(traindata) # 7255 rows for training i.e 77%
testdata <- subset(PriorAuth_num_Data, TransDate>as.Date("2013-08-04"))
nrow(testdata)#2138 rows for testing teh model i.e 23%
str(traindata)

# Building Modelz
###################################################################################################
# Logistic regression
###################################################################################################
LogReg <- glm(Target ~ ., data=traindata, family=binomial)
summary(LogReg)

step(glm(Target ~ .,data=traindata),direction = 'backward')

LogReg <- glm(formula = Target ~ DrugId + DrugSubClassId + NdcId + DrugGroupId + 
                DoctorID + RxGroupId + PcnId + com3Id + TransDate , data = traindata)

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
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/
                      ((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,])))
    )

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
plot(prf,colorize = TRUE, print.cutoffs.at=seq(0,1,by=0.1), text.adj=c(-0.2,1.7),main='AUC Value .69')

abline(a=0, b= 1)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc 

###################################################################################################
#Decision Tree
###################################################################################################
rm(list= ls()[!(ls() %in% c('traindata','testdata'))])

names(traindata)
PA_rpart <- rpart(Target ~ ., data=traindata, method="class")
plot(PA_rpart,main="Classification Tree for loan Class",margin=0.15,uniform=TRUE)
text(PA_rpart,use.n=T)
summary(PA_rpart)

conf.mat = table(traindata$Target, predict(PA_rpart, newdata=traindata, type="class"))
(conf.mat[2,2])/(conf.mat[2,1]+conf.mat[2,2])*100

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

conf.mat = table(testdata$Target, predict(PA_rpart, newdata=testdata, type="class"))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))


######### Trying C50 for good luck

traindata1<- traindata[-17]
str(traindata1)

PA_C50<- C5.0(Target~NdcId + com3Id +DoctorID,data = traindata1, rules=TRUE)
summary(PA_C50)
C5imp(PA_C50, pct=TRUE)


conf.mat = table(traindata1$Target, predict(PA_C50, newdata=traindata1, type="class"))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

conf.mat = table(testdata$Target, predict(PA_C50, newdata=testdata, type="class"))

cat("Accuracy : ",sum(diag(conf.mat))/sum(conf.mat))
cat("Recall : ",conf.mat[2,2]/sum(conf.mat[2,]))
cat("precision : ", conf.mat[2,2]/sum(conf.mat[,2]))
cat("F1 Score : ", 2*(conf.mat[2,2]/sum(conf.mat[,2])*conf.mat[2,2]/sum(conf.mat[2,]))/((conf.mat[2,2]/sum(conf.mat[,2])+conf.mat[2,2]/sum(conf.mat[2,]))))

