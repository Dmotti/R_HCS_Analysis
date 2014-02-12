#R script for ANALYSIS of 96 wells format HIGH CONTENT ANALYSIS

#Author: Dario Motti 
#Started 11/05/2013

#I'm still working on overexpressing microRNAs in different neuronal types,
#I've been trying for a long time with 48 wells plates, but it seems like 
#they had a lot of problems specially when I turned to the hippocampal cells
#and they were all dying.
#So I'm turning to a 96 wells plates with an easier layout to kjeep the media
#outof the incubator the least time possible.
#But to do that I need to rewrite my analysis code.


#Change directory to the 96 wells format:
setwd("C:/Users/dmotti/Desktop/R/96-wells format analysis")
#Start reading the data
Celldetails<-read.table("data.txt", sep="\t", header=TRUE)
head(Celldetails)

#Now we upload the matrix, that it's easier this time
Matrix<-read.table("matrix.txt", sep="\t", header=TRUE)
head(Matrix)

#Merging Matrix and Data to assign conditions to the wells
head(Matrix)
Datapre<-merge(Celldetails,Matrix, by.x="Well", by.y="Well", all.x=TRUE)
write.table(Datapre, "Output_Data.txt", sep="\t")

#Starting filtering only cells with neurites (NeuriteTotalLength != 0)
Data<-subset(Datapre, NeuriteTotalLengthCh2 != 0)
head(Data)

levels(Data$Plasmid)

#Read in Valid Neuron Count
Valid<-read.table("valid.txt", header=FALSE, sep="\t")
colnames(Valid) <- c(3:10)
rownames(Valid) <- c(paste(LETTERS[2:7]))
write.table(Valid, "Output_Counts.txt", sep="\t")

#CALCULATES NUMBER OF CELLS PER WELL WITH NEURITES

#create an empty vector
Grown<-c()
#Fill the vector with number of cells positive per well
for (i in levels(Data$Well)) {
  cells<-subset(Data, Well == i)
  Grown<-c(Grown,nrow(cells))
}
#Change the vector into a matrix, R designs matrices
#rows first, so I will have to transpose afterwords


#For plates 2 and 4
Grown<-append(Grown, c(0), after=(0))
Grown<-append(Grown, c(0,0,0), after=(6))
Grown<-append(Grown, c(0,0,0), after=(14))
Grown<-append(Grown, c(0,0,0), after=(22))
Grown<-append(Grown, c(0,0,0), after=(30))
Grown<-append(Grown, c(0,0,0), after=(38))
Grown<-append(Grown, c(0,0))





#USED ONLY WHEN YOU HAVE SOME 0s in Valid Neuron count
#Grown<-append(Grown, c(0,0))
#Grown<-append(Grown, 0, after=(40))
dim(Grown)<-c(8,6)
Grown<-t(Grown)

#Nedd to move the first column at the end because it's 
#actually B10, only than when reordered alphabetically it got moved
#at the beginning.
Grown<-Grown[,c(2,3,4,5,6,7,8,1)]

colnames(Grown) <- c(3:10)
rownames(Grown) <- c(paste(LETTERS[2:7]))
print(Grown)
write.table(Grown, file="grown_neurons.txt", sep="\t")
file.append("Output_Counts.txt", "grown_neurons.txt")
file.remove("grown_neurons.txt")

#Isolating cells from No Plasmid condition
empty<-subset(Data, Plasmid == "No Plasmid")
head(empty)

# Calculates 99% percentile of GFP intensity on the  untransfected cells
GFPthresh<-quantile(empty$AvgIntenCh4, c(.99))
GFPthresh
write.table(GFPthresh, "GFPintensity.txt", sep="\t")
file.append("Output_Counts.txt", "GFPintensity.txt")
file.remove("GFPintensity.txt")



#Filter cells for GFP intensity
filtered_GFP<-subset(Data, AvgIntenCh4 > GFPthresh)
head (filtered_GFP)

# Separates rows by Well

#create an empty vector
Positive_GFP<-c()
#Fill the vector with number of cells positive for GFPper well
for (i in levels(filtered_GFP$Well)) {
  cells<-subset(filtered_GFP, Well == i)
  Positive_GFP<-c(Positive_GFP,nrow(cells))
}
#Change the vector into a matrix, R designs matrices
#rows first, so I will have to transpose afterwords




#AGAIN ONLY FOR THE 0s

Positive_GFP<-append(Positive_GFP, c(0), after=(0))
Positive_GFP<-append(Positive_GFP, c(0,0,0), after=(6))
Positive_GFP<-append(Positive_GFP, c(0,0,0), after=(14))
Positive_GFP<-append(Positive_GFP, c(0,0,0), after=(22))
Positive_GFP<-append(Positive_GFP, c(0,0,0), after=(30))
Positive_GFP<-append(Positive_GFP, c(0,0,0), after=(38))
Positive_GFP<-append(Positive_GFP, c(0,0))


#Positive_GFP<-append(Positive_GFP, c(0,0))
#Positive_GFP<-append(Positive_GFP, 0, after=(40))


dim(Positive_GFP)<-c(8,6)
Positive_GFP<-t(Positive_GFP)
Positive_GFP
#B10 is moved at the beginning becuase of alphabetical order,
#I need to move the column of B10 at the end again:

Positive_GFP<-Positive_GFP[,c(2,3,4,5,6,7,8,1)]
Positive_GFP

colnames(Positive_GFP) <- c(3:10)
rownames(Positive_GFP) <- c(paste(LETTERS[2:7]))
print(Positive_GFP)
write.table(Positive_GFP, file="positive_neurons_GFP.txt", sep="\t")
file.append("Output_Counts.txt", "positive_neurons_GFP.txt")
file.remove("positive_neurons_GFP.txt")


#Calculate % of transfection
TransfPer<- (Positive_GFP/Grown) * 100

write.table(TransfPer, "Transfection.txt", sep="\t")
file.append("Output_Counts.txt", "Transfection.txt")
file.remove("Transfection.txt")



#NOW I HAVE TO DO THE SAME THING FOR mCHERRY!!!

# Calculates 99% percentile of mCherry intensity on the  untransfected cells
mCh_thresh<-quantile(empty$AvgIntenCh3, c(.99))
mCh_thresh
write.table(mCh_thresh, "mCherryintensity.txt", sep="\t")
file.append("Output_Counts.txt", "mCherryintensity.txt")
file.remove("mCherryintensity.txt")

#Filter cells for mCherry intensity
filtered_mCh<-subset(Data, AvgIntenCh3 > mCh_thresh)
head (filtered_mCh)

# Separates rows by Well

#create an empty vector
Positive_mCh<-c()
#Fill the vector with number of cells positive for mCherry per well
for (i in levels(filtered_mCh$Well)) {
  cells<-subset(filtered_mCh, Well == i)
  Positive_mCh<-c(Positive_mCh,nrow(cells))
}
#Change the vector into a matrix, R designs matrices
#rows first, so I will have to transpose afterwords

Positive_mCh<-append(Positive_mCh, c(0), after=(0))
Positive_mCh<-append(Positive_mCh, c(0,0,0), after=(6))
Positive_mCh<-append(Positive_mCh, c(0,0,0), after=(14))
Positive_mCh<-append(Positive_mCh, c(0,0,0), after=(22))
Positive_mCh<-append(Positive_mCh, c(0,0,0), after=(30))
Positive_mCh<-append(Positive_mCh, c(0,0,0), after=(38))
Positive_mCh<-append(Positive_mCh, c(0,0))


#Positive_mCh<-append(Positive_mCh, c(0,0))
#Positive_mCh<-append(Positive_mCh, 0, after=(40))


dim(Positive_mCh)<-c(8,6)
Positive_mCh<-t(Positive_mCh)
Positive_mCh
#B10 is moved at the beginning becuase of alphabetical order,
#I need to move the column of B10 at the end again:
Positive_mCh<-Positive_mCh[,c(2,3,4,5,6,7,8,1)]
Positive_mCh

colnames(Positive_mCh) <- c(3:10)
rownames(Positive_mCh) <- c(paste(LETTERS[2:7]))
print(Positive_mCh)
write.table(Positive_mCh, file="positive_neurons_mCh.txt", sep="\t")
file.append("Output_Counts.txt", "positive_neurons_mCh.txt")
file.remove("positive_neurons_mCh.txt")


#Calculate % of transfection
TransfPer<- (Positive_mCh/Grown) * 100

write.table(TransfPer, "Transfection.txt", sep="\t")
file.append("Output_Counts.txt", "Transfection.txt")
file.remove("Transfection.txt")


rm(Datapre)
rm(cells)
remove(Positive_GFP)
remove(Positive_mCh)
remove(Valid)
remove(TransfPer)
remove(i)
remove(GFPthresh)
remove(mCh_thresh)
remove(Grown)

#SECOND PART
#I need to recreate the Data file containing only the GFP positive cells in the
#wells transfected with EGFP, the mChery cells only from the wells transfected 
#with mCherry expressing vectors and the empty cells

#Create function for standard error

stderr<-function(x) sd(x)/sqrt(length(x))

#Let's start subsetting the GFP
GFP_subset<-subset(filtered_GFP, Plasmid != "S306A-2A-mCherry")
#GFP_pos_subset<-subset(GFP_subset, Plasmid != "Oxr1-2A-mCherry")
#GFP_subset<-subset(GFP_pos_subset, Plasmid != "S327A-2A-mCherry")


mCh_pos_subset_1<-subset(filtered_mCh, Plasmid == "S306A-2A-mCherry")
#mCh_pos_subset_2<-subset(filtered_mCh, Plasmid == "Oxr1-2A-mCherry")
#mCh_pos_subset_2<-subset(filtered_mCh, Plasmid == "S327A-2A-mCherry")


#Putting together all this subset shoud be enough to have all the needed data!

filtered<-rbind(GFP_subset,mCh_pos_subset_1,empty)

#Calculates stats for each well
#Prepare the matrix to store stats
results<-c()
results


#Calculate for CONDITION
for(i in levels(filtered$Plasmid)) {
  cells<-subset(filtered, Plasmid == i)
  stats<-c(i, mean(cells$NeuriteTotalLengthCh2),
           sd(cells$NeuriteTotalLengthCh2),
           stderr (cells$NeuriteTotalLengthCh2),
           min(cells$NeuriteTotalLengthCh2),
           quantile(cells$NeuriteTotalLengthCh2, c(.25))[["25%"]],
           median(cells$NeuriteTotalLengthCh2), 
           quantile(cells$NeuriteTotalLengthCh2, c(.75))[["75%"]],
           max(cells$NeuriteTotalLengthCh2), 
           sum(cells$NeuriteTotalLengthCh2),
           nrow(cells))
  results<-rbind(results, stats, deparse.level=0)
}
head(results)
colnames(results)<-c("Condition","Mean", "StDev", "StdErr", "Min", "1st Qr", "Median", "3rd Qr", "Max", "Sum", "Count")
rownames(results)<- c(levels(filtered$Plasmid))
write.table(results, "Output_Results_plasmid.txt", row.names=FALSE, sep="\t")



results_2<-c()
results_2


#Calculate for CONDITION
for(i in levels(filtered$Plasmid.1)) {
  cells<-subset(filtered, Plasmid.1 == i)
  stats<-c(i, mean(cells$NeuriteTotalLengthCh2),
           sd(cells$NeuriteTotalLengthCh2),
           stderr (cells$NeuriteTotalLengthCh2),
           min(cells$NeuriteTotalLengthCh2),
           quantile(cells$NeuriteTotalLengthCh2, c(.25))[["25%"]],
           median(cells$NeuriteTotalLengthCh2), 
           quantile(cells$NeuriteTotalLengthCh2, c(.75))[["75%"]],
           max(cells$NeuriteTotalLengthCh2), 
           sum(cells$NeuriteTotalLengthCh2),
           nrow(cells))
  results_2<-rbind(results_2, stats, deparse.level=0)
}
head(results_2)
colnames(results_2)<-c("Condition","Mean", "StDev", "StdErr", "Min", "1st Qr", "Median", "3rd Qr", "Max", "Sum", "Count")
rownames(results_2)<- c(levels(filtered$Plasmid.1))
write.table(results_2, "Output_Results_density.txt", row.names=FALSE, sep="\t")












