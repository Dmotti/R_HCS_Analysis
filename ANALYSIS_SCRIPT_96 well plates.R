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

#This is just a try
#Try number 2

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
GFP_pos_subset<-subset(GFP_subset, Plasmid != "Oxr1-2A-mCherry")
#GFP_subset<-subset(GFP_pos_subset, Plasmid != "S327A-2A-mCherry")


mCh_pos_subset_1<-subset(filtered_mCh, Plasmid == "S306A-2A-mCherry")
mCh_pos_subset_2<-subset(filtered_mCh, Plasmid == "Oxr1-2A-mCherry")
#mCh_pos_subset_2<-subset(filtered_mCh, Plasmid == "S327A-2A-mCherry")


#Putting together all this subset shoud be enough to have all the needed data!

filtered<-rbind(GFP_pos_subset,mCh_pos_subset_1, mCh_pos_subset_2,empty)

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


#How do I create a matrix that reports for each 
#plasmid conditions the average across each measured feature?
mat_vec<-c()

for(i in levels(filtered$Plasmid)) {
  cells_i<-subset(filtered, Plasmid == i)
  vec<-c(i, sapply(cells_i[,9:37], mean))
  mat_vec<-rbind(mat_vec, vec, deparse.level=0)
}
head(mat_vec)

write.table(mat_vec, "Output_Allfeatures_vectors.txt", row.names=FALSE, sep="\t")




#STOP HERE FOR 96-WELL PLATES!!!!!!!
#STOP HERE FOR 96-WELL PLATES!!!!!!!
#STOP HERE FOR 96-WELL PLATES!!!!!!!
#STOP HERE FOR 96-WELL PLATES!!!!!!!








#CALCULATING Z-SCORES - This part does not apply to this kind of plates since it's hard to get to 100 cells
#per well, even with the high transfection efficiency I have.
#But I'll leave it in, just don't use it

results_m<-data.frame(results)

head(results_m)

Controls<-subset(results_m, Condition == "?Luc_SIBR_GFP")
head(Controls)


#CREATE FUNCTION TO CALCULATE Z-SCORE
zscore<-function(x,y){
  (x - mean(as.numeric(paste(y)))) / sd(as.numeric(paste(y)))
}

#For loop to apply function to individual conditions
zscores<-c()
for (i in levels(results_m$Condition)){
  cond<-subset(results_m, Condition == i)
  meantemp<-cond$Mean
  a<-zscore(as.numeric(paste(meantemp)), Controls$Mean)
  zscores<-c(zscores,a)
}
print(zscores)

#Add Z-scores to the results table

results_fin<-cbind(results_m, zscores, deparse.level=0)
write.table(results_fin, "Output_Results_Z.txt", row.names=FALSE, sep="\t")
?rowMeans

#Make a list of results divided by Plasmid
#THE SPLIT FUNCTION IS VERY USEFUL!!!
#THE SPLIT FUNCTION IS VERY USEFUL!!!
#THE SPLIT FUNCTION IS VERY USEFUL!!!

splitresults <- split(results_fin, results_fin$Plasmid)
splitresults

#Create a function that make a vector with mean and strandard deviation
avgz <- function(x){
  a<-mean(x$zscores)
  b<-sd(x$zscores)
  v<-c(a,b)
}

#apply the function to the list of objects

avg_sd_list<-lapply(splitresults, avgz)
avg_sd_matrix<-do.call(rbind, avg_sd_list)
avg_sd_matrix
colnames(avg_sd_matrix) <-c("Mean", "Sd")
write.table(avg_sd_matrix, "Output_Zscores_results.txt", sep="\t")
#Nopl<-c("No Plasmid")
#avg_sd_noempty<-avg_sd_matrix[!rownames(avg_sd_matrix) %in% Nopl,]
#avg_sd_noempty

#THIS PART SETS THE PARAMETERS FOR THE GRAPHICAL OUTPUT
#UPPER AND LOWER BOUNDARIES OF THE GRAPHS,
#FIRST SET AS VECTORS TO BE FED TO THE BARPLOT FUNCTION!!!

upper<-c(avg_sd_matrix[1,1]+avg_sd_matrix[1,2], avg_sd_matrix[2,1]+avg_sd_matrix[2,2], avg_sd_matrix[3,1]+avg_sd_matrix[3,2],
         avg_sd_matrix[4,1]+avg_sd_matrix[4,2], avg_sd_matrix[5,1]+avg_sd_matrix[5,2], avg_sd_matrix[6,1]+avg_sd_matrix[6,2],
         avg_sd_matrix[7,1]+avg_sd_matrix[7,2], avg_sd_matrix[8,1]+avg_sd_matrix[8,2])

lower<-c(avg_sd_matrix[1,1]-avg_sd_matrix[1,2], avg_sd_matrix[2,1]-avg_sd_matrix[2,2], avg_sd_matrix[3,1]-avg_sd_matrix[3,2],
         avg_sd_matrix[4,1]-avg_sd_matrix[4,2], avg_sd_matrix[5,1]-avg_sd_matrix[5,2], avg_sd_matrix[6,1]-avg_sd_matrix[6,2],
         avg_sd_matrix[7,1]-avg_sd_matrix[7,2], avg_sd_matrix[8,1]-avg_sd_matrix[8,2])
pdf(file="Neurite Total Length.pdf", height = 8, width = 9)
par(mar=c(10, 4, 3, 2), font.axis=2, font.lab=2)
midpts<-barplot(avg_sd_matrix[,1], las=3, axes=FALSE, names.arg=FALSE, ylim=c((min(lower)-0.3), max(upper)),
                col=c("black", "white", "white", "grey", "grey", "grey", "grey", "grey"), 
                ylab="Z=Score", main = "Neurite Total Length")
axis(2)
axis(1, at=midpts, labels=rownames(avg_sd_matrix), las=3, cex.axis=0.8)
library(Hmisc)
errbar(midpts, avg_sd_matrix[,1], upper, lower, add=T, cex=0, lwd=1.5)
abline(h=2, col="red", lty=2, lwd=2)
abline(h=-2, col="red", lty=2, lwd=2)
dev.off()

tiff("Neurite Total Length.tiff", width=16, height=16, units="cm", res=600)
par(mar=c(10, 4, 3, 2), font.axis=2, font.lab=2)
midpts<-barplot(avg_sd_matrix[,1], las=3, axes=FALSE, names.arg=FALSE, ylim=c((min(lower)-0.3), max(upper)),
                col=c("black", "white", "white", "grey", "grey", "grey", "grey", "grey"), 
                ylab="Z=Score", main = "Neurite Total Length")
axis(2)
axis(1, at=midpts, labels=rownames(avg_sd_matrix), las=3, cex.axis=0.8)
library(Hmisc)
errbar(midpts, avg_sd_matrix[,1], upper, lower, add=T, cex=0, lwd=1.5)
abline(h=2, col="red", lty=2, lwd=2)
abline(h=-2, col="red", lty=2, lwd=2)
dev.off()






























#Calculate for WELL
results_1<-c()
results_1

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
  results_1<-rbind(results_1, stats, deparse.level=0)
}
head(results_1)
colnames(results_1)<-c("Condition","Mean", "StDev", "StdErr", "Min", "1st Qr", "Median", "3rd Qr", "Max", "Sum", "Count")
rownames(results_1)<- c(levels(filtered$Plasmid.1))
write.table(results_1, "Output_Results_well.txt", row.names=FALSE, sep="\t")



#EXTRA PART OF CODE FOR DENSITIES
#I'm wirting an extra part of the code to analyze the densities separately
#The first thing to do is to create a table of elements belonging to same densities
#starting from the filtered elements
density_1<-subset(filtered, Well %in% c("B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10",
                                        "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
                                        "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10")) 


#Calculates stats for each well
#Prepare the matrix to store stats
results<-c()
results


#Calculate for CONDITION
for(i in levels(density_1$Plasmid)) {
  cells<-subset(density_1, Plasmid == i)
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
rownames(results)<- c(levels(density_1$Plasmid))
write.table(results, "Output_Results_plasmid_Density_1.txt", row.names=FALSE, sep="\t")



#EXTRA PART OF CODE FOR DENSITIES
#And now we go for density 2
density_2<-subset(filtered, Well %in% c("E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10",
                                        "F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10",
                                        "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")) 


#Calculates stats for each well
#Prepare the matrix to store stats
results<-c()
results


#Calculate for CONDITION
for(i in levels(density_2$Plasmid)) {
  cells<-subset(density_2, Plasmid == i)
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
rownames(results)<- c(levels(density_2$Plasmid))
write.table(results, "Output_Results_plasmid_Density_2.txt", row.names=FALSE, sep="\t")


