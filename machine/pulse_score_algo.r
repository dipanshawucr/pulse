library(randomForest)

MM = 145  # Number of training examples
M = MM
N = 15929  # MM + Unlabelled data

# The number of iterations, different iterations are not dependent as this is not an iterative algorithm
numIterations = 5

# This should be adjusted using validation data, 0.6 works for our data set
cutoff = 0.6 # cutoff threshold for considering an isoform positive or negative

#directory containing the feature files
#EDIT UNDER NEEDS
setwd(".")


####===================================================================================
#    Initialization
####===================================================================================

Y = cbind(rep(0, N))

#This is the feature file with headers (postitive training data is the first M rows, except the header)
#This file is used to obtain a machine learning model (even the unlabelled data is part of this model as required by the semi-supervised learning set up)
#Each row represents an isoform and its features
#EDIT UNDER NEEDS
features = read.table('/home/wonjunetai/Documents/KimLab/pulse/input/param_files/ML_Training_data.txt', sep=',', header=TRUE)
N = dim(Y)[1]
print(as.character(N))

#define positive and negative sets here
pset = features[1:M,]
nset = features[(M+1):N,]

numPositive = length(pset[,1])
for (i in 1:numPositive){
	Y[i] = 1
}
for (i in (numPositive+1):N){
	Y[i] = 0
}

X = rbind(pset, nset)

d = cbind(Y, X)
names(d)[1] = "class"
d$class=as.factor(d$class)

# this segment of the code simply obtains a matrix of features (N by X+1) where X is the number of features,
# N is the number of rows (isoforms)
# The extra column contains 0 or 1 where 1 means it's training data, 0 means it's not.  This matrix is stored in d.
###===================================================================

# This matrix will hold extra data feature file, be sure the name is correct, the features should be the same
# features as the model generation feature file (same feature order, same number of features etc.)
# each row is an isoform with features
# EDIT UNDER NEEDS
extraData = read.table('../output/features/G41726.MCF7.5.bam/features_headers.txt', sep = ',', header = TRUE)
extraScores = rep(0, length(extraData[,1])) # This variable holds a vector of scores for these isoforms.
# This matrix is typically the validation set, or some set you want to explore using the model.
# If it is the validation set, be sure to exclude duplicates of model's isoforms

# this contains a vector of scores for the model's isoforms
grandVotes = rep(0, length(d[,1]))

# This initializes scores for feature importance, it's a length X vector (X = number of features) and will hold a decimal
# value for their relative importance
importance = rep(0, length(d[1,])-1)
names(importance) = names(d)[2:length(d[1,])] # This is the names of the features for the feature importance vector

iterationGrandVotes = 0


#This loop is the actual algorithm
for (j in 1:numIterations){
	train=d
	print("STARTING ITERATIONS..")
	print(j)
	print("first classifier...")

	#train first classifier
	rfFit <- randomForest(class ~ . , data=train, ntree=2000, importance=TRUE)

	#obtain initial scores for data points (out of bag score, see randomForest documentation)
	print("reweighing...")
	trainPredictions = rfFit$votes[,2]

	#compute c value
	cValue = sum(trainPredictions[1:MM])/sum(trainPredictions)
	print(cValue)

	for (i in 1:N){
		label = 0
		if (i <= MM){
			label = 1
		} else {
			probability = trainPredictions[i] # this is P(s=1 | x)

			#reformed probability for first classifier
			mProb = ((1-cValue)/cValue)*(probability/(1-probability))

			#relabel unlabelled data based on probability
			if (runif(1, min=0, max=1) <= mProb){
				label = 1
			} else {
				label = 0
			}
		}
		train[i,1] = label
	}

	print(train)
	# train second classifier using new labels
	print("second classifier...")
	secondShot <- randomForest(class ~ . , data=train, ntree=2000, importance=TRUE)

	# Second shot is the final model we get for this iteration
	# IMPORTANT, this is the place where we obtain prediction scores for extra data, we add it to extra scores so
	# we can average over the number of iterations at the end
	extraScores <- extraScores + predict(secondShot, extraData, type="vote")[,2]
	print("done...")
	#here we generate scores for isoforms in the model (we use the out of bag scores, see randomForest documentation)
	grandVotes <- grandVotes + secondShot$votes[,2]
	#this line obtains the feature importances of this iteration
	importance <- importance + secondShot$importance[,3]
}

#average scores over iterations
grandVotes = grandVotes / numIterations
importance = importance / numIterations

#These are the prediction scores for extra data features file averaged over numIterations
extraScores <- extraScores / numIterations
#write the scores
write(extraScores, "../output/machine/G41726.MCF7.5.bam/PULSE_output.txt", sep="\n")
#These three score vectors is essentially the output

print("done")