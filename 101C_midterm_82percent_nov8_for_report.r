# Code for model with 82.298% accuracy
# Removed 3 near-zero variance predictors below:
# "Recurrent_missense_fraction"      "Lost_start_and_stop_fraction"     "Canyon_genebody_hypermethylation"
# SEED 12, SEED 5
# Bootstrap - 1000 resamples, 500 NG resamp, correlation cutoff = 0.75, partitioned data set 70/30
library(readr)
library(mclust)
library(MASS)
library(class)
library(caret)
library(MLeval)
library(MLmetrics)
library(tidyverse)

# Load data into RStudio
training <- read_csv("C:/Users/kevin/Documents/Classes/STATS 101C/ucla-stats101c-lec4/training.csv")
test <- read_csv("C:/Users/kevin/Documents/Classes/STATS 101C/ucla-stats101c-lec4/test.csv")

# Generating boxplots of class against each predictor, commented out for brevity:
# for ( i in temp[1:length(temp)]) {
#   bp = ggplot(data = training, 
#               aes_string(x = names(training)[99], 
#                          y = names(training)[i], 
#                          group = names(training)[99])) + geom_boxplot()
#   print(bp)
# }

# Variables of interest that should be kept (after examining the box plots)
temp <- c(3, 4, 5, 7, 8, 9, 13, 15, 16, 17, 18, 20, 22, 23, 25, 26, 28, 29, 31,
          35, 36, 37, 38, 39, 40, 41, 43, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54,
          55, 57, 58, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
          78, 79, 80, 83, 87, 88, 89, 91, 92, 93, 94, 95, 96, 97, 98, 99)

# Finding variables that are correlated with each other and removing them
y <- training[,temp]
# Temporarily exclude the response variable column (72 in this new set), 
# and look at strictly predictor-to-predictor correlations
y2 <- y[,-72] 
corr_mat_y2 <- cor(y2)
corr_y2 <- findCorrelation(corr_mat_y2, cutoff = 0.75)

set.seed(12)
# Attempt to equalize the class imbalance, randomly sample 500 values from NGs
NGs <- y %>% filter(class == 0)
OGs <- y %>% filter(class == 1)
TSGs <- y %>% filter(class == 2)
NG.samp <- sample(1:dim(NGs)[1], 500) 

balanced.data <- rbind(NGs[NG.samp, ], OGs, TSGs)
balanced.data$class <- factor(balanced.data$class)
levels(balanced.data$class) <- c("NG", "OG", "TSG")

# Remove highly-correlated variables from each other first from finished model
balanced.data.new <- balanced.data[, -corr_y2]

##############################################################################
# Remove predictors with extremely low, near-zero variance
# recurrent_missense_fraction was causing warning messages to surface when running train()
# So this motivated further examination of the other predictors' variances.

# Examine whether each predictor is deemed to have a near-zero variance
all_ps = nearZeroVar(balanced.data.new[, -ncol(balanced.data.new)], names = TRUE,
                     freqCut = 19, uniqueCut = 10, saveMetrics = TRUE) # for viewing purposes
print(all_ps)

# Identify and possibly exclude other nonzero predictors:
# NOTE: 2 and 20 for freqCut and uniqueCut are considered an aggressive approach. 
#       The default parameter setting is 19 and 10, which is fairly conservative.
other_nzp = nearZeroVar(balanced.data.new[, -ncol(balanced.data.new)], names = TRUE,
            freqCut = 19, uniqueCut = 10, saveMetrics = FALSE) 
print(other_nzp)

# Unselect the near-zero variance predictors:
balanced.data.new = balanced.data.new %>% select(-other_nzp)

##############################################################################
# Test and training set for balanced.data.new
set.seed(5)
trainIndex <- createDataPartition(balanced.data.new$class, p = 0.7, list = FALSE)
gene.train <- balanced.data.new[trainIndex,]
gene.test <- balanced.data.new[-trainIndex,]
# Bootstrap with resampling set to 1000
train_control <- trainControl(method = "boot", number = 1000, summaryFunction = multiClassSummary,
                              classProbs = TRUE, savePrediction = TRUE)
# Multinomial logistic regression model. The LRfit line may take a while to load due to 1000 resamples
LRfit <- train(class ~ ., data = gene.train, method = "multinom", preProc = c("center", "scale"),
               trControl = train_control, trace = FALSE, metric = "Mean_Sensitivity")
print(LRfit)
# Predict and test out the multinomial logistic regression model
predLR <- predict(LRfit, newdata = gene.test)
confusionMatrix(data = predLR, reference = gene.test$class)

##############################################################################
## PREDICTION ON TEST DATA AND PROCESSING FOR SUBMISSION
saved_prediction_copy = read_csv("C:/Users/kevin/Documents/boot_predict_1000_with_3_nzvpredictors_removed_nov8_v2.csv")
# NOTE: saved_prediction_copy is a copy of the csv that obtained the 82.298% accuracy score.
# Comparison will be run against the csv we generate here to make sure the results are reproducible
test2 <- test[,temp]
test2 <- test2[,-corr_y2]

# Remove the 3 near-zero variance predictors
test2 <- test2 %>% select(-other_nzp) %>% as.data.frame()

# Apply multinomial regression model to the test data
predLR.test <- predict(LRfit, newdata = test2)

# Preparing and assembling the prediction data frame
prediction <- data.frame(test$id, predLR.test)
colnames(prediction) <- c("id", "class")
levels(prediction$class) <- c(0, 1, 2)

# CHECK: which(...) should return integer(0), and length(which(...)) should return 0
# meaning that the predictions are exactly the same for each observation in the test data
which(saved_prediction_copy$class != prediction$class)
length(which(saved_prediction_copy$class != prediction$class))

# Store file
write_csv(prediction, "boot_predict_1000_with_3_nzvpredictors_removed_nov8_final.csv")

# Print out the predictors used in the model
predictors_used = names(test2)
predictors_used = subset(predictors_used, predictors_used != "class")
print(predictors_used)
