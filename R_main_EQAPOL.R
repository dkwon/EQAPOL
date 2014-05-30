source('R_functions.R')
# Hierarchical DPGMM for each site, so common components (mus and sigmas, 14-dimensional) across samples within each site
# and only pi's (proportions) are different across samples
# For each site, cell subsets of interest were found by gating the common components (not the individual cells). 
# Analysis strategy: build a training data set by pooling data
# in the form
# Site; Subset; mu.vector 
# check if normalization is necessary 
# then perform classification for 
# 1. input: component centers (mus)
# 2. output: whether it belongs to a given cell subset

#####################################
### 
#####################################
############################################################################################
#### Perform normalization and affine transformation
#### Run this once and for all
#### Once this is done and transformed data sets are stored, directly go to IMPORT DATA below
############################################################################################
# repeat normalization for all files except for parameters, parameter maps, ungated (=no cell subset information) clusters: 
fnames = list.files("./EP4_orig/")
fnames = fnames[!grepl("(parameter)|(^[[:digit:]]{3}_mus.txt)", fnames)]
try.read = lapply(fnames, function(x) try(read.table(paste("./EP4_orig/", x, sep = "")), silent = T))
empty = sapply(try.read, function(x) inherits(x, "try-error")) #identifies the elements of the list with the try-error i.e. empty files 
fnames = fnames[!empty] #removes the elements identified above
invisible( # normalization for each site, each channel using 8 peak medians from 8-peak beads
lapply(fnames, norm.li, path.from = "EP4_orig/", path.peaks = "peaks_ep4/", path.to = "EP4_8peak/")
)
invisible( # normalization for each site, each channel using low/high peak medians from compensation beads
lapply(fnames, norm.li.comp, path.from = "EP4_orig/", data.peaks = comp_beads.EP4, path.to = "EP4_comp/")
)
invisible( # affine transformation on raw data for each site with optimal parameters 
lapply(fnames, affine, data.type = "orig", path.from = "EP4_orig/", fnames_affine = fnames_affine.EP4, path.to = "EP4_orig_affine/")
)
invisible( # affine transformation on comp-beads-normalized data for each site with optimal parameters (best affine transformation = least KL divergence)
lapply(fnames, affine, data.type = "comp", path.from = "EP4_comp/", fnames_affine = fnames_affine.EP4, path.to = "EP4_comp_affine/")
)
# repeat for EP5: 
## but first, rename original file names: 
## fnames.orig = list.files("./EP5_orig")
## fnames.orig = paste("./EP5_orig/", gsub("(EP5_)|(Time vs FSC::(.*::)*)", "", gsub("params", "parameters", gsub("FSC vs Time", "Time vs FSC", gsub("VS", "vs", gsub("IFNg2", "IFNg", fnames.orig))))), sep = "")
## file.rename(list.files("./EP5_orig", full.names = T), fnames.orig)
## fnames_peaks.orig = list.files("./peaks_ep5")
## fnames_peaks.orig = paste("./peaks_ep5/", gsub("EP5_", "", fnames_peaks.orig), sep = "")
## file.rename(list.files("./peaks_ep5", full.names = T), fnames_peaks.orig)
fnames = list.files("./EP5_orig")
fnames = fnames[!grepl("(parameters)|(^[[:digit:]]{3}_mus.txt)|(ids.txt)", fnames)]
try.read = lapply(fnames, function(x) try(read.table(paste("./EP5_orig/", x, sep = "")), silent = T))
empty = sapply(try.read, function(x) inherits(x, "try-error")) #identifies the elements of the list with the try-error i.e. empty files 
fnames = fnames[!empty] #removes the elements identified above
invisible(
lapply(fnames[!grepl("043", fnames)], norm.li, path.from = "EP5_orig/", path.peaks = "peaks_ep5/", path.to = "EP5_8peak/")
) # peaks_ep5/peak_medians_043 is missing
invisible(
lapply(fnames, norm.li.comp, path.from = "EP5_orig/", data.peaks = comp_beads.EP5, path.to = "EP5_comp/")
)
invisible(
lapply(fnames, affine, data.type = "orig", path.from = "EP5_orig/", fnames_affine = fnames_affine.EP5, path.to = "EP5_orig_affine/")
)
invisible(
lapply(fnames, affine, data.type = "comp", path.from = "EP5_comp/", fnames_affine = fnames_affine.EP5, path.to = "EP5_comp_affine/")
)

####### repeat for EP4 new: 
## but first, rename original file names: 
## fnames.orig = list.files("./EP4b_orig")
## fnames.orig = paste("./EP4b_orig/", gsub("(EP4_)|(Time vs FSC::(.*::)*)", "", gsub("params", "parameters", fnames.orig)), sep = "")
## file.rename(list.files("./EP4b_orig", full.names = T), fnames.orig)
fnames = list.files("./EP4b_orig")
fnames = fnames[!grepl("(parameters)|(^[[:digit:]]{3}_mus.txt)|(ids.txt)", fnames)]
try.read = lapply(fnames, function(x) try(read.table(paste("./EP4b_orig/", x, sep = "")), silent = T))
empty = sapply(try.read, function(x) inherits(x, "try-error")) #identifies the elements of the list with the try-error i.e. empty files 
fnames = fnames[!empty] #removes the elements identified above
invisible(
lapply(fnames, norm.li, path.from = "EP4b_orig/", path.peaks = "peaks_ep4/", path.to = "EP4b_8peak/")
) 
invisible(
lapply(fnames, norm.li.comp, path.from = "EP4b_orig/", data.peaks = comp_beads.EP4, path.to = "EP4b_comp/")
)
invisible(
lapply(fnames, affine, data.type = "orig", path.from = "EP4b_orig/", fnames_affine = fnames_affine.EP4, path.to = "EP4b_orig_affine/")
)
invisible(
lapply(fnames, affine, data.type = "comp", path.from = "EP4b_comp/", fnames_affine = fnames_affine.EP4, path.to = "EP4b_comp_affine/")
)

############################################
#### Import data: 
#### data.frame 'Mus_orig' (raw)
#### data.frame 'Mus_8peak' (normalized using 8-peak beads): NOT USED for affine transformation
#### data.frame 'Mus_comp' (normalized using compensation beads)
#### data.frame 'Mus_orig_affine' (affine transformed)
#### data.frame 'Mus_comp_affine' (normalized with comp beads THEN affine transformed)
############################################
Mus_orig.EP4 = import_mus('EP4_orig')
Mus_comp.EP4 = import_mus('EP4_comp')
#Mus_8peak.EP4 = import_mus('EP4_8peak')
#Mus_orig_affine.EP4 = import_mus('EP4_orig_affine')
Mus_comp_affine.EP4 = import_mus('EP4_comp_affine')

Mus_orig.EP5 = import_mus('EP5_orig')
Mus_comp.EP5 = import_mus('EP5_comp')
#Mus_8peak.EP5 = import_mus('EP5_8peak')
#Mus_orig_affine.EP5 = import_mus('EP5_orig_affine')
Mus_comp_affine.EP5 = import_mus('EP5_comp_affine')

Mus_orig.EP4b = import_mus('EP4b_orig')
Mus_comp.EP4b = import_mus('EP4b_comp')
#Mus_8peak.EP4b = import_mus('EP4b_8peak')
#Mus_orig_affine.EP4b = import_mus('EP4b_orig_affine')
Mus_comp_affine.EP4b = import_mus('EP4b_comp_affine')

##### scale the remaining 4 channels (scatter channels) "FSC-A", "FSC-H", "FSC-W", "SSC-A" to the (0,1) range by dividing by the maximum:
library(plyr)
# save the maximum for each channel
save.max = function(x) { # x: dataset
	sapply(x[, c("FSC-A", "FSC-H", "FSC-W", "SSC-A")], max)
}
scale.sc = function(x, maxima.from = x) { # x: dataset
	maxima = save.max(maxima.from); print(maxima)
	mutate(x, "FSC-A"=x$"FSC-A"/maxima["FSC-A"], "FSC-H"=x$"FSC-H"/maxima["FSC-H"], "FSC-W"=x$"FSC-W"/maxima["FSC-W"], "SSC-A"=x$"SSC-A"/maxima["SSC-A"]) 	
} 
Mus_orig_sc.EP4b = scale.sc(Mus_orig.EP4b)
Mus_comp_sc.EP4b = scale.sc(Mus_comp.EP4b)
Mus_comp_affine_sc.EP4b = scale.sc(Mus_comp_affine.EP4b)

Mus_orig_sc.EP5 = scale.sc(Mus_orig.EP5, Mus_orig.EP4b)
Mus_comp_sc.EP5 = scale.sc(Mus_comp.EP5, Mus_comp.EP4b)
Mus_comp_affine_sc.EP5 = scale.sc(Mus_comp_affine.EP5, Mus_comp_affine.EP4b)

###############################

unique_clusters = function(cell.subset, data = Mus_norm.st, plot.pairs = F) {
	CS = cell.subset
	df = data
	df = rbind(df[df$Subset != CS,], df[df$Subset == CS,])
	df = df[!duplicated(df[,-c(1:2)], fromLast = T), ] # 3584 rows
	if (plot.pairs) pairs(df[,3:10], pch = 21, cex = .1, col=c("grey", "red")[(df$Subset==CS) + 1])
	df
}

############## find best training model using cross-validation:  
bin_class_train = function(cell.subset="Singlets", train.set="Mus_orig.EP4", test.set="Mus_orig.EP4", method="knn", metric = "Kappa", path.to = "classification/", trControl = trainControl(method = "repeatedcv", repeats = 1), ...) {
	# create a file name
	cat(paste("\tMethod = ", method, "\n", sep = ""))
	data.name = gsub("Mus_", "", strsplit(train.set, "\\.")[[1]][1])
	wfile.name = paste("Data_", data.name, "-Model_", method, "-CS_", cell.subset, ".txt", sep = "")
	# modify data sets according to the given cell.subsets
	train.set = eval(as.name(train.set))
	test.set = eval(as.name(test.set))
	train.set = unique_clusters(cell.subset = cell.subset, data = train.set)
	test.set = unique_clusters(cell.subset = cell.subset, data = test.set)
	# train model using 10-fold cross-validation
	if (method=="glmnet") { # if glmnet, use default lambda values from glmnet
		alpha = seq(0.1, 1, length = 10)
		gr = do.call(rbind, lapply(alpha, function(x) data.frame(alpha = x, lambda = glmnet(as.matrix(train.set[,-c(1:2)]), train.set$Subset==cell.subset, alpha = x)$lambda)))
		model = train(train.set[,-c(1:2)], factor(train.set$Subset==cell.subset, levels = c(FALSE, TRUE)), method = method, metric = "Kappa", trControl = trControl, tuneGrid = gr)
	}
	else {
		model = train(train.set[,-c(1:2)], factor(train.set$Subset==cell.subset, levels = c(FALSE, TRUE)), method = method, metric = "Kappa", trControl = trControl, ...)
	}
	# write tuned parameters on a file 
	param = model$bestTune
	write.table(param, paste(path.to, "tuning-", wfile.name, sep = ""), row.names = F, quote = F)
	pred = predict(model, test.set[,-c(1:2)])
	label = test.set$Subset==cell.subset
	write.table(data.frame(pred = pred, label = label), paste(path.to, "prediction-", wfile.name, sep = ""), row.names = F, quote = F)
}
# For information about tuning parameters: 
# modelLookup(model = NULL)
# getModelInfo(model = NULL, regex = TRUE, ...)
getModelInfo("knn")$knn$grid
getModelInfo("rf")$rf$grid
getModelInfo("glmnet")$glmnet$grid
getModelInfo("svmRadialWeights")$svmRadialWeights$grid
getModelInfo("nb")$nb$grid

######################
# 5 data sets (orig, 8peak, comp, orig_affine, comp_affine)
# 23 cell subsets 
# 5 classifiers
# for each combination, save best training model (=parameters tuned by cross-validation with kappa as performance measure) and binary predictions for the test data using this best training model
# then for each combination, compute performance measures (and SE?) e.g., sens, spec, ppv, npv, kappa
######################
#### perform 10-fold cross-validation for tuning then perform predictions
#### (with scatter channels scaled) 
CS.vec = levels(Mus_orig_sc.EP4b$Subset)[-grep("(Time vs FSC)|(DP)|(DN)", levels(Mus_orig_sc.EP4b$Subset))]
# [1] "CD4+"           "CD4+IFNg+"      "CD4+IL2+"      
# [4] "CD4+TNFa+"      "CD8+"           "CD8+IFNg+"     
# [7] "CD8+IL2+"       "CD8+TNFa+"      "Singlets"      
# [10] "SSC vs aAmine-" "SSC vs CD3+"    "CD4+CD107+"    
# [13] "CD8+CD107+"
set.seed(100); lapply(CS.vec, function(CS) {
	cat(paste("Cell Subset = ", CS, "\n", sep = ""))
try(	Map(function(tr, te) {
bin_class_train(cell.subset = CS, train.set = tr, test.set = te, method = "nb", path.to = "classification_scaled/") #caret default: fL(Laplace correction) = 0; usekernel(distribution type) = FALSE/TRUE
bin_class_train(cell.subset = CS, train.set = tr, test.set = te, method = "knn", tuneLength = 10, path.to = "classification_scaled/") #caret default: k=5,7,9,...
bin_class_train(cell.subset = CS, train.set = tr, test.set = te, method = "glmnet", path.to = "classification_scaled/") #caret default: alpha=0.1+(k-1)*0.9/(n-1) (range 0.1~1); lambda = 0.1+(k-1)*2.9/(3*n-1) (range 0.1~3)
bin_class_train(cell.subset = CS, train.set = tr, test.set = te, method = "svmRadialWeights", tuneLength = 5, path.to = "classification_scaled/") #caret default: sigma = fixed; C = 0.25, 0.5, 1, 2, ...; Weight = 1,2,3,... 
bin_class_train(cell.subset = CS, train.set = tr, test.set = te, method = "rf", tuneLength = 7, path.to = "classification_scaled/") #caret default: mtry=2,3,4,...
}, list("Mus_orig_sc.EP4b", "Mus_comp_sc.EP4b", "Mus_comp_affine_sc.EP4b"), list("Mus_orig_sc.EP5", "Mus_comp_sc.EP5", "Mus_comp_affine_sc.EP5"))
)
}) 
# a lot of warnings: 
# 'Numerical 0 probability for all classes with observation xxx'

########### summarize results
### contingency table (as a vector) 
fnames = list.files('classification_scaled', pattern = 'prediction-Data')
fnames = fnames[!grepl("(DN)|(DP)", fnames)] # No DN/DP cell subsets in the test data EP5
class_results = lapply(fnames, function(x) {
	pred.df = read.table(paste('classification_scaled/', x, sep = ""), head = T)
	# frequencies in the confusion matrix (pred.df): TN, FP, FN, TP
	freq = as.vector(xtabs(~factor(pred, levels = c(FALSE, TRUE))+factor(label, levels = c(FALSE, TRUE)), pred.df))
	spl = strsplit(x, '-|\\.')[[1]]
	dat = gsub("Data_", "", spl[2])
	method = gsub("Model_", "", spl[3])
	CS = gsub("CS_", "", spl[4])
	data.frame(data = dat, method = method, cell.subset = CS, TN = freq[1], FP = freq[2], FN = freq[3], TP = freq[4])
})
class_results.table = Reduce(rbind, class_results)
# append performance measures using 'mutate' function
class_results.table = mutate(class_results.table, sens = TP/(FN+TP), spec = TN/(TN+FP), ppv = TP/(TP+FP), npv = TN/(TN+FN)) #from library(plyr)
class_results.table[order(class_results.table$cell.subset),]
write.csv(class_results.table, "table_results_scaled.csv", row.names = F)

### plot sensitivity against 1-specificity for each combination of method and data type
pdf("plot_ROC_scaled.pdf")
lapply(unique(class_results.table$cell.subset), function(CS) {
	x=class_results.table[class_results.table$cell.subset==CS,]
	plot(0, xlim = c(0,1), ylim = c(0,1), xlab = "1-Specificity", ylab = "Sensitivity", main = CS, type = "n")
	points(1-x$spec, x$sens, col = c(2:6)[as.numeric(x$method)], pch = c(1:3)[as.numeric(x$data)])
	legend.text = c("Method:", levels(x$method), "Data:", levels(x$data), NA, NA)
	legend.col = c(0, 2:6, 0, rep(1,3), 0, 0)
	legend.pch = c(0, rep(15,5), 0, 1:3, 0, 0)
	legend("topright", legend=legend.text,pch=legend.pch,col=legend.col,ncol=2)
})
dev.off()

######## call tuned parameters, retrieve the model, and perform predictions
######## NOTE: the following needs to be tested: 
### call best training model and perform prediction on test data
bin_class_test = function(cell.subset="Singlets", train.set="Mus_orig_sc.EP4b", test.set="Mus_orig_sc.EP5", method="knn", path.to = "classification_scaled/", trControl = trainControl(method = "none"), ...) {
	# create a file name
	data.name = gsub("Mus_", "", strsplit(train.set, "\\.")[[1]][1])
	wfile.name = paste("Data_", data.name, "-Model_", method, "-CS_", cell.subset, ".txt", sep = "")
	train.set = eval(as.name(train.set))
	test.set = eval(as.name(test.set))
	train.set = unique_clusters(cell.subset = cell.subset, data = train.set)
	test.set = unique_clusters(cell.subset = cell.subset, data = test.set)
	# read optimized parameter values 
	param = read.table(paste('classification_scaled/tuning', '-Data_', data.name, '-Model_', method, '-CS_', cell.subset, '.txt', sep = ""), head = T)
	# rebuild training model with the optimized parameters
	model = train(train.set[,-c(1:2)], factor(train.set$Subset==cell.subset, levels = c(FALSE, TRUE)), method = method, trControl = trainControl(method = "none"), tuneGrid = param)
	# perform prediction
	pred = predict(model, newdata = test.set[,-c(1:2)])
	# true class
	label = test.set$Subset==cell.subset
	write.table(data.frame(pred = pred, label = label), paste(path.to, "prediction-", wfile.name, sep = ""), row.names = F, quote = F)
	# performance measures
}
CS.vec.test = setdiff(levels(Mus_orig_sc.EP5$Subset), c("Time vs FSC", "CD8+IL2+"))
lapply(CS.vec.test, function(CS) {
	cat(paste("Cell Subset = ", CS, "\n", sep = ""))
try(	Map(function(tr, te) {
bin_class_test(cell.subset = CS, train.set = tr, test.set = te, method = "nb") 
bin_class_test(cell.subset = CS, train.set = tr, test.set = te, method = "knn") 
bin_class_test(cell.subset = CS, train.set = tr, test.set = te, method = "glmnet") 
bin_class_test(cell.subset = CS, train.set = tr, test.set = te, method = "svmRadialWeights") 
bin_class_test(cell.subset = CS, train.set = tr, test.set = te, method = "rf") 
}, list("Mus_orig_sc.EP4b", "Mus_comp_sc.EP4b", "Mus_comp_affine_sc.EP4b"), list("Mus_orig.EP5", "Mus_comp.EP5", "Mus_comp_affine.EP5"))
)
}) 

