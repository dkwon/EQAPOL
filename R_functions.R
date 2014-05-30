#library(RWeka) # kNN
#library('e1071') # naive Bayes and SVM 
#library('ada') # adaboost
library(glmnet) # penalized logistic regression
#library(randomForest) # random forest
#library(nnet) # neural network
library(caret) # confusionMatrix # for training/testing models, tuning model parameters
	#Model list in caret:
	#http://caret.r-forge.r-project.org/modelList.html

# order of channels for the peaks data: 
peaks.channels = unname(unlist(read.table("./peaks_ep4/peak_channels.txt", stringsAsFactors = F)))
# "FITC" "PE" "PerCP-Cy5.5" "PE-Cy7" "APC-H7" "APC" "BV421" "AquaAmine" 
# order of channels for the original data: 
orig.channels = unname(sapply(strsplit(unlist(read.table("./EP5_orig/001_parameters.txt", stringsAsFactors = F)), "_"), function(x) ifelse(length(x)<=2, paste(head(x,2), collapse = "-"), tail(x,1))))
# "PE" "APC-H7" "PE-Cy7" "PerCP-Cy5.5" "APC" "BV421" "FITC" "FSC-A" "FSC-H" "FSC-W" "SSC-A" "AquaAmine"

# exclude paramters.txt or parameter_map.txt and assign channel names (only for raw data)
assign.names = function(fname, path = "EP4_orig/") {
	# get the original data
	dat = read.table(paste(path, fname, sep = ""))
	# read off site id from the file name
	site = strsplit(fname, "_")[[1]][1]
	# get the vector of (simplified) channel names used for this site
	channels = strsplit(unlist(read.table(paste(path, site, "_parameters.txt", sep=""), stringsAsFactors = F)), "_") 
	channels = sapply(channels, function(x) ifelse(length(x)<=2, paste(head(x,2), collapse = "-"), tail(x,1)))
	# assign the original data channel [column] names
	names(dat) = channels
	list(site = site, data = dat)
}

# linear interpolation using 8 peak medians:
lin.int = function(x, ref, peak.ll = -45000, peak.ul = 110000) { # x: input value, ref: 8 peak medians (a vector of length 8)
	stopifnot(length(ref)==8)
	# multiply by 10^5 (only for the 8 non-scatter channels)
	x = x*1e+5
	if (x <= peak.ll) return(0) 
	else if (x >= peak.ul) return(9) # in the analysis, discard the entire cluster that has 0 or 9
	else {
		# reference vector {lower limit, peak1, peak2, ..., peak8, upper limit}
		ref = c(peak.ll, ref, peak.ul)
		#1. find the index of the smallest ref value to the right of the input value 
		i = min(which(ref > x))
		#2. linear interpolation
		li = (i-1) + (x-ref[i-1]) / (ref[i]-ref[i-1])
		stopifnot(li>=1, li<=10)
		#3. make li match up the peak number
		return(li - 1)
	}
}
# linear int-/extra- polation using low and high peak medians:
lin.int_ext = function(x, ref) { # x: input value, ref: compensation bead peak medians (a vector of length 2: low, high)
	stopifnot(length(ref)==2, ref[2]>ref[1])
	# multiply by 10^5 (only for the 8 non-scatter channels)
	x = x*1e+5
	# find a straight line that goes through the low and high peak medians
	# and perform linear interpolation or extrapolation; 
	# output value will be within [0,1] if interpolation and <0 or >1 if extrapolation
	(x-ref[1])/(ref[2]-ref[1])
}

# normalizing with 8-peak beads (linear interpolation):
norm.li = function(fname, path.from = "EP4_orig/", path.peaks = "peaks_ep4/", path.to = "EP4_8peak/") { 
	named.dat = assign.names(fname, path = path.from)
	site = named.dat$site
	dat = named.dat$data
	# get peak medians for this site: each file (=each lab) has a 8(channels)x8(peaks) matrix
	peaks = read.table(paste(path.peaks, "peak_medians_", site, ".txt", sep=""))
	rownames(peaks) = peaks.channels
	# perform linear interpolation on the original data using the peak medians
	li = lapply(peaks.channels, function(ch) sapply(dat[, ch], lin.int, ref = unlist(peaks[ch, ])))
	li = do.call(cbind, li)
	colnames(li) = peaks.channels
	rownames(li) = NULL
	# append untransformed columns (channels)
	li = cbind(li, dat[, setdiff(orig.channels, peaks.channels)])
	li = li[, orig.channels]
	# write a file
	write.table(li, file = paste(path.to, fname, sep = ""), quote = F, row.names = F)	
} # test: norm.li("036_CD4+_mus.txt")

# compensation bead peak medians data: visit, site, marker, low, high
comp_beads.EP4 = read.table("compbeads_ep4/comp_beads_low_high_peaks.txt", check.names = F, stringsAsFactors = F)
names(comp_beads.EP4) = c("visit", "site", "marker", "low", "high") #marker means channel
comp_beads.EP4$site = sprintf("%03d", comp_beads.EP4$site)
comp_beads.EP4$marker = sapply(strsplit(comp_beads.EP4$marker, split = "_"), '[', 2) 
comp_beads.EP5 = read.table("compbeads_ep5/comp_beads_low_high_peaks.txt", check.names = F, stringsAsFactors = F)
names(comp_beads.EP5) = c("visit", "site", "marker", "low", "high") #marker means channel
comp_beads.EP5$site = sprintf("%03d", comp_beads.EP5$site)
comp_beads.EP5$marker = sapply(strsplit(comp_beads.EP5$marker, split = "_"), '[', 2) 

# normalizing with compensation bead peaks (linear inter-/extra- polation):
norm.li.comp = function(fname, path.from = "EP4_orig/", data.peaks = comp_beads.EP4, path.to = "EP4_comp/") { 
	named.dat = assign.names(fname, path = path.from)
	site = named.dat$site
	dat = named.dat$data
	# get peak medians for this site: 8(channels)x2(peaks) matrix
	peaks = data.peaks[data.peaks$site == site, c("marker", "low", "high")]
	rownames(peaks) = peaks$marker # marker means channel
	peaks = peaks[, c("low", "high")]
	stopifnot(dim(peaks) == c(8,2))
	# perform linear inter-/extra- polation on the original data using the peak medians
	li = lapply(peaks.channels, function(ch) sapply(dat[, ch], lin.int_ext, ref = unlist(peaks[ch, ])))
	li = do.call(cbind, li)
	colnames(li) = peaks.channels
	rownames(li) = NULL
	# append untransformed columns (channels)
	li = cbind(li, dat[, setdiff(orig.channels, peaks.channels)])
	li = li[, orig.channels]
	# write a file
	write.table(li, file = paste(path.to, fname, sep = ""), quote = F, row.names = F)	
} # test: norm.li.comp("001_CD4+_mus.txt")

# affine parameters (in the e-mail from Cliburn):
aff_param = c('FCM_APC-H7_CD3_A', 'VIA_AquaAmine__A', 'FCM_PE-Cy7_CD4_A', 'FCM_PerCP-Cy5.5_CD8_A','FCM_APC_IFNg_A', 'FCM_BV421_IL2_A', 'FCM_FITC_TNFa_A', 'FCM_PE_CD107a_A')
aff_param = sapply(strsplit(aff_param, "_"), '[', 2)

fnames_affine.EP4 = list.files(path="affine_ep4", full.names = T)
fnames_affine.EP5 = list.files(path="affine_ep5", full.names = T)

affine = function(fname, data.type = c("orig", "comp"), path.from = "EP4_orig/", fnames_affine = fnames_affine.EP4, path.to = paste("EP4_", data.type, "_affine/", sep = "")) { 
	data.type = match.arg(data.type)
	if (data.type == "orig") {
		named.dat = assign.names(fname, path = path.from)
		site = named.dat$site
		dat = named.dat$data
	} else {
		site = strsplit(fname, "_")[[1]][1]
		dat = read.table(paste(path.from, fname, sep = ""), header = T, check.names = F)
	}
	# get parameters for affine transformation for this site: a = 8(channels)x8(channels) matrix, b = vector of length 8(channels)
	index.a=grep(paste("DKL_A_", site, "_[0-9]{4}_", data.type, ".txt", sep = ""), fnames_affine)
	index.b=grep(paste("DKL_B_", site, "_[0-9]{4}_", data.type, ".txt", sep = ""), fnames_affine)
	a=read.table(fnames_affine[index.a])
	b=unlist(read.table(fnames_affine[index.b]))
	stopifnot(dim(a)==c(8,8), length(b)==8)
	# perform affine transformation:
	# In Python, transformation was defined as 'xshift = numpy.dot(x,a)+b', no matter whether x is a vector (single cluster) of length 8 or an nx8 matrix (set of clusters):
	# Note: numpy.dot performs the usual matrix multiplication and x, if a vector, and b are row vecors. 
	# In R, 
	# if x is a vector of length 8 (single cluster), use 'colSums(x1*a)+b'  
	# if x is an nx8 matrix (set of clusters), use 'sweep(x[,aff_param]%*%a, 2, b, "+")' 
	dat.shift = sweep(as.matrix(dat[,aff_param])%*%as.matrix(a), 2, b, "+") 
	colnames(dat.shift) = aff_param
	# append untransformed columns (channels)
	dat.shift = cbind(dat.shift, dat[, setdiff(orig.channels, aff_param)]) # data.frame
	dat.shift = dat.shift[, orig.channels]
	# write a file
	write.table(dat.shift, file = paste(path.to, fname, sep = ""), quote = F, row.names = F)	
} # test: affine("004_CD4+TNFa+_mus.txt")

# import data
import_mus = function(path = "EP4_orig", is.rawdata = grepl("^EP((4b?)|5)_orig/?$", path)) {
	if (!grepl("/$", path)) path = paste(path, "/", sep="")
	cat(paste("Is raw data?\n", is.rawdata, "\n\n", sep = ""))
	# file names: Site_Subset_mus.txt 
	# each file corresponds to mus (no sigmas now) for each cell subset for each lab [site]; each file contains a nxk (mus) matrix where n=number of clusters, k=number of channels [colors], k=14 for raw data and 8 for processeced (beads-normalized or affine-transformed)
	# find all files 
	fnames = list.files(path)
	if (is.rawdata) {
		#exclude parameters, parameter maps, ungated (=no cell subset information) clusters
		fnames = fnames[!grepl("(parameter)|(^[[:digit:]]{3}_mus.txt)", fnames)]
	} 
	# data frame of file names for error checking
	#fnames.df = data.frame(t(sapply(strsplit(gsub(".txt", "", fnames), "_"), '['))) # data frame of file names
	fnames.df = data.frame(t(sapply(strsplit(gsub(" Amine", " aAmine", gsub(".txt", "", fnames)), "_"), '['))) 
	names(fnames.df) = c("Site", "Subset", "Stat")
	# ##### check errors/typos here: 
	# print(fnames.df)
	# print(sapply(fnames.df, summary))
	# print(sapply(fnames.df, function(x) length(unique(x))))
	# #####
	# read off mus (number of clusters x 8 or 14 channels) from each file
	# then group into a data frame (site, cell subset, mu_1, mu_2, ..., mu_k)
	# Note: in the normalized data, there is no empty files, so error checking part of the following code is doing is unnecessary 
	# Note: in the normalized data, there is column names (header) in each file
	if (is.rawdata) {
		try.read = lapply(fnames, function(x) try(data.frame(Site = strsplit(x, "_")[[1]][1],  Subset = gsub(" Amine", " aAmine", strsplit(x, "_")[[1]][2]), assign.names(x, path)$data[, orig.channels], check.names = F), silent = T))
	} else {
		try.read = lapply(fnames, function(x) try(data.frame(Site = strsplit(x, "_")[[1]][1],  Subset = gsub(" Amine", " aAmine", strsplit(x, "_")[[1]][2]), read.table(paste(path, x, sep = ""), header = T, check.names = F), check.names = F), silent = T))
	}
	empty = sapply(try.read, function(x) inherits(x, "try-error")) #identifies the elements of the list with the try-error i.e. empty files 
	Mus = do.call(rbind, try.read[!empty & grepl("mus", fnames)]) #removes the elements identified above
	### Optional: only for 8-peak beads normalization (interpolation)
	# Mus[,-(1:2)][Mus[,-(1:2)]<0]
	# Mus = Mus_norm[which(apply(Mus[,-(1:2)], 1, function(x) any(x<=0|x>=9))), ]
	print(lapply(Mus[,1:2], summary))
	#discard rows that have negative values for scatter channels
	neg.row = apply(Mus[, !orig.channels %in% peaks.channels], 1, function(x) any(x<0))
	Mus[!neg.row,]
}
