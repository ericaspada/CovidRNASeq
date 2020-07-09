# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152075
library(tidyverse)
library(ggplot2)

### Import data ###
filepath_1 = '/Users/ericaspada/Desktop/projects/COV/GSE152075_series_matrix.txt'
filepath_2 = '/Users/ericaspada/Desktop/projects/COV/GSE152075_raw_counts_GEO.txt'
metadata = read_tsv(filepath_1)
rawCounts = read_delim(filepath_2, " ")
head(metadata)
str(metadata)
head(rawCounts)
str(rawCounts)

### Clean and merge the raw data ### 
rawCounts_t = as.data.frame(t(rawCounts)) # transpose rawCounts to prepare for merging
rawCounts_t[1:25,1:25] # examine
gene_names = pull(rawCounts,Gene)
colnames(rawCounts_t) = gene_names
rawCounts_t = rawCounts_t[-1,]
rawCounts_t[1:25,1:25] # examine

print(metadata[1:25,1:3])
unique(metadata[])
metadata_keep = metadata[6:13,] # only rows 6-13 should be kept and merged into the rawCounts file
metadata_keep = as.data.frame(t(metadata_keep)) # transpose to prepare for merging
# TO DO: add viral load (n1_ct)

metadata_keep = rename(metadata_keep,age=V6)
metadata_keep = rename(metadata_keep,gender=V7)
metadata_keep = rename(metadata_keep,sequencingBatch=V8)
metadata_keep = rename(metadata_keep,covidPositive=V4) # rename columns

# unique(metadata_keep[,'sequencingBatch'])
# unique(metadata_keep[,'covidPositive'])
# unique(metadata_keep[,'age'])

# the age, gender, and positive values need to be parsed from text strings
metadata_final = metadata_keep %>%
  mutate(age=parse_number(as.character(age), na='age: Unknown')) %>% # extract age values
  mutate(gender=recode(gender, 'gender: F'='F', 'gender: M'='M', 'gender: not collected'='NC')) %>%
  mutate(covidPositive=recode(covidPositive,'sars-cov-2 positivity: pos'='Pos',
                              'sars-cov-2 positivity: neg'='Neg')) %>%
  select(age,gender,covidPositive) # keep only age, gender, covidPositive columns

metadata_final[1:25,1:3] # examine
unique(metadata_final$age)

# add a column for age group
hist(metadata_final$age, breaks=10) # look at how age is distributed
# break into groups of 30 and younger, 31 - 60, 61 and older

metadata_final = metadata_final %>%
  mutate(ageGroup = case_when(is.na(age) ~ 'Unknown',
                              age < 30 ~ '<30',
                              age >=30 & age <= 60 ~ '30-60',
                              TRUE ~ '>60'))

metadata_final[1:25,1:4] # examine
metadata_final$ageGroup # examine

metadata_final = metadata_final[-1,] # remove the first row
sample = colnames(metadata)[-1] # get the row names (sample labels)
sample
rownames(metadata_final) = sample # assign the row names

# merge the rawCounts_t and metadata_final dfs using the rownames
data = merge(metadata_final,rawCounts_t, by='row.names')
data[1:50,1:50] # examine
# check some samples manually to make sure merge worked correctly
str(data)
length(data)

# change count columns to numeric (first need to change to character)
data[6:length(data)] = sapply(data[6:length(data)], as.character)
data[6:length(data)] = sapply(data[6:length(data)], as.numeric)
dim(data)

# plot the raw counts for some samples using a histogram
sample1 = as.numeric(data[3,6:length(data)]) # a random observation
hist(sample1, breaks=300) # difficult to see
sample1_transformed = log2(sample1+1) # use log2(k+1) transformation to help visualize
hist(sample1_transformed, breaks=50)
sample2 = as.numeric(data[54,6:length(data)])
sample2_transformed = log2(sample2+1)
hist(sample2_transformed, breaks=50)

# remove genes that have >95% zero values
# 95% is arbitrary
pct_zeros = as.data.frame(sapply(data, function(x) {round(length(which(x==0))/length(x)*100,digits=2)}))
pct_zeros = cbind(row.names(pct_zeros), pct_zeros)
colnames(pct_zeros) = c('gene', 'pct')
pct_zeros = as_tibble(pct_zeros)

to_remove = pct_zeros %>%
  filter(pct>=95.00) %>%
  select(gene)
to_remove = pull(to_remove,gene) #pull creates a vector

to_remove #examine

data_cleaned = data %>%
  select(-all_of(to_remove))

# plot some samples from the cleaned dataset
sample1 = as.numeric(data_cleaned[22,6:length(data_cleaned)]) # a random observation
hist(sample1, breaks=300) # difficult to see
sample1_transformed = log2(sample1+1) # use log2(k+1) transformation to help visualize
hist(sample1_transformed, breaks=50)
sample2 = as.numeric(data_cleaned[11,6:length(data_cleaned)])
sample2_transformed = log2(sample2+1)
hist(sample2_transformed, breaks=50)

genes = colnames(data_cleaned[6:length(data_cleaned)])

### Normalization of raw read counts ###
# Step 1: Need to normalize the raw read counts to account for library size (also known as read depth: total number of reads per sample)
# This is to account for between sample effects
# Implement the median of ratios method which the DESeq2 package uses

data_norm_1 = data_cleaned[,6:19976] + 1 # add one to account for 0s
data_norm_1[1:25,1:25]
data_norm_1 = log(data_norm_1) # take the natural log
averages = sapply(data_norm_1, mean) # take the average of each column (gene) (this is the geometric average)
averages[1:50]
data_norm_2 = sweep(data_norm_1, 2, averages) # subtract the averages from the log counts
# sweep operates on a matrix either column or row-wise (subtraction is the default)
data_norm_2[16,4]
data_norm_1[16,4] - averages[4] # manual check
medians = apply(data_norm_2, 1, median) # calculate the median for each sample (each row) - this should return of vector of length 484
medians_e = exp(medians) # exponentiate the medians to convert back to reads - these are the scaling factors for each sample
data_norm = sweep(data_cleaned[,6:19976], 1, medians_e, '/') # divide the original read counts by the scaling factors

# look at the distribution of some samples
sample1 = as.numeric(data_norm[200,6:length(data_norm)]) # a random observation
hist(sample1, breaks=100) # difficult to see
sample1_transformed = log2(sample1+1) # use log2(k+1) transformation to help visualize
hist(sample1_transformed, breaks=50)
sample2 = as.numeric(data_norm[41,6:length(data_norm)])
hist(sample2, breaks=100)
sample2_transformed = log2(sample2+1)
hist(sample2_transformed, breaks=50)
# these appear more like negative binomial

### Transformation of the normalized read counts to prepare for PCA ###
data_norm_log2 = log2(data_norm + 1)
pca = prcomp(data_norm_log2, center=TRUE, scale.=TRUE) # run PCA with scale=true and center=true
summary(pca)

### Create bi-plots ###
biplot_covidPositive = ggplot(data=pca$x)+
  aes(x=pca$x[,1], y=pca$x[,2], col=data$covidPositive)+
  geom_point()  # plot the first two PCs and color on covidPositive
biplot_covidPositive

biplot_ageGroup = ggplot(data=pca$x)+
  aes(x=pca$x[,1], y=pca$x[,2], col=data$ageGroup)+
  geom_point()  # plot the first two PCs and color on ageGroup
biplot_ageGroup

biplot_gender = ggplot(data=pca$x)+
  aes(x=pca$x[,1], y=pca$x[,2], col=data$gender)+
  geom_point()  # plot the first two PCs and color on gender
biplot_gender



  
  




