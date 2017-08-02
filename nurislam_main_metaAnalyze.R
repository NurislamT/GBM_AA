source('metaAnalyze.R')
source('metaAnalyzePamr.R')

metaAnalysisName = 'nurislam_main' # name of file of processed expression data and prefix for output files of meta-analysis
studyMetadataFileName = 'nurislam_main_study_metadata.csv' # table of study metadata
sampleMetadataFileName = 'nurislam_sample_metadata.csv' # table of sample metadata
parentFolderPath = 'data' # path to folder that contains the expression data
denovo = FALSE # process the raw data de novo or load processed data
className = 'class' # column name in table of sample metadata that contains data to predict
classesTrain = c('AA', 'GBM') # values for multinomial classification
familyName = 'multinomial' # family option for glmnet
intercept = TRUE # intercept option for glmnet

# load the table of study metadata
studyMetadata = read.csv(studyMetadataFileName, header=TRUE, stringsAsFactors=FALSE)
rownames(studyMetadata) = studyMetadata[,'study']
studyMetadata[,'discovery'] = studyMetadata[,'discovery']==1
studyMetadata[,'validation'] = studyMetadata[,'validation']==1

# load the table of sample metadata
sampleMetadataTmp = read.csv(sampleMetadataFileName, header=TRUE, stringsAsFactors=FALSE)
rownames(sampleMetadataTmp) = sampleMetadataTmp[,'sample']
sampleMetadata = sampleMetadataTmp[sampleMetadataTmp[,'study'] %in% studyMetadata[,'study'],]

# load the expression data for all studies
if (denovo) {
	esetList = getStudyDataList(parentFolderPath, studyMetadata)
	saveRDS(esetList, file=paste0(metaAnalysisName, '.rds'))
} else {
	esetList = readRDS(paste0(metaAnalysisName, '.rds'))}

##############################################
###CHANGES OF ORIGINAL CODE#######
esetOrig=esetList[1][[1]]
featureDf = pData(featureData(esetOrig))
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character)
mapping = featureDf[,c('ID', 'Gene Symbol')]
mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
colnames(mapping) = c('probeSet', 'geneId')
exprsByGene = calcExprsByGene(esetOrig, mapping)
if (any(is.na(exprsByGene))) {
	warning(sprintf('Imputing missing expression values for study %s.', studyName))
	resultImputed = impute.knn(exprsByGene)
	exprsByGene = resultImputed$data}
colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
eset1 = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))

esetOrig=esetList[2][[1]]
featureDf = pData(featureData(esetOrig))
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character)
mapping = featureDf[,c('ID', 'Gene Symbol')]
mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
colnames(mapping) = c('probeSet', 'geneId')
exprsByGene = calcExprsByGene(esetOrig, mapping)
if (any(is.na(exprsByGene))) {
	warning(sprintf('Imputing missing expression values for study %s.', studyName))
	resultImputed = impute.knn(exprsByGene)
	exprsByGene = resultImputed$data}
colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
eset2 = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
	
esetOrig=esetList[3][[1]]
featureDf = pData(featureData(esetOrig))
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character)
mapping = featureDf[,c('ID', 'Gene Symbol')]
mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
colnames(mapping) = c('probeSet', 'geneId')
exprsByGene = calcExprsByGene(esetOrig, mapping)
if (any(is.na(exprsByGene))) {
	warning(sprintf('Imputing missing expression values for study %s.', studyName))
	resultImputed = impute.knn(exprsByGene)
	exprsByGene = resultImputed$data}
colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
eset3 = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
	
esetOrig=esetList[4][[1]]
featureDf = pData(featureData(esetOrig))
idx = sapply(featureDf, is.factor)
featureDf[idx] = lapply(featureDf[idx], as.character)
mapping = featureDf[,c('ID', 'Gene Symbol')]
mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
colnames(mapping) = c('probeSet', 'geneId')
exprsByGene = calcExprsByGene(esetOrig, mapping)
if (any(is.na(exprsByGene))) {
	warning(sprintf('Imputing missing expression values for study %s.', studyName))
	resultImputed = impute.knn(exprsByGene)
	exprsByGene = resultImputed$data}
colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
eset4 = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
	
esetList=c(eset1,eset2,eset3,eset4)
names(esetList)=c("GSE1993","GSE4271","GSE4290","GSE4412")
##############################END OF CHANGES###################

# select samples that are in sampleMetadata, extract the expression matrix
ematList = cleanStudyData(esetList, sampleMetadata)

# merge expression data and perform cross-study normalization of discovery datasets
ematMergedDiscoveryAllClasses = mergeStudyData(ematList[studyMetadata[studyMetadata[,'discovery'], 'study']],
															  sampleMetadata, covariateName=NA)

# select samples for discovery
discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
idxDiscovery = (sampleMetadata[,'study'] %in% discoveryStudyNames) & (sampleMetadata[,className] %in% classesTrain)

discoverySampleNames = rownames(sampleMetadata)[idxDiscovery]

# make table of sample metadata and matrix of expression data for discovery samples
sampleMetadataDiscovery = sampleMetadata[idxDiscovery,]
sampleMetadataDiscovery[,className] = factor(sampleMetadataDiscovery[,className], levels=classesTrain)
ematMergedDiscovery = ematMergedDiscoveryAllClasses[,sampleMetadataDiscovery[,'sample']]

# construct foldid and weights for glmnet
glmnetArgs = makeGlmnetArgs(sampleMetadataDiscovery)

# run leave-one-study-out cross-validation on discovery data
alphas = 0.9
cvFitList = crossValidateMerged(ematMergedDiscovery, sampleMetadata, yName=className, weights=glmnetArgs$weights,
										  foldid=glmnetArgs$foldid,alphas=alphas, family=familyName, lambda.min.ratio=0.001,
										  intercept=intercept, keep=TRUE)

# extract the glmnet fit
alpha = 0.9
cvFit = cvFitList[[alphas==alpha]]
fitResult = cvFit$glmnet.fit
lambda = cvFit$lambda.min

# plot the deviance vs. lambda
plotCvError(cvFit, metaAnalysisName=metaAnalysisName, size=0.4, ggplotArgs=list(theme_bw()))

# plot the coefficients vs. lambda
pdf(file=sprintf('%s_lambda_coef.pdf', metaAnalysisName), width=6, height=4.5)
plot(fitResult, xvar='lambda', xlim=c(-6, 0))
dev.off()

# write tables of coefficients
writeCoefficients(fitResult, lambda=lambda, metaAnalysisName=metaAnalysisName)

# set order of gene ids for plots
coefDf = makeCoefDf(coef(fitResult, s=lambda), decreasing=FALSE, classLevels=classesTrain)
coefDf = coefDf[coefDf[,'geneId']!='(Intercept)',]
geneIdOrder = coefDf[,'geneId']

# plot the coefficients for particular values of lambda
plotCoefficients(fitResult, lambda, classLevels=classesTrain, geneIdOrder=geneIdOrder,
					  metaAnalysisName=metaAnalysisName, width=7, height=10,
					  ggplotArgs=list(scale_fill_brewer(type='qual', palette=3), scale_y_continuous(breaks=c(-0.3, 0, 0.3)), theme_bw()))

# plot the expression of selected genes in the discovery set
annoNames = c('study', 'class')
annoLevels = list(studyMetadata[studyMetadata[,'discovery'],'study'], classesTrain)
names(annoLevels) = annoNames
annoColors = list(brewer.pal(4, 'Paired'))
names(annoColors) = 'class'
names(annoColors[[1]]) = classesTrain

plotExpressionHeatmapMerged(fitResult, lambda, ematMergedDiscovery, sampleMetadata, annoNames, annoLevels,
									 annoColors, clusterSamplesTogether=FALSE, geneIdOrder=geneIdOrder, className=className,
									 classLevels=classesTrain, metaAnalysisName=metaAnalysisName, width=10, height=7, fontsize_row=7)

# write confusion matrix for cross-validation
writeConfusionCrossValidation(cvFit, lambda=lambda, ematMerged=ematMergedDiscovery,
										sampleMetadata=sampleMetadata, className=className, classLevels=classesTrain,
										metaAnalysisName=metaAnalysisName)

# plot the class probabilities for cross-validation
plotClassProbsCrossValidation(cvFit, lambda, sampleMetadata, discoveryStudyNames, discoverySampleNames,
										className, classesTrain, metaAnalysisName, size=2, width=8, height=15,
										ggplotArgs=list(scale_color_brewer(type='qual', palette=3), theme_bw(), theme(legend.title=element_blank())))

# predict class for the validation datasets
predsList = predictValidationData(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
											 alpha=alpha, lambda=lambda, weights=glmnetArgs$weights, covariateName=NA,
											 className=className, familyName=familyName, intercept=intercept)

# write confusion matrix for validation
writeConfusionValidation(predsList, lambda=lambda, sampleMetadata=sampleMetadata,
								 className=className, classLevels=classesTrain, metaAnalysisName=metaAnalysisName)
writeConfusionValidationEach(predsList, lambda=lambda, sampleMetadata=sampleMetadata,
									  className=className, classLevels=classesTrain, metaAnalysisName=metaAnalysisName)

# plot class probabilities for validation
plotClassProbsValidation(predsList, sampleMetadata, className, classesTrain, metaAnalysisName, size=2, width=8, height=9,
								 ggplotArgs=list(scale_color_brewer(type='qual', palette=3), theme_bw(), theme(legend.title=element_blank())))


#################################################### run PAM for comparison
# cross-validation on discovery datasets
cvFitPam = crossValidateMergedPam(ematMerged=ematMergedDiscovery, sampleMetadata=sampleMetadata, foldid=glmnetArgs$foldid)

# plot the results of cross-validation
pdf(file=sprintf('%s_pamr_cv.pdf', metaAnalysisName), width=6, height=9)
pamr.plotcv(cvFitPam)
dev.off()

# write confusion matrix for cross-validation
writeConfusionCrossValidationPam(cvFitPam, 3, sampleMetadata, discoverySampleNames, classesTrain, className, metaAnalysisName)
writeConfusionCrossValidationPam(cvFitPam, 15, sampleMetadata, discoverySampleNames, classesTrain, className, metaAnalysisName)


#list of genes
foldid=glmnetArgs$foldid
ematMerged=ematMergedDiscovery
y = sampleMetadata[colnames(ematMerged), className]
x = t(scale(t(ematMerged), center=TRUE, scale=FALSE))
foldsTmp = foldid[colnames(ematMerged)]
folds = split(seq_along(foldsTmp), foldsTmp)
pamData = list(x=x, y=y, geneid=rownames(ematMerged))
pamTrain = pamr.train(pamData)
genes_list = pamr.listgenes(pamTrain,pamData, threshold=5.32939)[,'id']
genes=genes_list
genes[122]="SNRPN"
genes[76]= "TXNDC5"
genes[190]="MIR21"
genes[206]="FAM45A"
genes[57]="TMSB4X"
genes[68]="MYL12A"
geneTexts = sprintf('%s', genes_list)
names(geneTexts)=geneTexts

emat = ematMerged[genes_list,]
d = dist(t(emat))
co = order.optimal(d, hclust(d)$merge)
emat = emat[,co$order]

d = dist(emat)
co = order.optimal(d, hclust(d)$merge)
emat = emat[co$order,]
rownames(emat) = geneTexts[co$order]

annotation = sampleMetadata[colnames(ematMerged), annoNames, drop=FALSE]
for (annoName in annoNames) {
		if (!is.na(annoLevels[[annoName]][1])) {
		annotation[,annoName] = factor(annotation[,annoName], levels=annoLevels[[annoName]])}}
	
# heatmap

emat=ematMergedDiscovery
pdf(file=sprintf('hello.pdf'), width=width, height=height)
pheatmap(emat, color=colorRampPalette(rev(greenred(32)))(100),
breaks=seq(from=-3, to=3, length.out=101), cluster_rows=TRUE, cluster_cols=TRUE, treeheight_row=0,fontsize_row=2.75,
treeheight_col=0, show_colnames=FALSE, border_color=NA, annotation=annotation, annotation_colors=annoColors)
	dev.off()


	
# predict class for the validation datasets
predsListPam = predictValidationDataPam(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
													 threshold=cvFitPam$threshold[15], covariateName=NA, className=className)

# write confusion matrix for validation
writeConfusionValidationPam(predsListPam, threshold=cvFitPam$threshold[15], sampleMetadata=sampleMetadata,
									 className=className, classLevels=classesTrain, metaAnalysisName=metaAnalysisName)
writeConfusionValidationEachPam(predsListPam, threshold=cvFitPam$threshold[15], sampleMetadata=sampleMetadata,
										  className=className, classLevels=classesTrain, metaAnalysisName=metaAnalysisName)
