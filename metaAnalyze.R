require('GEOquery')
require('affy')
require('annotate')
require('org.Hs.eg.db')
require('foreach')
require('impute')
require('dplyr')
require('sva')
require('glmnet')
require('beeswarm')
require('reshape2')
require('ROCR')
require('cba')
require('pheatmap')
require('ggplot2')
require('gridExtra')
require('RColorBrewer')


fixCustomCdfGeneIds = function(geneIds) {
	return(sub('_at', '', geneIds))}


fixGeoSampleNames = function(sampleNames) {
	sampleNames = paste0(toupper(sampleNames), '_')
	regexResult = regexpr('^GSM[0-9]+[^0-9]', sampleNames)
	sampleNamesNew = mapply(function(sampleName, matchLength) substr(sampleName, 1, matchLength-1),
									sampleNames, attr(regexResult, 'match.length'))
	return(sampleNamesNew)}


fixCelSampleNames = function(sampleNames) {
	sampleNamesNew = gsub('\\.cel$', '', sampleNames, ignore.case=TRUE)
	return(sampleNamesNew)}


getGeneProbeMappingAffy = function(mappingFilePath) {
	mapping = read.table(mappingFilePath, sep='\t', header=TRUE, stringsAsFactors=FALSE)
	mappingUnique = unique(mapping[,c('Probe.Set.Name', 'Affy.Probe.Set.Name')])
	rownames(mappingUnique) = NULL
	colnames(mappingUnique) = c('geneId', 'probeSet')
	return(mappingUnique)}


getGeneProbeMappingDirect = function(featureDf, geneColname, probeColname='ID') {
	mapping = featureDf[,c(probeColname, geneColname)]
	mapping = mapping[apply(mapping, MARGIN=1, function(x) all(!is.na(x) & x!='')),]
	mapping = data.frame(lapply(mapping, as.character), stringsAsFactors=FALSE)
	colnames(mapping) = c('probeSet', 'geneId')
	return(mapping)}


getGeneProbeMappingAnno = function(featureDf, dbName, interName) {
	mappingProbeIntermediate = featureDf[!is.na(featureDf[,interName]) & featureDf[,interName]!='', c('ID', interName)]
	colnames(mappingProbeIntermediate) = c('probeSet', 'geneInter')
	mapTmp1 = eval(parse(text=dbName))
	mapTmp2 = mappedkeys(mapTmp1)
	mapTmp3 = as.list(mapTmp1[mapTmp2])
	geneId = do.call(c, mapTmp3)
	geneInter = do.call(c, mapply(function(inter, len) rep_len(inter, len), names(mapTmp3), sapply(mapTmp3, length),
											SIMPLIFY=FALSE))
	if (dbName=='org.Hs.egUNIGENE2EG') {
		geneInter = sub('Hs.', '', geneInter, fixed=TRUE)}
	mappingIdInter = data.frame(geneId, geneInter, stringsAsFactors=FALSE)
	mapping = merge(mappingIdInter, mappingProbeIntermediate, by='geneInter', sort=FALSE)}


calcExprsByGene = function(eset, mapping) {
	geneIds = unique(mapping[,'geneId'])
	exprsByGene = matrix(nrow=length(geneIds), ncol=ncol(eset), dimnames=list(geneIds, sampleNames(eset)))
	for (geneId in geneIds) {
		exprsTmp = exprs(eset)[mapping[mapping[,'geneId']==geneId, 'probeSet'],, drop=FALSE]
		if (nrow(exprsTmp)==1) {
			exprsByGene[geneId,] = exprsTmp
		} else {
			exprsByGene[geneId,] = rowMedians(t(exprsTmp), na.rm=TRUE)}}
	return(exprsByGene)}


getSupportedPlatforms = function() {
	return(c('GPL96','GPL570','GPL180', 'GPL885', 'GPL887', 'GPL962', 'GPL1053', 'GPL1291', 'GPL1293', 'GPL1390',
				'GPL1708', 'GPL5645', 'GPL6254', 'GPL6480', 'GPL6884', 'GPL6947', 'GPL7015'))}


getStudyData = function(parentFolderPath, studyName, studyDataType, platformInfo) {
	cat(sprintf('Loading study %s...\n', studyName))
	# load the data, convert to gene id, normalize and transform where necessary
	if (studyDataType %in% c('affy_geo', 'affy_custom')) {
		require(platformInfo, character.only=TRUE)
		cwd = setwd(file.path(parentFolderPath, studyName))
		eset = justRMA(cdfname=platformInfo)
		setwd(cwd)
		featureNames(eset) = fixCustomCdfGeneIds(featureNames(eset))
		if (studyDataType=='affy_geo') {
			sampleNames(eset) = fixGeoSampleNames(sampleNames(eset))
		} else {
			sampleNames(eset) = fixCelSampleNames(sampleNames(eset))}
		
	} else if (studyDataType=='affy_series_matrix') {	
		mapping = getGeneProbeMappingAffy(file.path(parentFolderPath, paste0(platformInfo, '_mapping.txt')))
		esetOrig = getGEO(filename=file.path(parentFolderPath, paste0(studyName, '_series_matrix.txt')))
		exprs(esetOrig)[exprs(esetOrig)<=0] = min(exprs(esetOrig)[exprs(esetOrig)>0])
		exprs(esetOrig) = log2(exprs(esetOrig))
		exprsByGene = calcExprsByGene(esetOrig, mapping)
		rownames(exprsByGene) = fixCustomCdfGeneIds(rownames(exprsByGene))
		colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
		eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
		
	} else if (studyDataType=='series_matrix') {
		supportedPlatforms = getSupportedPlatforms()
		
		if (!(platformInfo %in% supportedPlatforms)) {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}
		
		esetOrig = getGEO(filename=file.path(parentFolderPath, paste0(studyName,'_series_matrix.txt')))
		if (is.list(esetOrig) && length(esetOrig)==1) {
			esetOrig = esetOrig[[1]]}
		
		featureDf = pData(featureData(esetOrig))
		idx = sapply(featureDf, is.factor)
		featureDf[idx] = lapply(featureDf[idx], as.character)
		if (platformInfo=='GPL180') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE_SYM')
		} else if (platformInfo=='GPL885') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL887') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL962') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egUNIGENE2EG', interName='UNIGENE')
		} else if (platformInfo=='GPL1053') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='GENE')
		} else if (platformInfo=='GPL1291') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
		} else if (platformInfo=='GPL1293') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez Gene ID')
		} else if (platformInfo=='GPL1390') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_ACC')
		} else if (platformInfo=='GPL1708') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL5645') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egSYMBOL2EG', interName='Gene Name')
		} else if (platformInfo=='GPL6254') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egENSEMBL2EG', interName='ENSEMBL_GENE_ID')
		} else if (platformInfo=='GPL6480') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='GENE')
		} else if (platformInfo=='GPL6884') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
		} else if (platformInfo=='GPL6947') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
		} else if (platformInfo=='GPL7015') {
			mapping = getGeneProbeMappingAnno(featureDf, dbName='org.Hs.egREFSEQ2EG', interName='GB_LIST')
		} else if (platformInfo=='GPL96') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
		} else if (platformInfo=='GPL570') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='Entrez_Gene_ID')
		} else {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}
		
		exprsByGene = calcExprsByGene(esetOrig, mapping)
		if (any(is.na(exprsByGene))) {
			warning(sprintf('Imputing missing expression values for study %s.', studyName))
			resultImputed = impute.knn(exprsByGene)
			exprsByGene = resultImputed$data}
		colnames(exprsByGene) = fixGeoSampleNames(colnames(exprsByGene))
		eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
		
	} else if (studyDataType=='eset_rds') {
		esetOrig = readRDS(file.path(parentFolderPath, paste0(studyName, '.rds')))
		
		featureDf = pData(featureData(esetOrig))
		if (platformInfo=='ready') {
			return(esetOrig)
		} else if (platformInfo=='rosetta') {
			mapping = getGeneProbeMappingDirect(featureDf, geneColname='EntrezGene.ID', probeColname='probe')
		} else {
			warning(sprintf('Study %s not loaded, because platform %s is not currently supported.', studyName, platformInfo))
			return(NA)}
		
		exprsByGene = calcExprsByGene(esetOrig, mapping)
		if (any(is.na(exprsByGene))) {
			warning(sprintf('Imputing missing expression values for study %s.', studyName))
			resultImputed = impute.knn(exprsByGene)
			exprsByGene = resultImputed$data}
		eset = ExpressionSet(assayData=exprsByGene, phenoData=phenoData(esetOrig), experimentData=experimentData(esetOrig))
		
	} else {
		warning(sprintf('Study %s not loaded, because data type %s is not currently supported.', studyName, studyDataType))
		eset = NA}
	return(eset)}


getStudyDataList = function(parentFolderPath, studyMetadata) {
	esetList = foreach(studyName=rownames(studyMetadata)) %do% {
		if (any(is.na(studyMetadata[studyName,]))) {
			NA
		} else {
			getStudyData(parentFolderPath, studyName, studyMetadata[studyName, 'studyDataType'],
							 studyMetadata[studyName, 'platformInfo'])}}
	names(esetList) = rownames(studyMetadata)
	return(esetList[!is.na(esetList)])}


cleanStudyData = function(esetList, sampleMetadata) {
	# select relevant samples, convert to matrix
	ematList = foreach(studyName=names(esetList)) %do% {
		keepIdx = colnames(esetList[[studyName]]) %in% sampleMetadata[sampleMetadata[,'study']==studyName, 'sample']
		exprs(esetList[[studyName]])[,keepIdx]}
	names(ematList) = names(esetList)
	return(ematList)}


makeMatchSampleMapping = function(metadata, subStudyNames, matchColname) {
	metadataNow = metadata[metadata[,'study'] %in% subStudyNames, c('study', 'sample', matchColname)]
	metadataNow = metadataNow[order(metadataNow[,'study'], decreasing=is.unsorted(subStudyNames)), c('sample', matchColname)]
	headFunc = function(x) x[1]
	mappingDf = metadataNow %>% group_by_(.dots=list(matchColname)) %>% summarise_each(funs(headFunc)) %>%
		data.frame(check.names=FALSE)
	mapping = mappingDf[,'sample']
	names(mapping) = mappingDf[,matchColname]
	return(mapping)}


mergeMatchStudyData = function(ematAtomicList, studyMetadataAtomic, sampleMetadataAtomic, matchColname,
										 mergeFunc=function(x) mean(x, na.rm=TRUE)) {
	ematList = list()
	sampleMetadataList = list()
	
	for (matchStudyName in unique(studyMetadataAtomic[,'matchStudy'])) {
		if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)==1) {
			ematList[[matchStudyName]] = ematAtomicList[[matchStudyName]]
			sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[sampleMetadataAtomic[,'study']==matchStudyName,]
			
		} else if (sum(studyMetadataAtomic[,'matchStudy']==matchStudyName)>1) {
			atomicStudyNames = studyMetadataAtomic[studyMetadataAtomic[,'matchStudy']==matchStudyName, 'study']
			edfListNow = list()
			for (atomicStudyName in atomicStudyNames) {
				edf = data.frame(rownames(ematAtomicList[[atomicStudyName]]), ematAtomicList[[atomicStudyName]])
				rownames(edf) = NULL
				colnames(edf) = c('geneId', sampleMetadataAtomic[colnames(edf)[2:ncol(edf)], matchColname])
				edfListNow[[atomicStudyName]] = edf}
			
			edfBound = suppressWarnings(rbind_all(edfListNow))
			edfMerged = edfBound %>% group_by_(.dots=list('geneId')) %>% summarise_each(funs(mergeFunc)) %>%
				data.frame(check.names=FALSE)
			rownames(edfMerged) = edfMerged[,'geneId']
			edfMerged = edfMerged[,-1]
			
			mapping = makeMatchSampleMapping(sampleMetadataAtomic, atomicStudyNames, matchColname)
			colnames(edfMerged) = mapping[colnames(edfMerged)]
			ematList[[matchStudyName]] = as.matrix(edfMerged)
			
			idx = (sampleMetadataAtomic[,'study'] %in% atomicStudyNames) &
				(sampleMetadataAtomic[,'sample'] %in% colnames(edfMerged))
			sampleMetadataList[[matchStudyName]] = sampleMetadataAtomic[idx,]
			sampleMetadataList[[matchStudyName]][,'study'] = matchStudyName}}
	
	headFunc = function(x) x[1]
	studyMetadata = studyMetadataAtomic %>% group_by_(.dots=list('matchStudy')) %>% summarise_each(funs(headFunc)) %>%
		data.frame(check.names=FALSE)
	studyMetadata = studyMetadata[,colnames(studyMetadata)!='study']
	colnames(studyMetadata)[colnames(studyMetadata)=='matchStudy'] = 'study'
	rownames(studyMetadata) = studyMetadata[,'study']
	studyMetadata[,'matchStudy'] = studyMetadata[,'study']
	
	sampleMetadata = suppressWarnings(rbind_all(sampleMetadataList)) %>% data.frame(check.names=FALSE)
	rownames(sampleMetadata) = sampleMetadata[,'sample']
	
	result = list(ematList, studyMetadata, sampleMetadata)
	names(result) = c('ematList', 'studyMetadata', 'sampleMetadata')
	return(result)}


mergeStudyData = function(ematList, sampleMetadata, batchColname='study', covariateName=NA,
								  batchCorrection=TRUE, parPrior=TRUE) {
	# merge data and perform cross-study normalization
	geneIds = Reduce(intersect, lapply(ematList, function(x) rownames(x)))
	ematList2 = foreach(studyName=names(ematList)) %do% {ematNow = ematList[[studyName]][geneIds,]}
	if (batchCorrection) {
		# if both one-color and two-color data is present, ComBat can fail catastrophically, if data is not scaled beforehand
		ematListScaled = lapply(ematList2, function(emat) (emat - mean(emat)) / sd(emat))
		ematMerged = do.call(cbind, ematListScaled)
		if (is.na(covariateName)) {
			covariateInfo = model.matrix(~rep_len(1, ncol(ematMerged)))
		} else {
			covariateInfo = model.matrix(~sampleMetadata[colnames(ematMerged), covariateName])}
		if (length(unique(sampleMetadata[colnames(ematMerged), batchColname]))>1) {
			ematMergedNorm = ComBat(ematMerged, batch=sampleMetadata[colnames(ematMerged), batchColname],
											mod=covariateInfo, par.prior=parPrior)
		} else {
			ematMergedNorm = ematMerged}
		
	
		return(ematMergedNorm)
	} else {
		return(do.call(cbind, ematList2))}}


makeGlmnetArgs = function(metadata, foldidColname='study') {
	# construct foldid and weights for glmnet
	foldid = as.numeric(factor(metadata[,foldidColname], labels=1:length(unique(metadata[,foldidColname]))))
	names(foldid) = rownames(metadata)
	weights = length(unique(foldid)) /
		do.call("c", sapply(sapply(unique(foldid), function(x) sum(foldid==x)), function(n) rep_len(n, n), simplify=FALSE))
	names(weights) = rownames(metadata)
	return(list(foldid=foldid, weights=weights))}


crossValidateMerged = function(ematMerged, sampleMetadata, weights, alphas, nFolds=10, foldid=NA, nRepeats=3,
										 yName='class', clinVarColnames=NA, ...) {
	args = list(...)
	if (!is.null(args[['family']]) & args[['family']]=='cox') {
		y = as.matrix(sampleMetadata[colnames(ematMerged), yName])
		colnames(y) = c('time', 'status')
	} else {
		y = sampleMetadata[colnames(ematMerged), yName]}
	
	if (is.na(clinVarColnames[1])) {
		x = scale(t(ematMerged), center=TRUE, scale=FALSE)
	} else {
		clinVarTmp = data.frame(lapply(sampleMetadata[colnames(ematMerged), clinVarColnames], factor))
		clinDummy = model.matrix(~ 0 + ., data=clinVarTmp)
		x = cbind(scale(t(ematMerged), center=TRUE, scale=FALSE), clinDummy)}
	
	if (is.na(foldid[1])) {
		cvFitList = list()
		for (ii in 1:nRepeats) {
			foldid = sample(rep(seq(nFolds), length=ncol(ematMerged)))
			cvFitList[[ii]] = foreach(alpha=alphas) %do% {
				cv.glmnet(x, y, weights=weights[colnames(ematMerged)], foldid=foldid, alpha=alpha, standardize=FALSE, ...)}}
	} else {
		cvFitList = foreach(alpha=alphas) %do% {
			cv.glmnet(x, y, weights=weights[colnames(ematMerged)], foldid=foldid[colnames(ematMerged)], alpha=alpha,
						 standardize=FALSE, ...)}}
	return(cvFitList)}



predictValidationData = function(ematList, studyMetadata, sampleMetadata, discoverySampleNames, classesTrain,
											alpha, lambda, weights, batchColname='study', covariateName=NA, className='class',
											familyName='binomial', predType='response', intercept=TRUE) {
	
	discoveryStudyNames = studyMetadata[studyMetadata[,'discovery'], 'study']
	validationStudyNames = studyMetadata[studyMetadata[,'validation'], 'study']
	
	predsList = foreach(validationStudyName=validationStudyNames) %do% {
		idxValidation = sampleMetadata[,'study']==validationStudyName & (sampleMetadata[,className] %in% classesTrain)
		validationSampleNames = sampleMetadata[idxValidation, 'sample']
		
		ematListNow = ematList[c(discoveryStudyNames, validationStudyName)]
		ematMergedDiscVal = mergeStudyData(ematListNow, sampleMetadata, batchColname=batchColname, covariateName=covariateName)
            ematMergedDiscVal=ematMergedDiscVal[c('CDKN3','CHI3L1','COL4A2','DCN','DLL3','FABP7','IGFBP2','IGFBP3','LAMB1','LGALS1','LGALS3','NAMPT','PLAT','PTTG1','TIMP1','TOP2A'),]
		ematMergedDisc = ematMergedDiscVal[,discoverySampleNames]
		
		fitResult = glmnet(t(ematMergedDisc), sampleMetadata[discoverySampleNames, className], alpha=alpha, lambda=lambda,
								 weights=weights[discoverySampleNames], family=familyName, standardize=FALSE, intercept=intercept)
		preds = predict(fitResult, newx=t(ematMergedDiscVal[,validationSampleNames]), s=lambda, type=predType)}
	
	names(predsList) = validationStudyNames
	return(predsList)}


calcPredictionAuc = function(predsList, sampleMetadata, className='class') {
	nLambda = ncol(predsList[[1]])
	auc = matrix(nrow=nLambda, ncol=length(predsList))
	colnames(auc) = names(predsList)
	for (validationStudyName in names(predsList)) {
		pred = prediction(predsList[[validationStudyName]],
								matrix(rep(sampleMetadata[rownames(predsList[[validationStudyName]]), className], nLambda), ncol=nLambda),
								label.ordering=levels(sampleMetadata[,className]))
		auc[,validationStudyName] = sapply(performance(pred, 'auc')@y.values, function(x) x[[1]])}
	return(auc)}


plotCvError = function(cvFit, metaAnalysisName='metaAnalysis', size=0.4, width=5, height=3, ggplotArgs=NA) {
	df = data.frame(log(cvFit$lambda), cvFit$cvm, cvFit$cvlo, cvFit$cvup)
	colnames(df) = c('logLambda', 'cvm', 'cvlo', 'cvup')
	p = ggplot(data=df) + geom_pointrange(aes_string(x='logLambda', y='cvm', ymax='cvup', ymin='cvlo'), size=size) +
		geom_vline(xintercept=log(cvFit$lambda.min), color='blue', linetype='dashed') + labs(x='log(lambda)', y=cvFit$name)
	if (!is.na(ggplotArgs[1])) {
		for (ggplotArg in ggplotArgs) {
			p = p + ggplotArg}}
	ggsave(sprintf('%s_cv_lambda_error.pdf', metaAnalysisName), plot=p, width=width, height=height)}


makeCoefDf = function(coefResult, decreasing=TRUE, classLevels=NA) {
	if (is.list(coefResult)) {
		coefResultNonzero = foreach(coefSparse=coefResult) %do% {
			x = data.frame(rownames(coefSparse)[(coefSparse@i)+1], coefSparse[(coefSparse@i)+1], stringsAsFactors=FALSE)
			colnames(x) = c('geneId', 'coefficient')
			return(x)}
		names(coefResultNonzero) = names(coefResult)
		if (!is.na(classLevels[1])) {
			coefResult = coefResult[classLevels]
			coefResultNonzero = coefResultNonzero[classLevels]}
		
		for (ii in 1:length(coefResult)) {
			colnames(coefResultNonzero[[ii]])[2] = names(coefResult)[ii]}		
		coefDf = Reduce(function(x, y) merge(x, y, by='geneId', all=TRUE), coefResultNonzero)
		idx = do.call(order, c(coefDf[,2:ncol(coefDf)], list(decreasing=decreasing)))
		coefDf = coefDf[idx,]
		coefDf[is.na(coefDf)] = 0		
		
	} else {
		coefDf = data.frame(names(coefResult[(coefResult@i)+1,]), coefResult[(coefResult@i)+1,], stringsAsFactors=FALSE)
		colnames(coefDf) = c('geneId', 'coefficient')
		coefDf = coefDf[order(coefDf[,'coefficient'], decreasing=decreasing),]}
	rownames(coefDf) = NULL
	return(coefDf)}


plotCoefficients = function(fitResult, lambda, classLevels=NA, decreasing=FALSE, geneIdOrder=NA,
									 metaAnalysisName='metaAnalysis', width=4, height=10, ggplotArgs=NA) {
	coefResult = coef(fitResult, s=lambda)
	coefDf = makeCoefDf(coefResult, decreasing=decreasing, classLevels=classLevels)
	coefDf = coefDf[coefDf[,'geneId']!='(Intercept)',]
	
	if (!is.na(geneIdOrder[1])) {
		rownames(coefDf) = coefDf[,'geneId']
		coefDf = coefDf[geneIdOrder,]
		rownames(coefDf) = NULL}
	
	if (ncol(coefDf)==2) {
		geneSymbols = getSYMBOL(coefDf[,'geneId'], 'org.Hs.eg')
		coefDf[,'geneId'] = factor(coefDf[,'geneId'], levels=rev(coefDf[,'geneId']),
											labels=sprintf('%s (%s)', rev(geneSymbols), rev(coefDf[,'geneId'])))
		p = ggplot(data=coefDf) + geom_bar(aes_string(x='geneId', y='coefficient'), stat='identity')
		
	} else {
		if (is.na(classLevels[1])) {
			classLevels = colnames(coefDf)[2:ncol(coefDf)]}
		coefDfMolten = melt(coefDf, id.vars='geneId', variable.name='class', value.name='coefficient')
		coefDfMolten[,'class'] = factor(coefDfMolten[,'class'], levels=classLevels)
		
		geneIds = coefDf[,'geneId']
		geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')
		coefDfMolten[,'geneId'] = factor(coefDfMolten[,'geneId'], levels=rev(geneIds),
													labels=sprintf('%s (%s)', rev(geneSymbols), rev(geneIds)))
		p = ggplot(data=coefDfMolten) + facet_wrap(as.formula('~ class'), ncol=ncol(coefDf)-1) +
			geom_bar(aes_string(x='geneId', y='coefficient', fill='class'), stat='identity') + guides(fill=FALSE)}
	
	p = p + coord_flip() + labs(x='', y='Coefficient')
	if (!is.na(ggplotArgs[1])) {
		for (ggplotArg in ggplotArgs) {
			p = p + ggplotArg}}
	ggsave(filename=sprintf('%s_lambda_%.3g_gene_coef.pdf', metaAnalysisName, lambda), plot=p, width=width, height=height)}


writeCoefficients = function(fitResult, lambda, metaAnalysisName='metaAnalysis', decreasing=TRUE) {
	coefResult = coef(fitResult, s=lambda)
	coefDf = makeCoefDf(coefResult, decreasing=decreasing)
	coefDf1 = coefDf[c(which(coefDf[,'geneId']=='(Intercept)'), which(coefDf[,'geneId']!='(Intercept)')),]
	write.csv(coefDf1, file=sprintf('%s_lambda_%.3g_coef.csv', metaAnalysisName, lambda), quote=FALSE, row.names=FALSE)}


writeConfusionCrossValidation = function(cvFit, lambda, ematMerged, sampleMetadata, className='class',
													  classLevels=NA, metaAnalysisName='metaAnalysis') {
	if (is.na(classLevels[1])) {
		classLevels = names(cvFit$glmnet.fit$beta)}
	
	cvProbs = cvFit$fit.preval[,,which.min(abs(cvFit$lambda - lambda))]
	rownames(cvProbs) = colnames(ematMerged)
	colnames(cvProbs) = names(cvFit$glmnet.fit$beta)
	preds = colnames(cvProbs)[apply(cvProbs, MARGIN=1, function(x) which.max(x))]
	predsFactor = factor(preds, levels=classLevels)
	trueClasses = factor(sampleMetadata[colnames(ematMerged), className], levels=classLevels)
	confus = table(trueClasses, predsFactor)
	write.csv(confus, file=sprintf('%s_cv_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)}


writeConfusionValidation = function(predsList, lambda, sampleMetadata, className='class',
												classLevels=NA, metaAnalysisName='metaAnalysis') {
	if (is.na(classLevels[1])) {
		classLevels = colnames(predsList[[1]])}
	
	predsProb = do.call(rbind, lapply(predsList, function(x) x[,,1]))
	predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
	predsFactor = factor(predsClass, levels=classLevels)
	trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
	confus = table(trueClasses, predsFactor)
	write.csv(confus, file=sprintf('%s_val_lambda_%.3g_confusion.csv', metaAnalysisName, lambda), quote=FALSE)}


writeConfusionValidationEach = function(predsList, lambda, sampleMetadata, className='class',
													 classLevels=NA, metaAnalysisName='metaAnalysis') {
	if (is.na(classLevels[1])) {
		classLevels = colnames(predsList[[1]])}
	
	for (validationStudyName in names(predsList)) {
		predsProb = predsList[[validationStudyName]][,,1]
		predsClass = colnames(predsProb)[apply(predsProb, MARGIN=1, function(x) which.max(x))]
		predsFactor = factor(predsClass, levels=classLevels)
		trueClasses = factor(sampleMetadata[rownames(predsProb), className], levels=classLevels)
		confus = table(trueClasses, predsFactor)
		write.csv(confus, file=sprintf('%s_val_%s_lambda_%.3g_confusion.csv', metaAnalysisName, validationStudyName,
												 lambda), quote=FALSE)}}


plotValidationBeeswarm = function(predsList, sampleMetadata, metaAnalysisName, className='class', plotColors=1:2) {
	for (validationStudyName in names(predsList)) {
		dfNow = data.frame(predsList[[validationStudyName]], sampleMetadata[rownames(predsList[[validationStudyName]]), className])
		colnames(dfNow) = c('response', 'class')
		levels(dfNow[,'class']) = levels(sampleMetadata[,className])
		
		pdf(file=sprintf('%s_%s_beeswarm.pdf', metaAnalysisName, validationStudyName), width=4, height=6)
		beeswarm(as.formula('response ~ class'), data=dfNow, pch=16, col=plotColors, cex=0.7, corral='wrap', ylim=c(0, 1),
					xlab='', ylab=sprintf('Estimated probability of %s', levels(dfNow[,'class'])[2]), main=validationStudyName)
		dev.off()}}


plotValidationRoc = function(predsList, sampleMetadata, metaAnalysisName, className='class') {
	for (validationStudyName in names(predsList)) {
		pdf(file=sprintf('%s_%s_roc.pdf', metaAnalysisName, validationStudyName), width=5, height=4)
		pred = prediction(predsList[[validationStudyName]], sampleMetadata[rownames(predsList[[validationStudyName]]), className],
								label.ordering=levels(sampleMetadata[,className]))
		perfRoc = performance(pred, 'tpr', 'fpr')
		perfAuc = performance(pred, 'auc')
		plot(perfRoc, main=validationStudyName)
		text(x=0.85, y=0.05, labels=sprintf('AUC = %.3g', perfAuc@y.values[[1]]), cex=0.9)
		dev.off()}}


makeExprDfSafely = function(ematList, geneIds, scaled=TRUE) {
	exprList = foreach(geneId=geneIds) %do% {
		x = sapply(ematList, function(emat) ifelse(geneId %in% rownames(emat), list(emat[geneId,]), list(rep_len(NA, ncol(emat)))))
		if (scaled) {
			x = lapply(x, scale)}
		x = unlist(x)}
	exprDf = data.frame(exprList, unlist(sapply(ematList, function(x) colnames(x))), row.names=NULL)
	colnames(exprDf) = c(paste0('g', geneIds), 'sample')
	return(exprDf)}


plotExpressionBeeswarm = function(ematList, sampleMetadata, geneIds, className='class', metaAnalysisName='metaAnalysis',
											 scaled=TRUE, plotColors=1:2, width=10, height=6, pch=16, cex=0.6, cexLeg=1, legLoc='topleft',
											 main=NA) {
	geneIds = as.character(geneIds)
	geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')
	names(geneSymbols) = geneIds
	
	exprDf = makeExprDfSafely(ematList, geneIds, scaled)
	metadataExpr = merge(sampleMetadata, exprDf, by='sample', sort=FALSE)
	metadataExpr[,'study'] = factor(metadataExpr[,'study'], levels=unique(metadataExpr[,'study']))
	pwColors = as.character(factor(metadataExpr[,className], labels=plotColors))
	
	for (geneId in geneIds) {
		pdf(file=sprintf('%s_%s_beeswarm.pdf', metaAnalysisName, geneId), width=width, height=height)
		if (is.na(main[1])) {
			main = sprintf('%s (%s)', geneSymbols[geneId], geneId)}
		beeswarm(as.formula(sprintf('g%s ~ study', geneId)), data=metadataExpr, pwcol=pwColors, pch=pch, cex=cex,
					corral='wrap', xlab='Dataset', ylab='Expression', main=main)
		legend(legLoc, legend=levels(metadataExpr[,className]), pch=pch, cex=cexLeg, col=plotColors)
		dev.off()}}


plotExpressionHeatmapMerged = function(fitResult, lambda, ematMerged, sampleMetadata, annoNames, annoLevels, annoColors,
													clusterSamplesTogether=FALSE, geneIdOrder=NA, className='class', classLevels=NA,
													metaAnalysisName='metaAnalysis', width=8, height=8, ...) {
	coefResult = coef(fitResult, s=lambda)
	coefDf = makeCoefDf(coefResult)
	geneIds = coefDf[coefDf[,'geneId']!='(Intercept)', 'geneId']
	geneSymbols = getSYMBOL(geneIds, 'org.Hs.eg')
	geneTexts = sprintf('%s',geneIds)
	geneTexts[10]="HLA-DRB1"	
	names(geneTexts) = geneIds
	emat = ematMerged[geneIds,]
	
      # order the samples
	if (clusterSamplesTogether) {
		d = dist(t(emat))
		co = order.optimal(d, hclust(d)$merge)
		emat = emat[,co$order]
	} else {
		if (is.na(classLevels[1])) {
			classLevels = unique(sampleMetadata[colnames(ematMerged), className])}
		ematSmallList = foreach(classLevel=classLevels) %do% {
			x = emat[,colnames(emat) %in% sampleMetadata[sampleMetadata[,className]==classLevel, 'sample']]
			d = dist(t(x))
			co = order.optimal(d, hclust(d)$merge)
			x = x[,co$order]}
		emat = do.call(cbind, ematSmallList)}

	# order the genes
	if (is.na(geneIdOrder[1])) {
		d = dist(emat)
		co = order.optimal(d, hclust(d)$merge)
		emat = emat[co$order,]
		rownames(emat) = geneTexts[co$order]
	} else {
		emat = emat[geneIdOrder,]
		rownames(emat) = geneTexts[geneIdOrder]}
	
	# scale the matrix
	emat = t(scale(t(emat)))
	emat[emat>3] = 3
	emat[emat<(-3)] = -3
	
	annotation = sampleMetadata[colnames(ematMerged), annoNames, drop=FALSE]
	for (annoName in annoNames) {
		if (!is.na(annoLevels[[annoName]][1])) {
		annotation[,annoName] = factor(annotation[,annoName], levels=annoLevels[[annoName]])}}
	
	pdf(file=sprintf('%s_lambda_%.3g_heatmap.pdf', metaAnalysisName, lambda), width=width, height=height)
	pheatmap(emat, color=colorRampPalette(rev(greenred(32)))(100),
				breaks=seq(from=-3, to=3, length.out=101), cluster_rows=FALSE, cluster_cols=TRUE, treeheight_row=0,fontsize_row=2.75,
				treeheight_col=0, show_colnames=FALSE, border_color=NA, annotation=annotation, annotation_colors=annoColors, ...)
	dev.off()}



plotClassProbsCrossValidation = function(cvFit, lambda, sampleMetadata, discoveryStudyNames, discoverySampleNames,
													  className, classesTrain, metaAnalysisName, size=2, width=8, height=12, ggplotArgs=NA) {
	cvProbs = cvFit$fit.preval[,,which.min(abs(cvFit$lambda - lambda))]
	pList = list()
	for (discoveryStudyName in discoveryStudyNames) {
		discoverySampleNamesNow = discoverySampleNames[sampleMetadata[discoverySampleNames, 'study']==discoveryStudyName]
		df = data.frame(cvProbs[sampleMetadata[discoverySampleNames, 'study']==discoveryStudyName,])
		colnames(df) = names(cvFit$glmnet.fit$beta)
		df[,'study'] = discoveryStudyName
		df[,'sample'] = discoverySampleNamesNow
		df[,'trueClass'] = factor(sampleMetadata[discoverySampleNamesNow, className], levels=classesTrain)
		df[,'trueClassProb'] = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))
		
		df = df[order(df[,'trueClass'], -df[,'trueClassProb']),]
		df = do.call(rbind, lapply(classesTrain, function(x) df[df[,'trueClass']==x,]))
		
		idxTmp = c()
		for (classTrain in classesTrain) {
			if (any(df[,'trueClass']==classTrain)) {
				idxTmp = c(idxTmp, 1:(sum(df[,'trueClass']==classTrain)))}}
		df[,'idx'] = idxTmp
		
		dfMolten = melt(df, measure.vars=classesTrain, variable.name='probClass', value.name='prob')
		p = ggplot(dfMolten) + facet_grid(study ~ trueClass, scales='free_x', space='free_x') +
			geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass'), size=size) +
			labs(x='Sample', y='Probability') + theme(legend.title=element_blank())
		if (!is.na(ggplotArgs[1])) {
			for (ggplotArg in ggplotArgs) {
				p = p + ggplotArg}}
		pList[[discoveryStudyName]] = p}
	
	g = do.call(arrangeGrob, c(pList, list(nrow=length(discoveryStudyNames))))
	ggsave(sprintf('%s_cv_class_probs.pdf', metaAnalysisName), plot=g, width=width, height=height)}


plotClassProbsValidation = function(predsList, sampleMetadata, className, classesTrain, metaAnalysisName,
												size=2, width=8, height=9, ggplotArgs=NA) {
	pList = list()
	for (validationStudyName in names(predsList)) {
		df = data.frame(predsList[[validationStudyName]][,,1])
		df[,'study'] = sampleMetadata[rownames(df), 'study']
		df[,'sample'] = rownames(df)
		df[,'trueClass'] = factor(sampleMetadata[rownames(df), className], levels=classesTrain)
		df[,'trueClassProb'] = apply(df, MARGIN=1, function(x) as.numeric(x[x['trueClass']]))
		
		df = df[order(df[,'trueClass'], -df[,'trueClassProb']),]
		df = do.call(rbind, lapply(classesTrain, function(x) df[df[,'trueClass']==x,]))
		
		idxTmp = c()
		for (classTrain in classesTrain) {
			if (any(df[,'trueClass']==classTrain)) {
				idxTmp = c(idxTmp, 1:(sum(df[,'trueClass']==classTrain)))}}
		df[,'idx'] = idxTmp
		rownames(df) = NULL
		
		dfMolten = melt(df, measure.vars=classesTrain, variable.name='probClass', value.name='prob')
		p = ggplot(dfMolten) + facet_grid(study ~ trueClass, scales='free_x', space='free_x') +
			geom_point(aes_string(x='idx', y='prob', color='probClass', shape='probClass'), size=size) +
			labs(x='Sample', y='Probability') + theme(legend.title=element_blank())
		if (!is.na(ggplotArgs[1])) {
			for (ggplotArg in ggplotArgs) {
				p = p + ggplotArg}}
		pList[[validationStudyName]] = p}
	
	g = do.call(arrangeGrob, c(pList, list(nrow=length(predsList))))
	ggsave(sprintf('%s_val_class_probs.pdf', metaAnalysisName), plot=g, width=width, height=height)}
