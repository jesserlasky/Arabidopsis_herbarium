
#########################################
#########################################
#########################################

#library(vcfR)

#library(SCAT)


#go over to plink to make a file for GWAS applications
#/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS
library(snpStats)
library(RColorBrewer)
#library(GMMAT)
#library(MCMCglmm)
#library(ape)
#library(GWASTools)
#library(GENESIS)
library(sommer)

eps <- 10 * .Machine$double.eps

corfun <- function(x, y, lat, lon){ #x is SNP, y is year
	tmp <- glm(cbind(x, 2-x) ~ y + lat + lon, family = 'binomial')
	if(tmp$converged & ! (any(tmp$fitted.values > 1 - eps) || any(tmp$fitted.values < eps))){
	summary(tmp)$coefficients[2,4]
	}else{
		NA
	}

	}



corfun_deu <- function(x, y, lat, lon){ #x is SNP, y is year
	lat2 <- lat^2
	lon2 <- lon^2
	tmp <- glm(cbind(x, 2-x) ~ y + lat + lon + lat2 + lon2, family = 'binomial')
	if(tmp$converged & ! (any(tmp$fitted.values > 1 - eps) || any(tmp$fitted.values < eps))){
	summary(tmp)$coefficients[2,4]
	}else{
		NA
	}

	}


wfun <- function(x, y) { #year is y, snp is x
	wilcox.test(y ~ x, method = 'sp')[[3]]
}



sommerfun <- function(SNP, env_gen_data, kinship){

		DT <- data.frame(snp = SNP, env_gen_data)

		sm_mod <- mmer(fixed = snp ~ Long + Lat + Year, random = ~vsr(genotype, Gu = kinship), data = DT, tolParInv = 1e-1)

#		tmp <- summary(sm_mod)$betas[summary(sm_mod)$betas[,'Effect'] == 'Year', 't.value']

#		p <- pt(tmp, df = summary(sm_mod)$groups[1] - 5 - 1)

		tmp <- wald.test(b = sm_mod$Beta$Estimate, Sigma = sm_mod$VarBeta, Terms = 4)

		tmp$result$chi2['P']

		}


#sommer gives t statistics, but not sure about calculating df

#tz <- rt(100000, df = 24)
#pz <- pt(tz, df = 24)
#hist(pz)

#pz[tz > 0] <- 2 * (1 - pz[tz > 0])
#pz[tz <= 0] <- 2 * (pz[tz <= 0])
#hist(pz)
#plot(tz, pz) #looks good. but, 


tabs1 <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/Supplementary_Table1_FINAL24_JRL.csv', as.is = T)

g1001clust <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/1001genomes-accessions and 1001genomes-admixture.csv')

g1001clust$luaName <- tabs1$Sample.ID[match(g1001clust$id, tabs1$Ecotype)]





###########Eurasia GWAS ###########
f <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.fam', header = FALSE, as.is = TRUE)



SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)


#write.table(data.frame(SNPm$fam[,1],SNPm$fam[,1:4], phenotype = tabs1$Year[match(rownames(SNPm$fam), tabs1$Sample.ID)]), '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.fam', col.names = F, row.names = F, sep = '\t')

#gemmaCov <- data.frame(1, scale(tabs1[match(rownames(EURsnp), tabs1$Sample.ID), c('Long', 'Lat')]))

#write.table(gemmaCov, '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_cov.txt', col.names = F, row.names = F, sep = '\t')


kin <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/plink.rel')
kinid <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/plink.rel.id')
colnames(kin) <- kinid$V1
rownames(kin) <- kinid$V1


##### Eurasia
###
EURsnp <- as(SNPm$genotypes, 'numeric')




##wilcox
wil <- apply(EURsnp, 2, wfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])


#logistic
options(warn=0, error=NULL)


logcor <- apply(EURsnp, 2, corfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)], lat = tabs1$Lat[match(rownames(EURsnp), tabs1$Sample.ID)], lon = tabs1$Long[match(rownames(EURsnp), tabs1$Sample.ID)])


######### Mixed model ###### sommer #### sommer lets me put the SNP as the response, but gemma works better and it still makes sense i think with year as response and lat/lon covariates

kinship2 <- as.matrix(kin[rownames(EURsnp), rownames(EURsnp)])

env_test <- data.frame(scale(tabs1[match(rownames(EURsnp), tabs1$Sample.ID), c('Long', 'Lat', "Year")]))

env_test$genotype = rownames(EURsnp)

sommerfun((EURsnp[,10]), env_test, kinship2)


varz <- apply(EURsnp, 2, var, na.rm = T)


somcor <- apply(EURsnp, 2, sommerfun, env_gen_data = env_test, kinship = kinship2)


#for(i in 1e4:ncol(EURsnp)) sommerfun(EURsnp[,i], env_test, kinship2)


rez <- data.frame(rs = colnames(EURsnp) , ps_wilcox = wil, ps_logist = logcor, ps_mixed = somcor)


#save(rez, file = '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/EurasiaTimeGWAS_11Dec2024.RDA')


######### Germany GWAS ####




SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)




kin <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.rel')
kinid <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.rel.id')
colnames(kin) <- kinid$V1
rownames(kin) <- kinid$V1

EURsnp <- as(SNPm$genotypes, 'numeric')




##wilcox
wil <- apply(EURsnp, 2, wfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])


#logistic
options(warn=0, error=NULL)


logcor <- apply(EURsnp, 2, corfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)], lat = tabs1$Lat[match(rownames(EURsnp), tabs1$Sample.ID)], lon = tabs1$Long[match(rownames(EURsnp), tabs1$Sample.ID)])


######### Mixed model ###### sommer #### sommer lets me put the SNP as the response, but gemma works better and it still makes sense i think with year as response and lat/lon covariates

kinship2 <- as.matrix(kin[rownames(EURsnp), rownames(EURsnp)])

env_test <- data.frame(scale(tabs1[match(rownames(EURsnp), tabs1$Sample.ID), c('Long', 'Lat', "Year")]))

env_test$genotype = rownames(EURsnp)

sommerfun((EURsnp[,10]), env_test, kinship2)



somcor <- apply(EURsnp, 2, sommerfun, env_gen_data = env_test, kinship = kinship2)




rez <- data.frame(rs = colnames(EURsnp) , ps_wilcox = wil, ps_logist = logcor, ps_mixed = somcor)


#save(rez, file = '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/GermanyTimeGWAS_17Dec2024.RDA')


####### iberia

SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)




kin <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.rel')
kinid <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.rel.id')
colnames(kin) <- kinid$V1
rownames(kin) <- kinid$V1

EURsnp <- as(SNPm$genotypes, 'numeric')




##wilcox
wil <- apply(EURsnp, 2, wfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])


#logistic
options(warn=0, error=NULL)


logcor <- apply(EURsnp, 2, corfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)], lat = tabs1$Lat[match(rownames(EURsnp), tabs1$Sample.ID)], lon = tabs1$Long[match(rownames(EURsnp), tabs1$Sample.ID)])


######### Mixed model ###### sommer #### sommer lets me put the SNP as the response, but gemma works better and it still makes sense i think with year as response and lat/lon covariates

kinship2 <- as.matrix(kin[rownames(EURsnp), rownames(EURsnp)])

env_test <- data.frame(scale(tabs1[match(rownames(EURsnp), tabs1$Sample.ID), c('Long', 'Lat', "Year")]))

env_test$genotype = rownames(EURsnp)

sommerfun((EURsnp[,10]), env_test, kinship2)



somcor <- apply(EURsnp, 2, sommerfun, env_gen_data = env_test, kinship = kinship2)




rez <- data.frame(rs = colnames(EURsnp) , ps_wilcox = wil, ps_logist = logcor, ps_mixed = somcor)


save(rez, file = '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/IberiaTimeGWAS_14Dec2024.RDA')




####### Norway  ##### haven't done yet - going to wait on LUa 

SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = NULL)




kin <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway.rel')
kinid <- read.table('/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_norway.rel.id')
colnames(kin) <- kinid$V1
rownames(kin) <- kinid$V1

EURsnp <- as(SNPm$genotypes, 'numeric')




##wilcox
wil <- apply(EURsnp, 2, wfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])


#logistic
options(warn=0, error=NULL)


logcor <- apply(EURsnp, 2, corfun, y= tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)], lat = tabs1$Lat[match(rownames(EURsnp), tabs1$Sample.ID)], lon = tabs1$Long[match(rownames(EURsnp), tabs1$Sample.ID)])


######### Mixed model ###### sommer #### sommer lets me put the SNP as the response, but gemma works better and it still makes sense i think with year as response and lat/lon covariates

kinship2 <- as.matrix(kin[rownames(EURsnp), rownames(EURsnp)])

env_test <- data.frame(scale(tabs1[match(rownames(EURsnp), tabs1$Sample.ID), c('Long', 'Lat', "Year")]))

env_test$genotype = rownames(EURsnp)

sommerfun((EURsnp[,10]), env_test, kinship2)



somcor <- apply(EURsnp, 2, sommerfun, env_gen_data = env_test, kinship = kinship2)




rez <- data.frame(rs = colnames(EURsnp) , ps_wilcox = wil, ps_logist = logcor, ps_mixed = somcor)


save(rez, file = '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/NorwayTimeGWAS_14Dec2024.RDA')








########################################################
############################################
##### Enrichment testing Loop. Make into a table to start - easier than figure.

#regional flowering time GWAS here 
#/storage/group/jrl35/default/Diana/data_analyses/AT_FT_GWAS_Lua
#moved to /media/elgon/Jesse/HerbariumSeq/AT_FT_GWAS_Lua/

##### Flowering time SNPs
ft16 <- read.csv('~/Dropbox/jesse/WaveletQstFst/Data/FT16.all.gwas.csv')

ft10 <- read.csv('~/Dropbox/jesse/WaveletQstFst/Data/FT10.all.gwas.csv')

ft10 <- ft10[order(ft10$p_wald),]
ft16 <- ft16[order(ft16$p_wald),]

ft10$rs2 <- paste0('Chr', ft10$chr, '_', ft10$ps)

ft16$rs2 <- paste0('Chr', ft16$chr, '_', ft16$ps)


thinft10 <- ft10[1,]
thinft16 <- ft16[1,]

i <- 2
while(nrow(thinft10) < 1000){
	samechr <- thinft10[thinft10$chr == ft10$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(ft10$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thinft10 <- rbind(thinft10, ft10[i,])
	}else{
		thinft10 <- rbind(thinft10, ft10[i,])
	}
	i <- i + 1
}



i <- 2
while(nrow(thinft16) < 1000){
	samechr <- thinft16[thinft16$chr == ft16$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(ft16$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thinft16 <- rbind(thinft16, ft16[i,])
	}else{
		thinft16 <- rbind(thinft16, ft16[i,])
	}
	i <- i + 1
}




####germination SNPs
gpc1 <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/MartinezBerdeja/PC1_germEE_output_lmm.assoc.txt', header = T)
gpc2 <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/MartinezBerdeja/PC2_germEE_output_lmm.assoc.txt', header = T)

gpc1$rs2 <- paste0('Chr', gpc1$rs)

gpc2$rs2 <- paste0('Chr', gpc2$rs)

gpc1 <- gpc1[order(gpc1$p_wald),]
gpc2 <- gpc2[order(gpc2$p_wald),]


thingpc1 <- gpc1[1,]
thingpc2 <- gpc2[1,]

i <- 2
while(nrow(thingpc1) < 1000){
	samechr <- thingpc1[thingpc1$chr == gpc1$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(gpc1$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thingpc1 <- rbind(thingpc1, gpc1[i,])
	}else{
		thingpc1 <- rbind(thingpc1, gpc1[i,])
	}
	i <- i + 1
}


i <- 2
while(nrow(thingpc2) < 1000){
	samechr <- thingpc2[thingpc2$chr == gpc2$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(gpc2$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thingpc2 <- rbind(thingpc2, gpc2[i,])
	}else{
		thingpc2 <- rbind(thingpc2, gpc2[i,])
	}
	i <- i + 1
}

#NOW do d13C SNPs

d13C <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/d13C_gemma.assoc.txt', header = T)

d13C$rs <- sub('rs_', '', d13C$rs)

d13C$chr <- sub('_.*', '', d13C$rs)
d13C$ps <- sub('*_', '', d13C$rs)

d13C$rs2 <- paste0('Chr', d13C$rs)

d13C <- d13C[order(d13C$p_wald),]

thind13C <- d13C[1,]

i <- 2
while(nrow(thind13C) < 1000){
	samechr <- thind13C[thind13C$chr == d13C$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(d13C$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thind13C <- rbind(thind13C, d13C[i,])
	}else{
		thind13C <- rbind(thind13C, d13C[i,])
	}
	i <- i + 1
}



#NOW do C:N SNPs

CN <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/CNratio_gemma.assoc.txt', header = T)

CN$rs <- sub('rs_', '', CN$rs)

CN$chr <- sub('_.*', '', CN$rs)
CN$ps <- sub('*_', '', CN$rs)

CN$rs2 <- paste0('Chr', CN$rs)

CN <- CN[order(CN$p_wald),]

thinCN <- CN[1,]

i <- 2
while(nrow(thinCN) < 1000){
	samechr <- thinCN[thinCN$chr == CN$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(CN$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thinCN <- rbind(thinCN, CN[i,])
	}else{
		thinCN <- rbind(thinCN, CN[i,])
	}
	i <- i + 1
}




###########

#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/EurasiaTimeGWAS_11Dec2024.RDA')
#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/GermanyTimeGWAS_17Dec2024.RDA')
#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/IberiaTimeGWAS_14Dec2024.RDA')
##load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/NorwayTimeGWAS_14Dec2024.RDA')

pop <- 'norway'


n_null <- 1e3 #number of permutations

lowtailz2 <- c()

#chose 25 or 100 or 250
for(topNsnp in c(25,100,250)){
#topNsnp <- 1e2 #signal in the top N GWAS snps (thinned)


#chose 0.05, 0.1, 0.25
#quantilez <- 0.1 #test statistic - the quantile of p-values
for(quantilez in c(0.05, 0.1, 0.25)){

nullz <- matrix(ncol = 18, nrow = n_null)


ft10nor <- ft10[ft10$rs2 %in% rez$rs,]
ft16nor <- ft16[ft16$rs2 %in% rez$rs,]
gpc1nor <- gpc1[gpc1$rs2 %in% rez$rs,]
gpc2nor <- gpc2[gpc2$rs2 %in% rez$rs,]
d13Cnor <- d13C[d13C$rs2 %in% rez$rs,]
CNnor <- CN[CN$rs2 %in% rez$rs,]

rezFT10 <- rez[rez$rs %in% ft10nor$rs2,]
rezFT16 <- rez[rez$rs %in% ft16nor$rs2,]
rezgpc1 <- rez[rez$rs %in% gpc1nor$rs2,]
rezgpc2 <- rez[rez$rs %in% gpc2nor$rs2,]
rezd13C <- rez[rez$rs %in% d13Cnor$rs2,]
rezCN <- rez[rez$rs %in% CNnor$rs2,]


sum(thinft10$rs2 %in% rezFT10$rs)
sum(thinft16$rs2 %in% rezFT16$rs)
sum(thingpc1$rs2 %in% rezgpc1$rs)
sum(thingpc2$rs2 %in% rezgpc2$rs)
sum(thind13C$rs2 %in% rezd13C$rs)
sum(thinCN$rs2 %in% rezCN$rs)

thinft10nor <- thinft10[thinft10$rs2 %in% rezFT10$rs,]
thinft16nor <- thinft16[thinft16$rs2 %in% rezFT16$rs,]
thingpc1nor <- thingpc1[thingpc1$rs2 %in% rezgpc1$rs,]
thingpc2nor <- thingpc2[thingpc2$rs2 %in% rezgpc2$rs,]
thind13Cnor <- thind13C[thind13C$rs2 %in% rezd13C$rs,]
thinCNnor <- thinCN[thinCN$rs2 %in% rezCN$rs,]


for(i in 1:n_null){

	permed10 <- sample(2:nrow(rezFT10), 1)

	nullz[i,1] <- quantile(rezFT10$ps_wilcox[c(permed10:nrow(rezFT10), 1:(permed10 - 1))][match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez)

	nullz[i,2] <- quantile(rezFT10$ps_logist[c(permed10:nrow(rezFT10), 1:(permed10 - 1))][match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez)

	nullz[i,3] <- quantile(rezFT10$ps_mixed[c(permed10:nrow(rezFT10), 1:(permed10 - 1))][match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez)



	permed16 <- sample(2:nrow(rezFT16), 1)

	nullz[i,4] <- quantile(rezFT16$ps_wilcox[c(permed16:nrow(rezFT16), 1:(permed16 - 1))][match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez)

	nullz[i,5] <- quantile(rezFT16$ps_logist[c(permed16:nrow(rezFT16), 1:(permed16 - 1))][match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez)

	nullz[i,6] <- quantile(rezFT16$ps_mixed[c(permed16:nrow(rezFT16), 1:(permed16 - 1))][match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez)


	permedgpc1 <- sample(2:nrow(rezgpc1), 1)

	nullz[i,7] <- quantile(rezgpc1$ps_wilcox[c(permedgpc1:nrow(rezgpc1), 1:(permedgpc1 - 1))][match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez)

	nullz[i,8] <- quantile(rezgpc1$ps_logist[c(permedgpc1:nrow(rezgpc1), 1:(permedgpc1 - 1))][match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez)

	nullz[i,9] <- quantile(rezgpc1$ps_mixed[c(permedgpc1:nrow(rezgpc1), 1:(permedgpc1 - 1))][match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez)





	permedgpc2 <- sample(2:nrow(rezgpc2), 1)

	nullz[i,10] <- quantile(rezgpc2$ps_wilcox[c(permedgpc2:nrow(rezgpc2), 1:(permedgpc2 - 1))][match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez)

	nullz[i,11] <- quantile(rezgpc2$ps_logist[c(permedgpc2:nrow(rezgpc2), 1:(permedgpc2 - 1))][match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez)

	nullz[i,12] <- quantile(rezgpc2$ps_mixed[c(permedgpc2:nrow(rezgpc2), 1:(permedgpc2 - 1))][match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez)



	permedd13C <- sample(2:nrow(rezd13C), 1)

	nullz[i,13] <- quantile(rezd13C$ps_wilcox[c(permedd13C:nrow(rezd13C), 1:(permedd13C - 1))][match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez)

	nullz[i,14] <- quantile(rezd13C$ps_logist[c(permedd13C:nrow(rezd13C), 1:(permedd13C - 1))][match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez)

	nullz[i,15] <- quantile(rezd13C$ps_mixed[c(permedd13C:nrow(rezd13C), 1:(permedd13C - 1))][match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez)

	permedCN <- sample(2:nrow(rezCN), 1)

	nullz[i,16] <- quantile(rezCN$ps_wilcox[c(permedCN:nrow(rezCN), 1:(permedCN - 1))][match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez)

	nullz[i,17] <- quantile(rezCN$ps_logist[c(permedCN:nrow(rezCN), 1:(permedCN - 1))][match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez)

	nullz[i,18] <- quantile(rezCN$ps_mixed[c(permedCN:nrow(rezCN), 1:(permedCN - 1))][match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez)



}


lowtailz <- c( #lower tail for temporal p-values. if there is enrichment, should be low

sum(quantile(rezFT10$ps_wilcox[match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez) > nullz[,1]),

sum(quantile(rezFT10$ps_logist[match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez) > nullz[,2]),

sum(quantile(rezFT10$ps_mixed[match(thinft10nor$rs2[1:topNsnp], rezFT10$rs)], na.rm = T, quantilez) > nullz[,3]),


sum(quantile(rezFT16$ps_wilcox[match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez) > nullz[,4]),

sum(quantile(rezFT16$ps_logist[match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez) > nullz[,5]),

sum(quantile(rezFT16$ps_mixed[match(thinft16nor$rs2[1:topNsnp], rezFT16$rs)], na.rm = T, quantilez) > nullz[,6]),



sum(quantile(rezgpc1$ps_wilcox[match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez) > nullz[,7]),

sum(quantile(rezgpc1$ps_logist[match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez) > nullz[,8]),

sum(quantile(rezgpc1$ps_mixed[match(thingpc1nor$rs2[1:topNsnp], rezgpc1$rs)], na.rm = T, quantilez) > nullz[,9]),



sum(quantile(rezgpc2$ps_wilcox[match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez) > nullz[,10]),

sum(quantile(rezgpc2$ps_logist[match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez) > nullz[,11]),

sum(quantile(rezgpc2$ps_mixed[match(thingpc2nor$rs2[1:topNsnp], rezgpc2$rs)], na.rm = T, quantilez) > nullz[,12]),


sum(quantile(rezd13C$ps_wilcox[match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez) > nullz[,13]),

sum(quantile(rezd13C$ps_logist[match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez) > nullz[,14]),

sum(quantile(rezd13C$ps_mixed[match(thind13Cnor$rs2[1:topNsnp], rezd13C$rs)], na.rm = T, quantilez) > nullz[,15]),


sum(quantile(rezCN$ps_wilcox[match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez) > nullz[,16]),

sum(quantile(rezCN$ps_logist[match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez) > nullz[,17]),

sum(quantile(rezCN$ps_mixed[match(thinCNnor$rs2[1:topNsnp], rezCN$rs)], na.rm = T, quantilez) > nullz[,18])



	) / n_null

names(lowtailz) <- c('FT10 Wilcox', 'FT10 logistic', 'FT10 mixed', 'FT16 Wilcox', 'FT16 logistic', 'FT16 mixed', 'GermPC1 Wilcox', 'GermPC1 logistic', 'GermPC1 mixed', 'GermPC2 Wilcox', 'GermPC2 logistic', 'GermPC2 mixed', 'd13C Wilcox', 'd13C logistic', 'd13C mixed', 'CN Wilcox', 'CN logistic', 'CN mixed')

lowtailz <- data.frame('p value quantile' = quantilez, nPerm = n_null, pop = pop, topN_SNP = topNsnp, t(lowtailz))


#write.csv(lowtailz, paste0('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_', pop, '_', quantilez, '_', topNsnp, 'nSNP_', n_null, 'nullP_Dec2024.csv'))


lowtailz2 <- rbind(lowtailz2, lowtailz)

}
}



write.csv(lowtailz2, paste0('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_combined_',pop, '_', n_null, 'nullP_Dec2024.csv'))

lowtailz3 <- lowtailz2
lowtailz3[,-(1:4)][lowtailz3[,-(1:4)] > 0.5] <- 1 - lowtailz3[,-(1:4)][lowtailz3[,-(1:4)] > 0.5]



### FDR

global <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_combined_1000nullP_Dec2024.csv')

global_p <- global[,-(1:5)]
global_p[global_p > 0.5] <- 1 - global_p[global_p > 0.5]
global_p <- global_p * 2

sum(p.adjust(unlist(global_p),method = 'BH') < 0.05) # = 0

germany <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_combined_germany_1000nullP_Dec2024.csv')


germany_p <- germany[,-(1:5)]
germany_p[germany_p > 0.5] <- 1 - germany_p[germany_p > 0.5]
germany_p <- germany_p * 2

sum(p.adjust(unlist(germany_p),method = 'BH') < 0.05) # = 0


iberia <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_combined_iberia_1000nullP_Dec2024.csv')

iberia_p <- iberia[,-(1:5)]
iberia_p[iberia_p > 0.5] <- 1 - iberia_p[iberia_p > 0.5]
iberia_p <- iberia_p * 2

sum(p.adjust(unlist(iberia_p),method = 'BH') < 0.05) # = 1
#d13 logistic, not clear how important. but p is zero which isn't right  - NEED more permz



norway <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_combined_norway_1000nullP_Dec2024.csv')


norway_p <- norway[,-(1:5)]
norway_p[norway_p > 0.5] <- 1 - norway_p[norway_p > 0.5]
norway_p <- norway_p * 2

sum(p.adjust(unlist(norway_p + 1e-3),method = 'BH') < 0.05) # = 0




###################################################
####################### regional FT
regiun <- c('germany', 'iberia', 'scandinavia')[3]

ft16R <- read.table(paste0('/media/elgon/Jesse/HerbariumSeq/AT_FT_GWAS_Lua/FT16_', regiun, '_gemma.assoc.txt'), header = T)

ft10R <- read.table(paste0('/media/elgon/Jesse/HerbariumSeq/AT_FT_GWAS_Lua/FT10_', regiun, '_gemma.assoc.txt'), header = T)

ft10R <- ft10R[order(ft10R$p_wald),]
ft16R <- ft16R[order(ft16R$p_wald),]

ft10R$rs2 <- sub('rs_', 'Chr', ft10R$rs)
ft16R$rs2 <- sub('rs_', 'Chr', ft16R$rs)


ft10R$chr <- sub('rs_', '', ft10R$rs)
ft10R$chr <- sub('_.*', '', ft10R$chr)

ft16R$chr <- sub('rs_', '', ft16R$rs)
ft16R$chr <- sub('_.*', '', ft16R$chr)

ft10R$ps <- sub('rs_._', '', ft10R$rs)
ft16R$ps <- sub('rs_._', '', ft16R$rs)


thinft10R <- ft10R[1,]
thinft16R <- ft16R[1,]

i <- 2
while(nrow(thinft10R) < 1000){
	samechr <- thinft10R[thinft10R$chr == ft10R$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(ft10R$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thinft10R <- rbind(thinft10R, ft10R[i,])
	}else{
		thinft10R <- rbind(thinft10R, ft10R[i,])
	}
	i <- i + 1
}



i <- 2
while(nrow(thinft16R) < 1000){
	samechr <- thinft16R[thinft16R$chr == ft16R$chr[i],]
	if(nrow(samechr) > 0){
		tmpD <- dist(c(ft16R$ps[i], samechr$ps))
		if(sum(unlist(tmpD) < 25e3) == 0) thinft16R <- rbind(thinft16R, ft16R[i,])
	}else{
		thinft16R <- rbind(thinft16R, ft16R[i,])
	}
	i <- i + 1
}



load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/NorwayTimeGWAS_14Dec2024.RDA')
#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/GermanyTimeGWAS_17Dec2024.RDA')
#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/IberiaTimeGWAS_14Dec2024.RDA')

pop <- 'norway'


n_null <- 1e4 #number of permutations



lowtailz2 <- c()

#chose 25 or 100 or 250
for(topNsnp in c(25,100,250)){
#topNsnp <- 1e2 #signal in the top N GWAS snps (thinned)


#chose 0.05, 0.1, 0.25
#quantilez <- 0.1 #test statistic - the quantile of p-values
for(quantilez in c(0.05, 0.1, 0.25)){


#topNsnp <- 1e2 #signal in the top N GWAS snps (thinned)


#quantilez <- 0.05 #test statistic - the quantile of p-values


nullz <- matrix(ncol = 6, nrow = n_null)


ft10Rnor <- ft10R[ft10R$rs2 %in% rez$rs,]
ft16Rnor <- ft16R[ft16R$rs2 %in% rez$rs,]

rezFT10R <- rez[rez$rs %in% ft10Rnor$rs2,]
rezFT16R <- rez[rez$rs %in% ft16Rnor$rs2,]

sum(thinft10R$rs2 %in% rezFT10R$rs)
sum(thinft16R$rs2 %in% rezFT16R$rs)

thinft10Rnor <- thinft10R[thinft10R$rs2 %in% rezFT10R$rs,]
thinft16Rnor <- thinft16R[thinft16R$rs2 %in% rezFT16R$rs,]


for(i in 1:n_null){

	permed10 <- sample(2:nrow(rezFT10R), 1)

	nullz[i,1] <- quantile(rezFT10R$ps_wilcox[c(permed10:nrow(rezFT10R), 1:(permed10 - 1))][match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez)

	nullz[i,2] <- quantile(rezFT10R$ps_logist[c(permed10:nrow(rezFT10R), 1:(permed10 - 1))][match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez)

	nullz[i,3] <- quantile(rezFT10R$ps_mixed[c(permed10:nrow(rezFT10R), 1:(permed10 - 1))][match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez)



	permed16 <- sample(2:nrow(rezFT16R), 1)

	nullz[i,4] <- quantile(rezFT16R$ps_wilcox[c(permed16:nrow(rezFT16R), 1:(permed16 - 1))][match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez)

	nullz[i,5] <- quantile(rezFT16R$ps_logist[c(permed16:nrow(rezFT16R), 1:(permed16 - 1))][match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez)

	nullz[i,6] <- quantile(rezFT16R$ps_mixed[c(permed16:nrow(rezFT16R), 1:(permed16 - 1))][match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez)



}


lowtailz <- c( #lower tail for temporal p-values. if there is enrichment, should be low

sum(quantile(rezFT10R$ps_wilcox[match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez) > nullz[,1]),

sum(quantile(rezFT10R$ps_logist[match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez) > nullz[,2]),

sum(quantile(rezFT10R$ps_mixed[match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)], na.rm = T, quantilez) > nullz[,3]),


sum(quantile(rezFT16R$ps_wilcox[match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez) > nullz[,4]),

sum(quantile(rezFT16R$ps_logist[match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez) > nullz[,5]),

sum(quantile(rezFT16R$ps_mixed[match(thinft16Rnor$rs2[1:topNsnp], rezFT16R$rs)], na.rm = T, quantilez) > nullz[,6])




	) / n_null

names(lowtailz) <- c('FT10 Wilcox', 'FT10 logistic', 'FT10 mixed', 'FT16 Wilcox', 'FT16 logistic', 'FT16 mixed')

lowtailz <- data.frame('p value quantile' = quantilez, nPerm = n_null, pop = pop, topN_SNP = topNsnp, region = regiun, t(lowtailz))




lowtailz2 <- rbind(lowtailz2, lowtailz)

}
}





write.csv(lowtailz2, paste0('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_regGWAS_', pop, '_combined', n_null, 'nullP_Dec2024.csv'))



####

iberia <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_regGWAS_iberia_combined10000nullP_Dec2024.csv')

iberia_p <- iberia[,-(1:6)]
iberia_p[iberia_p > 0.5] <- 1 - iberia_p[iberia_p > 0.5]
iberia_p <- iberia_p * 2

sum(p.adjust(unlist(iberia_p+1e-4),method = 'BH') < 0.05) # = 1

iberia_p <= max(unlist(iberia_p+1e-4)[p.adjust(unlist(iberia_p+1e-4),method = 'BH') < 0.05]) #max pvalue fdr sig

germany <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_regGWAS_germany_combined10000nullP_Dec2024.csv')

germany_p <- germany[,-(1:6)]
germany_p[germany_p > 0.5] <- 1 - germany_p[germany_p > 0.5]
germany_p <- germany_p * 2

sum(p.adjust(unlist(germany_p+1e-4),method = 'BH') < 0.05) # = 1

germany_p <= max(unlist(germany_p+1e-4)[p.adjust(unlist(germany_p+1e-4),method = 'BH') < 0.05]) #max pvalue fdr sig



norway <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/lowerTailQTL_regGWAS_norway_combined10000nullP_Dec2024.csv')

norway_p <- norway[,-(1:6)]
norway_p[norway_p > 0.5] <- 1 - norway_p[norway_p > 0.5]
norway_p <- norway_p * 2

sum(p.adjust(unlist(norway_p+1e-4),method = 'BH') < 0.05) # = 1

#top 25
################################### how about some go enrichments






