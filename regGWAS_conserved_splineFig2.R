



library(snpStats)



#requiring the Package 
require(gam)
#ISLR package contains the 'Wage' Dataset
require(ISLR)




tabs1 <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/Supplementary_Table1_FINAL24_JRL.csv', as.is = T)


####################### regional FT
regiun <- c('germany', 'iberia', 'scandinavia')[2]


ft10R <- read.table(paste0('/media/elgon/Jesse/HerbariumSeq/AT_FT_GWAS_Lua/FT10_', regiun, '_gemma.assoc.txt'), header = T)

ft10R <- ft10R[order(ft10R$p_wald),]

ft10R$rs2 <- sub('rs_', 'Chr', ft10R$rs)


ft10R$chr <- sub('rs_', '', ft10R$rs)
ft10R$chr <- sub('_.*', '', ft10R$chr)

ft10R$ps <- sub('rs_._', '', ft10R$rs)


thinft10R <- ft10R[1,]

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


load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/IberiaTimeGWAS_14Dec2024.RDA')

topNsnp <- 25


ft10Rnor <- ft10R[ft10R$rs2 %in% rez$rs,]

rezFT10R <- rez[rez$rs %in% ft10Rnor$rs2,]

thinft10Rnor <- thinft10R[thinft10R$rs2 %in% rezFT10R$rs,]


thinft10Rnor$rs2[1:topNsnp] #need to get the snp alleles for iberia


####

ib_snpz <- thinft10Rnor$rs2[1:topNsnp]


sigSNP_ib <- thinft10Rnor$rs2[1:topNsnp][rezFT10R$ps_mixed[match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)] < 0.05 ]


rsnp_ib <- sample(rezFT10R$rs, 25)

sigSNP_ib_random <- rsnp_ib[rezFT10R[match(rsnp_ib, rezFT10R$rs), 'ps_mixed'] < 0.05]


###

####################### regional FT
regiun <- c('germany', 'iberia', 'scandinavia')[1]


ft10R <- read.table(paste0('/media/elgon/Jesse/HerbariumSeq/AT_FT_GWAS_Lua/FT16_', regiun, '_gemma.assoc.txt'), header = T)

ft10R <- ft10R[order(ft10R$p_wald),]

ft10R$rs2 <- sub('rs_', 'Chr', ft10R$rs)


ft10R$chr <- sub('rs_', '', ft10R$rs)
ft10R$chr <- sub('_.*', '', ft10R$chr)

ft10R$ps <- sub('rs_._', '', ft10R$rs)


thinft10R <- ft10R[1,]

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


load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/GermanyTimeGWAS_17Dec2024.RDA')

topNsnp <- 100


ft10Rnor <- ft10R[ft10R$rs2 %in% rez$rs,]

rezFT10R <- rez[rez$rs %in% ft10Rnor$rs2,]

thinft10Rnor <- thinft10R[thinft10R$rs2 %in% rezFT10R$rs,]


thinft10Rnor$rs2[1:topNsnp] #need to get the snp alleles for iberia


####

ge_snpz <- thinft10Rnor$rs2[1:topNsnp]


sigSNP_ge <- thinft10Rnor$rs2[1:topNsnp][rezFT10R$ps_mixed[match(thinft10Rnor$rs2[1:topNsnp], rezFT10R$rs)] < 0.05 ]


rsnp_ge <- sample(rezFT10R$rs, topNsnp)

sigSNP_ge_random <- rsnp_ge[rezFT10R[match(rsnp_ib, rezFT10R$rs), 'ps_mixed'] < 0.05]


################
######   conserved ####
library(ape)

dagenes <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/Yuxings_ConservedList_TableS1c.csv')

#read in gene info
allgeneinfo <- read.gff('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/Araport11_GFF3_genes_transposons.current.gff')
allgeneinfo <- allgeneinfo[allgeneinfo$type == 'gene',]

allgeneinfo$chr <- sub('Chr', '', allgeneinfo$seqid)

allgeneinfo$gene <- iconv(allgeneinfo$attributes, 'latin1')

allgeneinfo$gene <- sub('^ID=', '', allgeneinfo$gene)
allgeneinfo$gene <- sub(';.*', '', allgeneinfo$gene)

allgeneinfo <- allgeneinfo[!grepl('ATM', allgeneinfo$gene),]
allgeneinfo <- allgeneinfo[!grepl('ATC', allgeneinfo$gene),]

allgeneinfo$len <- (allgeneinfo$end - allgeneinfo$start)

allgeneinfo <- allgeneinfo[grep('locus_type=protein_coding', allgeneinfo$attributes),]

allgeneinfo$cons <- dagenes$Conserved[match(allgeneinfo$gene, dagenes$Ath_Gene_Locus)]


d2 <- data.frame(dagenes[,c('Ath_Gene_Locus', 'Conserved')], allgeneinfo[match(dagenes$Ath_Gene_Locus, allgeneinfo$gene), c('chr', 'start', 'end', 'len')])

lens <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/ath_gene_stat.csv')

d2$len <- lens$Longest_AA_Length[match(d2$Ath_Gene_Locus, lens$Gene_ID)]




load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/EurasiaTimeGWAS_11Dec2024.RDA')
#load('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/IberiaTimeGWAS_14Dec2024.RDA')

rez$chr <- sub('Chr', '', rez$rs)
rez$chr <- as.numeric(sub('_.*$', '', rez$chr))
rez$bp <- as.numeric(sub('^.*_', '', rez$rs))



#######
#also test GENIC enrichment. add in conservation of the gene
tmprez <- c()
for(i in 1:5){
		tmp <- rez[which(rez$chr == i),]
		tmp$genic <- 'F'
		tmpGenes <- allgeneinfo[allgeneinfo$chr == i,]
#separate ones conserved vs not
	for(j in 1:nrow(tmpGenes)){
		genz <- which(tmp$bp >= tmpGenes$start[j] & tmp$bp <= tmpGenes$end[j] )
		if(tmpGenes$cons[j] == 'yes') tmp$genic[genz] <- 'TC'
		if(tmpGenes$cons[j] == 'no') tmp$genic[genz] <- 'TN'

	}
	tmprez <- rbind(tmprez, tmp)
	}

rez <- tmprez

########### pick 1000 random snps




#########################

pdf('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/gamcheck.pdf')





#Chr3_941890 is legit ft gene


SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = ib_snpz)


EURsnp <- as(SNPm$genotypes, 'numeric')



for(i in 1:topNsnp){



plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(0,1))

	points(jitter(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)]), (EURsnp[,i]/2))


	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(colnames(EURsnp)[i] %in% sigSNP) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 2)
	if(!colnames(EURsnp)[i] %in% sigSNP) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7))


	axis(1)
	axis(2)



}

dev.off()








pdf('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/af_curvesFig2.pdf', height = 13, width = 21, pointsize = 25)

par(mfcol = c(2,3))




#Chr3_941890 is legit ft gene




SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = ib_snpz)


EURsnp <- as(SNPm$genotypes, 'numeric')



plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.3,1))

for(i in 1:length(ib_snpz)){




	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic

	pdic <- pdic - pdic[1]

	if(colnames(EURsnp)[i] %in% sigSNP_ib) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 2)
	if(!colnames(EURsnp)[i] %in% sigSNP_ib) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7))

}

	axis(1)
	axis(2)


mtext('Iberia', line = 2, font = 2)

mtext('A. Top 25 SNPs flowering time 10ºC', line = 0.5, font = 1, at  = 1940, cex = 0.75)


text(1930, -0.25, 'Enriched in temporal change', pos = 4)

title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))



	#random







SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_iberia.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = rsnp_ib)


EURsnp <- as(SNPm$genotypes, 'numeric')



plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.3,1))
for(i in 1:length(rsnp_ib)){



	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic

	pdic <- pdic - pdic[1]

	if(colnames(EURsnp)[i] %in% sigSNP_ib_random) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 2)
	if(!colnames(EURsnp)[i] %in% sigSNP_ib_random) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7))

}

	axis(1)
	axis(2)

mtext('B. Random 25 SNPs', line = 0.5, font = 1, at  = 1920, cex = 0.75)

legend('bottomright', lwd = 2, col = gray(c(0, 0.7)), legend = c('p<0.05', expression(p>=0.05)))

title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))


SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = ge_snpz)


EURsnp <- as(SNPm$genotypes, 'numeric')



plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.3,1))



for(i in 1:length(ge_snpz)){




	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic

	pdic <- pdic - pdic[1]

	if(colnames(EURsnp)[i] %in% sigSNP_ge) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 2)
	if(!colnames(EURsnp)[i] %in% sigSNP_ge) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7))

}

	axis(1)
	axis(2)


 title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))

mtext('Germany', line = 2, font = 2)

mtext('C. Top 100 SNPs flowering time 16ºC', line = 0.5, font = 1, at  = 1890, cex = 0.75)


text(1830, -0.25, 'Depleted in temporal change', pos = 4)




	#random







SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS_germany.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = rsnp_ge)


EURsnp <- as(SNPm$genotypes, 'numeric')



plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.3,1))
for(i in 1:length(rsnp_ge)){



	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic

	pdic <- pdic - pdic[1]

	if(colnames(EURsnp)[i] %in% sigSNP_ge_random) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 2)
	if(!colnames(EURsnp)[i] %in% sigSNP_ge_random) lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7))

}

	axis(1)
	axis(2)

legend('bottomright', lwd = 2, col = gray(c(0, 0.7)), legend = c('p<0.05', expression(p>=0.05)), bg = 'white')


mtext('D. Random 100 SNPs', line = 0.5, font = 1, at  = 1855, cex = 0.75)

	title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))



snp_Pick <- sample(rez$rs[rez$genic == 'TC'], 5e3)



SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = snp_Pick)

EURsnp <- as(SNPm$genotypes, 'numeric')

sigz <- snp_Pick[rez$ps_mixed[match(snp_Pick, rez$rs)] < 0.01]




plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.5,1))



snptrendz <- c()

for(i in 1:length(snp_Pick)){


if(!colnames(EURsnp)[i] %in% sigz){
	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic
#	if(mean(pdic[1:57]) > mean(pdic[58:114])) pdic <- 1 - pdic #standardize so year 1 is reference allele

	pdic <- pdic - pdic[1]

	 lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7), lwd = 0.01)
#	 snptrendz <- rbind(snptrendz, pdic)
	}
	}



for(i in 1:length(snp_Pick)){


if(colnames(EURsnp)[i] %in% sigz){
	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic
#	if(mean(pdic[1:57]) > mean(pdic[58:114])) pdic <- 1 - pdic #standardize so year 1 is reference allele

	pdic <- pdic - pdic[1]

	lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 0.6)
	 snptrendz <- rbind(snptrendz, pdic)

	}
	}


	

	axis(1)
	axis(2)



mtext('Eurasia', line = 2, font = 2)

mtext('E. 5000 SNPs from genic region of conserved', line = 0.5, font = 1, at  = 1918, cex =0.75)

title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))

text(1830, -0.45, 'Depleted in temporal change', pos = 4)


######################################


snptrendz2<- c()

snp_Pick <- sample(rez$rs[rez$genic == 'F'], 5e3)



SNPm <- read.plink(bed = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bed',  bim = '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.bim', fam= '/media/elgon/Jesse/HerbariumSeq/stem_out_GWAS.fam', na.strings = c("0", "-9"), sep = "." , select.subjects = NULL, select.snps = snp_Pick)

EURsnp <- as(SNPm$genotypes, 'numeric')

sigz <- snp_Pick[rez$ps_mixed[match(snp_Pick, rez$rs)] < 0.01]




plot.new()
plot.window(xlim = (range(tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])), ylim = c(-0.5,1))



for(i in 1:length(snp_Pick)){

if(!colnames(EURsnp)[i] %in% sigz){
	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic
#	if(mean(pdic[1:57]) > mean(pdic[58:114])) pdic <- 1 - pdic #standardize so year 1 is reference allele

	pdic <- pdic - pdic[1]

	lines(min(tdat$year):max(tdat$year), pdic, col = gray(0.7), lwd = 0.01)
#	 snptrendz2 <- rbind(snptrendz2, pdic)

	}
	}

for(i in 1:length(snp_Pick)){

if(colnames(EURsnp)[i] %in% sigz){
	tdat <- data.frame(snp = EURsnp[,i]/2, year = tabs1$Year[match(rownames(EURsnp), tabs1$Sample.ID)])

	tmp <- gam(snp ~ s(year, bs = 'ts'), family = 'binomial', data = tdat)

	pdic <- predict(tmp, data.frame(year = min(tdat$year):max(tdat$year)), type = 'response')

	if(pdic[1] > pdic[length(pdic)]) pdic <- 1 - pdic #standardize so year 1 is reference allele
#	if(mean(pdic[1:57]) > mean(pdic[58:114])) pdic <- 1 - pdic #standardize so year 1 is reference allele

	pdic <- pdic - pdic[1] #standardize so starts at zero

	 lines(min(tdat$year):max(tdat$year), pdic, col = gray(0), lwd = 0.6)
	 	 snptrendz2 <- rbind(snptrendz2, pdic)

	}

	}

	 






	axis(1)
	axis(2)



title(xlab = 'Year', ylab = expression(Delta~freq.~of~old~allele))

mtext('F. 5000 intergenic SNPs', line = 0.5, font = 1, at  = 1865, cex = 0.75)


legend('bottomright', lwd = 2, col = gray(c(0, 0.7)), legend = c('p<0.01', expression(p>=0.01)), bg = 'white')



dev.off()








