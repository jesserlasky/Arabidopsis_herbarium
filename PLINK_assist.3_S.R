
library(qqman)
library(RColorBrewer)
library(geosphere)
library(phyclust)
library(maps)


mepal <- colorRampPalette(brewer.pal(9, 'RdBu')[-(4:6)])



#
admx <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/1001genomes-accessions and 1001genomes-admixture.csv')

tabs1 <- read.csv('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Data/Supplementary_Table1_FINAL24_JRL.csv', as.is = T)

luafam <- read.table('/media/elgon/Lua/aDNA_2020/mpileup_FINAL_24/hap_grab_24/stem_out.tfam')





tabs1$herbarium <- NA
tabs1$herbarium[! tabs1$Herbarium.Fresh.tissue.Data.publicly.avaliable %in% c('Data Publicly avaliable', 'Fresh Tissue')] <- 'Y'

tabs1$herbarium[tabs1$Data.Set == 'AK' & grepl('hb', tabs1$Sample.ID)] <- 'Y'




badz1001 <- read.csv('~/Dropbox/jesse/Arabidopsis/CAREER/1001GRegMap_stocks_in_question.csv')

badz <- c(tabs1$Sample.ID[!is.na(tabs1$Ecotype) & tabs1$Ecotype %in% badz1001$ecotype_id])



badQCt <- data.frame( badz, 1, 0, 0, 0, -9)



###########################
### genetic distance
pdist <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_1001Gbadrem.mdist')
pdistID <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_1001Gbadrem.mdist.id', as.is = T)


geoD1 <- as.matrix(dist(tabs1[match(pdistID$V1, tabs1$Sample.ID) ,c('Long', 'Lat')]))

rownames(geoD1) <- pdistID$V1
colnames(geoD1) <- pdistID$V1

#bad ones
badz <- c(badz, 'NOR_1968', 'RUS_2003b_1001G', 'ETH_1838_S52', 'ETH_2018i')



#make regions for colors
tabs1$region <- NA
tabs1$region[tabs1$Country %in% c('Algeria', 'Morocco')] <- 'NW Africa'
tabs1$region[tabs1$Country %in% c('Austria', 'Belgium', 'Czechia', 'Denmark', 'France', 'Germany', 'Poland', 'Slovakia', 'Switzerland', 'The Netherlands')] <- 'C Eur'
tabs1$region[tabs1$Country %in% c('Cyprus', 'Italy', 'Lebanon')] <- 'Ib'
tabs1$region[tabs1$Country %in% c('Democratic Republic of the Congo', 'Eritrea', 'Ethiopia', 'Kenia', 'Rwanda', 'Tanzania', 'Uganda', 'South Africa')] <- 'E & S Africa'
tabs1$region[tabs1$Country %in% c('Great Britain')] <- 'GB'
tabs1$region[tabs1$Country %in% c('Finland', 'Norway', 'Sweden')] <- 'Scand'
tabs1$region[tabs1$Country %in% c('Spain', 'Portugal')] <- 'Ib'

unique(tabs1$Country[is.na(tabs1$region)])

tabs1$region[is.na(tabs1$region)] <- 'East'



###########################



######### pick morrocan near herbarium
dz <- as.matrix(dist(tabs1[,c('Long', 'Lat')]))
md <- min(dz[which(tabs1$Sample.ID == 'MAR_1981'),][dz[which(tabs1$Sample.ID == 'MAR_1981'),] > 0],na.rm = T)
tabs1$Sample.ID[which( dz[which(tabs1$Sample.ID == 'MAR_1981'),] == md)]







#keep MAR_Taz16 MAR_Bbe0 MAR_Til2 

#MAR_Khe's MAR_Agl's, MAR_Elk's, MAR_Azr's

badz <- c(badz, tabs1$Sample.ID[grep('MAR_', tabs1$Sample.ID)])
badz <- badz[- grep('MAR_2007_Taz16', badz)]
#badz <- badz[- grep('MAR_Bbe0', badz)]
badz <- badz[- grep('MAR_2007_Til2', badz)]
badz <- badz[- grep('MAR_1981', badz)]
badz <- badz[- grep('MAR_2007_Khe', badz)]
badz <- badz[- grep('MAR_2007_Agl', badz)]
badz <- badz[- grep('MAR_2007_Elk', badz)]
badz <- badz[- grep('MAR_2007_Azr', badz)]


#pick norwegian modern to match herbarium

norHerb <- tabs1[tabs1$Country == 'Norway' & tabs1$Herbarium.Fresh.tissue.Data.publicly.avaliable == 'University of Oslo' & tabs1$Sample.ID %in% pdistID$V1,]

norNew <- tabs1[tabs1$Country == 'Norway' & tabs1$Herbarium.Fresh.tissue.Data.publicly.avaliable == 'Fresh Tissue',]

addNor <- c()
for(i in 1:nrow(norHerb)){
	tmp <- as.character(norNew$Sample.ID)[which(geoD1[norHerb$Sample.ID[i], as.character(norNew$Sample.ID)] == min(geoD1[norHerb$Sample.ID[i], as.character(norNew$Sample.ID)], na.rm = T))]
	if(!tmp[1] %in% addNor){
		 addNor <- c(addNor, tmp[1])
		 	}else{
		 		tmpD <- geoD1[norHerb$Sample.ID[i], as.character(norNew$Sample.ID)]
		 		tmp <- as.character(norNew$Sample.ID)[which(geoD1[norHerb$Sample.ID[i], as.character(norNew$Sample.ID)] == min(geoD1[norHerb$Sample.ID[i], as.character(norNew$Sample.ID)][! as.character(norNew$Sample.ID) %in% addNor], na.rm = T))[1]]
		 	 addNor <- c(addNor, tmp[1])

		 	}
}



norNew$Site <- sub('_.$', '', norNew$Sample.ID)
norNew$Site <- sub('_..$', '', norNew$Site)

addNorSite <- sub('_.$', '', addNor)
addNorSite <- sub('_..$', '', addNorSite)


for(i in 1:length(unique(norNew$Site[!norNew$Site %in% addNorSite]))){
 	 addNor <- c(addNor, norNew$Sample.ID[norNew$Site == unique(norNew$Site[!norNew$Site %in% addNorSite])[i]][1])
	}

badz <- c(badz, norNew$Sample.ID[!norNew$Sample.ID %in% addNor])


badQCt <- data.frame( badz, 1, 0, 0, 0, -9)

#write.table(badQCt, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', col.names = F, row.names = F, quote = F)



#########
#remove african ones and the madeiran for eurasian structure analyses
badQCt <- read.table( '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', as.is = T)

badz <- badQCt[,1]

badz <- c(badz, tabs1$Sample.ID[tabs1$region %in% c('NW Africa', 'E & S Africa')], 'PRT_1954')

badQCt <- data.frame( badz, 1, 0, 0, 0, -9)

#write.table(badQCt, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_eurasianFocus.txt', col.names = F, row.names = F, quote = F)

########




pdist <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_allrange.mdist')
pdistID <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_allrange.mdist.id', as.is = T)

############ combined figure ###############

#panels:
#map
#NJ tree (large)
#global MDS plot
#eurasia MDS plot
#Structure plot K = 3


source('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Code/MyPlotNJ.R')


#admix
imss <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_allrange.imiss', header = T, as.is = T)


tbl <- read.table("~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_allrange_pruned.6.Q")

rownames(tbl) <- imss$FID

win <- rep(NA, nrow(tbl))

for(i in 1:nrow(tbl)) win[i] <- which(tbl[i,] == max(tbl[i,]))

imss$win <- win

dominant <- tbl[1:nrow(tbl),win]


tbl2 <- c()

for(i in unique(win)) tbl2 <- rbind(tbl2, tbl[win == i,][order(-tbl[win == i,i]),])


win2 <- rep(NA, nrow(tbl2))

for(i in 1:nrow(tbl2)) win2[i] <- which(tbl2[i,] == max(tbl2[i,]))


tbl3 <- c()

tbl_region <- tabs1$region[match(rownames(tbl), tabs1$Sample.ID)]

for(i in c("E & S Africa", 'NW Africa', 'Ib', 'GB','C Eur', 'Scand',  'East')){
	tmp <- table(win[tbl_region == i])
	tdom <- as.numeric(names(tmp)[order(-tmp)][1])
	tbl3 <- rbind(tbl3, tbl[tbl_region == i,][order(-tbl[tbl_region == i, tdom]),])
	}



####

tabs1b <- tabs1[match(pdistID$V1, tabs1$Sample.ID),]







pdf('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_NJ_mds_regions_1Dec2024.pdf', height = 7.5, pointsize = 14)

#par(mfrow = c(2,2))
par(mar = c(3,3,1,2))

par(plt = c(0.05,0.45,0.63,0.999))


plot(tabs1b$Long, tabs1b$Lat, cex =0, col = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1b$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], pch = c(19, 17)[(tabs1b$herbarium %in% c('Y'))+1], lwd = 0.5, bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '')

map(add = T, col = gray(0.8), lwd = 0.5)

mtext('A.', at = -20, line = -1.5)



points(tabs1b$Long, tabs1b$Lat, cex =0.4, bg = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1b$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], pch = c(21, NA)[(tabs1b$herbarium %in% c('Y'))+1], lwd = 0.25, col = 1)

points(tabs1b$Long, tabs1b$Lat, cex =0.4, bg = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1b$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], pch = c(NA, 24)[(tabs1b$herbarium %in% c('Y'))+1], lwd = 0.5, col = 1)


legend('bottomright', c('Herbarium', 'Field/Stockcenter'), pch = c(17, 19), cex = 0.5, bty = 'n')


par(new = T)
par(plt = c(0.45,0.99,0.35,0.999))

#NJ tree
plotnj_JRL(nj(as.matrix(pdist)), show.tip.label = F, edge.color = gray(0.5), edge.width = 0.1, meRotate = -70)# cex =0.15, show.tip.label = T, col = gray(0.5), lwd = 0.5)#, type = 'fan') #, tip.color = terrain.colors(200)[(combP$year-1866)])

#geo tip labels
tiplabels(pch = c(19, 17)[(tabs1b$herbarium %in% c('Y'))+1], col = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1b$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))])#, bg = pdistID$bgcol)



tiplabels(text = pdistID$V1, frame = 'none', bg = NA, cex = 0.1, col = gray(0.2))


mtext('B.', line = -1.5, at = 0.005)

text(0.025, 0.042, 'Colors in (A-D)', cex = 0.5)
text(0.026, 0.041, ' = geographic groups', cex = 0.5)


##admix
par(new = T)
par(plt = c(0.6,0.95,0.07,0.28))





admx_sp <- table(tbl_region)[c("E & S Africa", 'NW Africa', 'Ib', 'GB','C Eur', 'Scand',  'East')]

mp <- barplot(t(as.matrix(tbl3)), col= c(brewer.pal(11, 'PiYG')[c(2,4)], gray(0.3), gray(0.7), brewer.pal(11, 'BrBG')[c(5,2)])[1:6][1:max(win)], xlab="Individual #", ylab="", border=NA,  cex.axis = 0.001, cex.lab = 0.001, las = 2, cex = 0.25, yaxt ='n', space = c(rep(0, admx_sp[1] ), 4, rep(0, admx_sp[2] - 1 ), 4, rep(0, admx_sp[3] - 1), 4, rep(0, admx_sp[4] - 1), 4, rep(0, admx_sp[5] - 1), 4, rep(0, admx_sp[6] - 1), 4, rep(0, admx_sp[7] - 1)), legend.text = F, names.arg =rep('', nrow(tbl3)))#, 


points(mp[which(rownames(tbl3) %in% tabs1b$Sample.ID[tabs1b$herbarium %in% c('Y')])], rep(0, length(which(rownames(tbl3) %in% tabs1b$Sample.ID[tabs1b$herbarium %in% c('Y')]))), pch = 17, cex = 0.25)


mtext('S & E', at = 10, cex = 0.5, line = 0.5)
mtext('Afr.', at = 20, cex = 0.5, line = 0)

mtext('N', at = 70, cex = 0.5, line = 0.5)
mtext('Afr.', at = 72, cex = 0.5, line = 0)

mtext('Iber.', at = 110, cex = 0.5, line = 0.5)

mtext('Brit.', at = 159, cex = 0.5, line = 0.5)

mtext('Central', at = 240, cex = 0.5, line = 0.5)
mtext('Europe', at = 240, cex = 0.5, line = 0)

mtext('Scand.', at = 350, cex = 0.5, line = 0.5)

mtext('E Europe', at = 440, cex = 0.5, line = 0.5)
mtext('& Asia', at = 430, cex = 0.5, line = 0.)

mtext('Ancestry', line = 2, side = 2, cex = 0.8)

mtext('Colors = ADMIXTURE clusters', line = 0, side = 1, cex = 0.6)
mtext('Geographic groups as in (A)', line = 0.8, side = 1, cex = 0.6)

mtext('E.', line = 1, side = 3, at = -50)

#axis(1)
axis(2, at = c(0,1))

par(xpd = T)

polygon(c(-150, 500, 500, -150), c(-0.3, -0.3, 1.35, 1.35))


par(xpd = F)

p_mds <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/globalMDS_2.mds', header = T)


all.equal(tabs1b[,1], p_mds[,1])



par(new = T)
par(plt = c(0.15,0.36,0.35,0.57))





plot(-p_mds$C2, -p_mds$C1, cex = 0.8, pch = c(21, 24)[(tabs1b$herbarium %in% c('Y'))+1], bg = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1b$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', col = 1)


mtext('C.', line = 0, side = 3, at = -0.02)

#text(-p_mds$C1, -p_mds$C2,  p_mds$FID, cex = 0.25)

text(-0.012, 0.007, 'N Africa', cex = 0.6)

text(-0.01, -0.015, 'S & E Africa', cex = 0.6)
text(0.002, -0.0015, 'Eurasia', cex = 0.6)

mtext('MDS 1', side = 2, line = 0.5, cex = 0.8)
mtext('MDS 2', side = 1, line = 0.5, cex = 0.8)
mtext('Global', side = 3, line = 0)





p_mds_eur <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/EurasiaMDS.mds', header = T)

tabs1_eur <- tabs1[match(p_mds_eur$FID, tabs1$Sample.ID),]



par(new = T)
par(plt = c(0.15,0.4,0.05,0.23))





plot(-p_mds_eur$C1, p_mds_eur$C2, cex = 0.8, pch = c(21, 24)[(tabs1_eur$herbarium %in% c('Y'))+1], bg = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1_eur$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', col = 1)


text(-0.035, 0, 'Iberian relict', cex = 0.6, pos = 4)
text(-0.005, -0.016, 'Asia', cex = 0.6, pos = 4)
text(-0.015, 0.012, 'Norway', cex = 0.6, pos = 4)

#text(-p_mds_eur$C1, p_mds_eur$C2,  p_mds_eur$FID, cex = 0.25)



mtext('MDS 2', side = 2, line = 0.5, cex = 0.8)
mtext('MDS 1', side = 1, line = 0.5, cex = 0.8)
mtext('Eurasia', side = 3, line = 0)


mtext('D.', line = 0, side = 3, at = -0.035)





dev.off()



####################
#keep african ones african structure analyses
badQCt <- read.table( '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', as.is = T)

badz <- badQCt[,1]


badz <- unique(c(badz, tabs1$Sample.ID[!tabs1$region %in% c('E & S Africa')]))


badQCt <- data.frame( badz, 1, 0, 0, 0, -9)

#write.table(badQCt, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_seAfricanFocus.txt', col.names = F, row.names = F, quote = F)



pdist2T <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_seAfrica.mdist')
rownames(pdist2T) <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_seAfrica.mdist.id', as.is = T)[,1]





#admix
imss <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_seAfrica.imiss', header = T, as.is = T)

tbl <- read.table("~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_seAfrica_pruned.4.Q")

tblOld <- read.table("~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_Afr_pruned.4.Q")

tbl5 <- read.table("~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/stem_outLowMiss_seAfrica_pruned.5.Q") #to check the schimper


rownames(tbl) <- imss$FID
rownames(tbl5) <- imss$FID


win <- rep(NA, nrow(tbl))

for(i in 1:nrow(tbl)) win[i] <- which(tbl[i,] == max(tbl[i,]))

imss$win <- win

dominant <- tbl[1:nrow(tbl),win]

#tbl2 <- tbl[order(win),]
tbl2 <- c()

for(i in unique(win)) tbl2 <- rbind(tbl2, tbl[win == i,][order(-tbl[win == i,i]),])


win2 <- rep(NA, nrow(tbl2))

for(i in 1:nrow(tbl2)) win2[i] <- which(tbl2[i,] == max(tbl2[i,]))



plot(tabs1$Long[match(imss$FID, tabs1$Sample.ID)], tabs1$Lat[match(imss$FID, tabs1$Sample.ID)])
text(tabs1$Long[match(imss$FID, tabs1$Sample.ID)], tabs1$Lat[match(imss$FID, tabs1$Sample.ID)], win, col = 2, cex = 2)
 map(add = T)

#1 is SE of rift
#2 is SE Africa
#3 is nw of rift
#4 is elgon/mt kenya



pdf('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/plink_nj_EAfrica_comb_3Dec2024.pdf', width = 7, pointsize = 15)


par(plt = c(0.05,0.55,0.45,0.95))



plotnj_JRL(nj(as.matrix(pdist2T)), show.tip.label = F, edge.color = gray(0.5), edge.width = 0.1, meRotate = -150)# cex =0.15, show.tip.label = T, col = gray(0.5), lwd = 0.5)#, type = 'fan') #, tip.color = terrain.colors(200)[(combP$year-1866)])



tiplabels(pch = c(19,17)[(rownames(pdist2T) %in% tabs1$Sample.ID[tabs1$herbarium %in% c('Y')]) + 1], col = brewer.pal(11, 'Spectral')[1])#, col = c("red","blue","darkmagenta","forestgreen","hotpink","black", "brown","darkcyan")[match('E & S Africa', c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], bg = 'white')






##admix
par(new = T)
par(plt = c(0.45,0.95,0.07,0.28))


mp <- barplot(t(as.matrix(tbl2)), col= c(brewer.pal(8, 'Accent')[7],brewer.pal( 8, 'Dark2')[6:8])[1:max(win)], xlab="Individual #", ylab="", border=NA,  cex.axis = 0.001, cex.lab = 0.001, las = 2, cex = 0.05, yaxt ='n', 
		space = c(rep(0, sum(win2 == 2) ), 1, rep(0, sum(win2 == 1) - 1 ), 1, rep(0, sum(win2 == 3) - 1), 1, rep(0, sum(win2 == 4) - 1)), names.arg = rep('', nrow(tbl2)))

points(mp[which(rownames(tbl2) %in% tabs1$Sample.ID[tabs1$herbarium %in% c('Y')])], rep(0, length(which(rownames(tbl2) %in% tabs1$Sample.ID[tabs1$herbarium %in% c('Y')]))), pch = 17, cex = 0.5)


mtext('SE Africa,', at = 8, cex = 0.5, line = 1)
mtext('Simien,', at = 8, cex = 0.5, line = 0.5)
mtext('low elev. Horn', at = 8, cex = 0.5)

mtext('Ethiopia', at = 28, cex = 0.5, line = 0.5)
mtext('SE of Rift', at = 28, cex = 0.5)

mtext('Choke', at = 40, cex = 0.5, line = 0.5)
mtext('', at = 35, cex = 0.5)

mtext('Mt Elgon', at = 50, cex = 0.5, line = 0.5)
mtext('& Mt Kenya', at = 50, cex = 0.5)







p_mds_afr <- read.table('~/Dropbox/jesse/Arabidopsis/HerbariumSeq/Results/seAfricaMDS.mds', header = T)

tabs1_afr <- tabs1[match(p_mds_afr$FID, tabs1$Sample.ID),]



par(new = T)
par(plt = c(0.1,0.4,0.05,0.33))





plot(p_mds_afr$C1, -p_mds_afr$C2, cex = 0.8, pch = c(21, 24)[(tabs1_afr$herbarium %in% c('Y'))+1], bg = brewer.pal(11, 'Spectral')[c(2,8,9,1,7,11,3,10)][match(tabs1_afr$region, c('NW Africa', 'C Eur', 'E Med', 'E & S Africa', 'GB', 'Scand', 'Ib', 'East'))], xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', col = 1)






mtext('MDS 2', side = 2, line = 0.5, cex = 0.8)
mtext('MDS 1', side = 1, line = 0.5, cex = 0.8)



legend('topleft', c('Herbarium', 'Field/Stockcenter'), pch = c(17, 19), cex = 0.5, bty = 'n')








dev.off()





## for regional gwas 


######### 
badQCt <- read.table( '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', as.is = T)

badz <- badQCt[,1]

#make one for germany

badz2ger <- c(badz, tabs1$Sample.ID[which(tabs1$Country != 'Germany')])

badQCt2ger <- data.frame( badz2ger, 1, 0, 0, 0, -9)
write.table(badQCt2ger, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forGermanyFocus.txt', col.names = F, row.names = F, quote = F)





## spain & portugal
badQCt <- read.table( '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', as.is = T)

badz <- badQCt[,1]

badz2 <- c(badz, tabs1$Sample.ID[which(tabs1$Long < -10)])

badz2iberia <- c(badz2, tabs1$Sample.ID[which(! tabs1$Country %in% c('Spain', 'Portugal'))])

badQCt2iberia <- data.frame( badz2iberia, 1, 0, 0, 0, -9)
write.table(badQCt2iberia, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forIberiaFocus.txt', col.names = F, row.names = F, quote = F)




#make one for norway

badQCt <- read.table( '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps3_allrange.txt', as.is = T)

badz <- badQCt[,1]

badz2nor <- c(badz, tabs1$Sample.ID[which(tabs1$Country != 'Norway')])

badQCt2nor <- data.frame( badz2nor, 1, 0, 0, 0, -9)
write.table(badQCt2nor, '~/Dropbox/jesse/Arabidopsis/HerbariumSeq/badsamps_3_forNorwayFocus.txt', col.names = F, row.names = F, quote = F)


























