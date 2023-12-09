install.packages("Rtsne")
install.packages("ggplot2")


#引用包
library(Rtsne)
library(ggplot2)
setwd("C:\\Users\\Kangtao Wang\\Desktop\\pyroptosis and PC\\model and validation\\03.PCA t-SNE")    #设置工作目录

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
	#读取输入文件,提取数据
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=rt[c(3:(ncol(rt)-2))]
	risk=rt[,"risk"]

    #PCA分析
	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
	#绘制PCA图
	pdf(file=pcaFile, height=4.5, width=5.5)       #保存输入出文件
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
	

	#t-SNE分析
	tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
	tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	
	#绘制tSNE图
	pdf(file=tsneFile, height=4.5, width=5.5)       #保存输入出文件
	p=ggplot(data = tsne, aes(tSNE1, tSNE2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
}
bioPCA(inputFile="TCGA.txt", pcaFile="TCGA.PCA.pdf", tsneFile="TCGA.t-SNE.pdf")
bioPCA(inputFile="PACA.txt", pcaFile="PACA.PCA.pdf", tsneFile="PACA.t-SNE.pdf")
bioPCA(inputFile="GSE57495.txt", pcaFile="GSE57495.PCA.pdf", tsneFile="GSE57495.t-SNE.pdf")