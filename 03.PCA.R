install.packages("Rtsne")
install.packages("ggplot2")


#���ð�
library(Rtsne)
library(ggplot2)
setwd("C:\\Users\\Kangtao Wang\\Desktop\\pyroptosis and PC\\model and validation\\03.PCA t-SNE")    #���ù���Ŀ¼

bioPCA=function(inputFile=null, pcaFile=null, tsneFile=null){
	#��ȡ�����ļ�,��ȡ����
	rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
	data=rt[c(3:(ncol(rt)-2))]
	risk=rt[,"risk"]

    #PCA����
	data.pca=prcomp(data, scale. = TRUE)
	pcaPredict=predict(data.pca)
	PCA = data.frame(PC1 = pcaPredict[,1], PC2 = pcaPredict[,2],risk=risk)	
	#����PCAͼ
	pdf(file=pcaFile, height=4.5, width=5.5)       #����������ļ�
	p=ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = risk)) +
		scale_colour_manual(name="Risk",  values =c("red", "blue"))+
	    theme_bw()+
	    theme(plot.margin=unit(rep(1.5,4),'lines'))+
	    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	print(p)
	dev.off()
	

	#t-SNE����
	tsneOut=Rtsne(data, dims=2, perplexity=10, verbose=F, max_iter=500,check_duplicates=F)
	tsne=data.frame(tSNE1 = tsneOut$Y[,1], tSNE2 = tsneOut$Y[,2],risk=risk)	
	#����tSNEͼ
	pdf(file=tsneFile, height=4.5, width=5.5)       #����������ļ�
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