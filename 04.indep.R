#install.packages('survival')


library(survival)       #引用包
setwd("C:\\Users\\Kangtao Wang\\Desktop\\pyroptosis and PC\\model and validation\\04.indep")     #设置工作目录

############绘制森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#读取输入文件
	rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrLow[hrLow<0.001]=0.001
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#输出图形
	pdf(file=forestFile, width=6.5, height=4.8)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#绘制森林图左边的临床信息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	#绘制森林图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex=2 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,LOGindex^a1)
	dev.off()
}
############绘制森林图函数############

#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
	
	#数据合并
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	
	#单因素独立预后分析
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

	#多因素独立预后分析
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
	rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

#调用函数，进行独立预后分析
indep(riskFile="GSE57495.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox.txt",
      multiOutFile="multiCox.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")