
#install.packages("pheatmap")


#引用包
library(pheatmap)
setwd("C:\\Users\\Kangtao Wang\\Desktop\\pyroptosis and PC\\model and validation\\02.riskplot")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null)
  {
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
	rt$riskScore[rt$riskScore>quantile(rt$riskScore,0.99)]=quantile(rt$riskScore,0.99)
	rt$risk=factor(rt$risk, levels=c("low", "high"))
	rt=rt[order(rt$riskScore),]      #按照风险得分对样品排序
		
	#绘制风险曲线
	riskClass=rt[,"risk"]
	lowLength=length(riskClass[riskClass=="low"])
	highLength=length(riskClass[riskClass=="high"])
	lowMax=max(rt$riskScore[riskClass=="low"])
	line=rt[,"riskScore"]
	pdf(file=riskScoreFile, width=6, height=5)
	plot(line, type="p", pch=20,
		 xlab="Patients (increasing risk socre)", ylab="Risk score",
		 col=c(rep("blue",lowLength),rep("red",highLength)) )
	abline(h=lowMax,v=lowLength,lty=2)
	legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	dev.off()
		
	#绘制生存状态图
	color=as.vector(rt$fustat)
	color[color==1]="red"
	color[color==0]="blue"
	pdf(file=survStatFile, width=6, height=5)
	plot(rt$futime, pch=19,
		 xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
		 col=color)
	legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
	abline(v=lowLength,lty=2)
	dev.off()
}

#调用函数，绘制风险曲线
bioRiskPlot(inputFile="TCGA.txt",
            riskScoreFile="TCGA.riskScore.pdf",
            survStatFile="TCGA.survStat.pdf")

bioRiskPlot(inputFile="PACA.txt",
            riskScoreFile="PACA.riskScore.pdf",
            survStatFile="PACA.survStat.pdf")

bioRiskPlot(inputFile="GSE57495.txt",
            riskScoreFile="GSE57495.riskScore.pdf",
            survStatFile="GSE57495.survStat.pdf")