library(tidyverse)
library(tibble)
library(data.table)
library(survival)
library(survminer)
library(cowplot)
library(timeROC)
library(pROC)

# riskScore ---------------------------------------------------------------

setwd("C:\\Users\\Kangtao Wang\\Desktop\\pyroptosis and PC\\model and validation\\01.ROC")
surv_expr <- fread("r.tcga.paad._surv.expr.txt", data.table = F)
rownames(surv_expr) <- surv_expr$id
surv_expr <- surv_expr[, -1]
surv_expr[1:5, 1:5]

names(surv_expr)[1:2] <- c("time", "status")

## month
surv_expr$time <- surv_expr$time / 30
surv_expr <- surv_expr %>% dplyr::filter(time != 0)

coef <- read.csv("Coef_result_lasso.csv") %>% dplyr::filter(final_gene %in% intersect(.$final_gene, names(surv_expr)))
final_gene <- coef$final_gene
Active.Coefficients <- coef$Active.Coefficients
expr_train <- surv_expr[, final_gene]
risk.score <- c()
for (i in 1:nrow(expr_train)) {
  risk.score[i] <- sum(Active.Coefficients * expr_train[i, ])
}

df <- surv_expr %>% dplyr::select(time, status) %>% mutate(riskScore = risk.score)
rownames(df) <- rownames(surv_expr)
df$group <- ifelse(df$riskScore >= median(df$riskScore), "high", "low")
try<-cbind(df,expr_train)
write.csv(try, "riskScore.csv")

# survplot ----------------------------------------------------------------


surv_df <- df
sfit <- survfit(Surv(time, status) ~ group, data = surv_df)
diff <- survdiff(formula = Surv(time, status) ~ group, data = surv_df, rho = 0)
pval <- pchisq(diff$chisq, length(diff$n) - 1, lower.tail = FALSE)

gg <- ggsurvplot(sfit,
                 risk.table = TRUE,
                 conf.int = F, # 置信区间
                 palette = c("#DC0000FF","#0000ff"),
                 pval.method = F,
                 pval = ifelse(pval < 0.001, "p < 0.001", paste0("p = ", pval %>% round(digits = 3))),
                 ncensor.plot = F,
                 surv.median.line = "hv",
                 linetype = "strata", xlab = "Time(Months)",
                 legend = c(0.8, 0.95),
                 legend.title = "",
                 legend.labs = c("high-risk", "low-risk")
)
gg2 <- plot_grid(gg$plot + theme(legend.text = element_text(size = 12)), gg$table, ncol = 1, rel_widths = 1, rel_heights = c(1, .4), align = "hv")
pdf("survplot.pdf", width = 5, height = 6)
print(gg2)
dev.off()

# timeROC -----------------------------------------------------------------


mycol <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")

rt <- surv_df
rt$time <- rt$time / 12

ROC.DSST <- timeROC(
  T = rt$time, delta = rt$status,
  marker = rt$riskScore, cause = 1,
  weighting = "marginal",
  times = c(1, 3, 5), ROC = TRUE
)

pdf("timeROC.pdf", 4, 4)
plot(c(0, 1), c(0, 1), type = "l", xlab = "False positive rate", ylab = "Ture positive rate", main = "")
for (i in 1:length(ROC.DSST$times)) {
  los <- lowess(ROC.DSST$FP[, i], y = ROC.DSST$TP[, i], f = 1 / 3, iter = 100)
  los$x <- c(0, los$x, 1)
  los$y <- c(0, los$y, 1)
  lines(los, col = mycol[i], lwd = 3)
}
lb <- "year(s)"
legend("bottomright", paste0(
  ROC.DSST$times,
  "-", lb,
  " "
  # ,' AUC='
  , round(ROC.DSST$AUC, 3)
), bty = "n", lwd = 3, col = mycol[1:3])
dev.off()


pdf("ROC.pdf", width = 4, height = 4)
rocobj <- roc(rt$status, rt$riskScore)
plot(rocobj, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
dev.off()
