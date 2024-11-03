library(rms)
library(PredictABEL)
library(ggDCA)
library(rmda)

rm(list = ls())
gc()

#load('./Rdata_08_ROC_analysis.RData')
setwd('/data/home/jiangminghui/KC/Project-0035/09_Nomogram_construction/')

GSE77938_exp <-readRDS('../00_raw_data/combat_data.rds')
Candidate_biomarker <-readRDS('../07_ValidateExpressions/01_Candidate_biomarker.rds')

data <- GSE77938_exp[Candidate_biomarker,] %>%t() %>%as.data.frame()
data$Group <- ifelse(grepl("^KC", colnames(GSE77938_exp)), 1,0)
saveRDS(data,'./01_data.rds')
write.csv(data ,'./01_data.csv',row.names = T)

dd <- datadist(data)
options(datadist = 'dd')

cat(Candidate_biomarker)
# 构建逻辑回归模型
model <- lrm(Group ~ ADAM8 +FAM135A +LCORL+ PINK1+ SFI1, data = data)

# 绘制列线图
nomogram <- nomogram(model, fun = plogis, fun.at = c(0.1, 0.5, 0.9), lp = FALSE)

pdf('./01_Nomogram_model.pdf',width = 10,height = 5)
plot(nomogram)
dev.off()

png('./01_Nomogram_model.png',width = 10,height = 5,units = 'in',res = 300)
plot(nomogram)
dev.off()

# predicted_probabilities <- predict(model, type = "fitted")
# 
# # Create a data frame with the actual and predicted probabilities
# calibration_data <- data.frame(
#   actual = data$Group,
#   predicted = predicted_probabilities
# )
# 
# pdf('./01_Calibration_curve.pdf', width = 10, height = 5)
# plotCalibration(data = calibration_data, cOutcome = 1, predRisk = calibration_data$predicted)
# dev.off()
# 
# png('./01_Calibration_curve.png', width = 10, height = 5, units = 'in', res = 300)
# plotCalibration(data = calibration_data, cOutcome = 1, predRisk = calibration_data$predicted)
# dev.off()
# 

fit1 <- lrm(Group ~ GADD45B + GATA6 + RBM24, data = data,x = T,y = T)
cali <-calibrate(fit1,B =1000)

pdf('./02_Calibration_Plot.pdf',width = 5,height = 5)
plot(cali,
     xlim = c(0,1),
     xlab = "Nomogram Predicted Probability",
     ylab = "Actual Probability",
     legend = FALSE,
     subtitles = FALSE)

abline(0,1,col = "#04198C",lty = 2,lwd = 2) 

lines(cali[,c("predy","calibrated.orig")], 
      type = "l",lwd = 2,col="#D91E1E",pch =16)

lines(cali[,c("predy","calibrated.corrected")],
      type = "l",lwd = 2,col="#4E7300",pch =16)

legend(0.55,0.35,  ##legend是绘制图例的函数，有兴趣的可以去深入了解这个函数
       c("Apparent","Ideal","Bias-corrected"), #表示曲线名称的集合
       lty = c(2,1,1), #lty表示线条类型，1代表实线，2至6都是虚线，虚的程度不一样
       lwd = c(2,1,1), #lwd表示宽度，以默认值的相对大小来表示宽度
       col = c("#04198C","#D91E1E","#4E7300"),  #给曲线添加颜色，对应上面c("Ax","Ix","Bx")
       bty = "n") # bty为o表示加边框，，注意是字母，不是数字0。bty可以取6种字符，分别为“o”、“l”、“7”、“c”、“u”、“]”

dev.off()

png('./02_Calibration_Plot.png',width = 5,height = 5,units = 'in',res = 300)

plot(cali,
     xlim = c(0,1),
     xlab = "Nomogram Predicted Probability",
     ylab = "Actual Probability",
     legend = FALSE,
     subtitles = FALSE)

abline(0,1,col = "#04198C",lty = 2,lwd = 2) 

lines(cali[,c("predy","calibrated.orig")], 
      type = "l",lwd = 2,col="#D91E1E",pch =16)

lines(cali[,c("predy","calibrated.corrected")],
      type = "l",lwd = 2,col="#4E7300",pch =16)

legend(0.55,0.35,  ##legend是绘制图例的函数，有兴趣的可以去深入了解这个函数
       c("Apparent","Ideal","Bias-corrected"), #表示曲线名称的集合
       lty = c(2,1,1), #lty表示线条类型，1代表实线，2至6都是虚线，虚的程度不一样
       lwd = c(2,1,1), #lwd表示宽度，以默认值的相对大小来表示宽度
       col = c("#04198C","#D91E1E","#4E7300"),  #给曲线添加颜色，对应上面c("Ax","Ix","Bx")
       bty = "n") # bty为o表示加边框，，注意是字母，不是数字0。bty可以取6种字符，分别为“o”、“l”、“7”、“c”、“u”、“]”
dev.off()


form.bestglm<-as.formula(Group ~GADD45B + GATA6 + RBM24)


DCA.1<- decision_curve(formula=form.bestglm,
                       family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals =0.95,
                       study.design = 'cohort',
                       data = data)
saveRDS(DCA.1,'./03_DCA_data.rds')
DCA.1$derived.data
head(DCA.1$derived.data[,c("thresholds","NB","sNB","cost.benefit.ratio")])

pdf('./03_DCA_curve.pdf',width = 5,height = 6)
plot_decision_curve(DCA.1,
                    curve.names= c("Model Bestglm"),
                    xlim=c(0,0.8),
                    cost.benefit.axis =TRUE,
                    col = "#E64B35B2",
                    confidence.intervals =FALSE,
                    standardize = FALSE)
dev.off()

png('./03_DCA_curve.png',width = 5,height = 6,units = 'in',res = 300)
plot_decision_curve(DCA.1,
                    curve.names= c("Model Bestglm"),
                    xlim=c(0,0.8),
                    cost.benefit.axis =TRUE,
                    col = "#E64B35B2",
                    confidence.intervals =FALSE,
                    standardize = FALSE)
dev.off()

# plot_clinical_impact(DCA.1,
#                      population.size = 1000,
#                      cost.benefit.axis = TRUE,
#                      n.cost.benefits= 8,
#                      col =c("#E64B35B2","#4DBBD5B2"),
#                      confidence.intervals= FALSE)



# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", "Rdata_",dir_name, ".RData"))

