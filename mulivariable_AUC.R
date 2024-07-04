library(ROCR); library(MASS)
setwd("F:/OneDrive - University of Nebraska Medical Center/Radha/HS-marker-data-Heat Map/HS Network/Latest/New Data")
data <- read.table("2-HS Individual-multivariate linear regression analysis-Aug 3 2018-To Nitish.txt", header = TRUE, row.names = 2)
data[,c("AUC","Index", "CI_lower", "CI_upper")] <- NULL
colnames(data) <- gsub("X", "", colnames(data))
colnames(data) <- gsub(".AVG_Beta", "", colnames(data))
data1 <- as.data.frame(t(data))
data1$Status <- ifelse(grepl("HS", rownames(data1)), 1, 0)

#model.all <- glm(Status ~ cg15475055+cg10455757+cg11332822+cg05861071+cg08235057+cg27261837+cg22081084+cg12875892+cg00419759+cg07368817+cg15731816+cg06707168+cg09028157+cg15358743+cg17420288+cg21308111+cg24782169+cg26112564+cg01152906+cg04317977+cg01366085+cg14875093+cg12506373+cg06745235+cg22430734+cg11087358+cg06509239+cg04436474+cg05165044+cg25642048+cg24891034+cg09273779+cg16151558+cg18125479+cg06084907+cg27495336+cg08204369+cg14022721+cg04579211+cg21774377+cg01311654+cg19872996+cg23104823+cg00302479+cg14522731+cg00717218+cg19191984+cg02209504+cg25390230+cg01353569+cg12344490+cg12393716+cg18880390+cg02773019+cg03705947+cg10671167+cg06074314+cg02908587+cg00904138+cg07935924+cg07869135+cg19545348+cg24738006+cg25929399+cg06310300+cg12220753+cg10659575+cg02529035+cg02756545+cg16672449+cg04691540+cg15877295+cg20284440+cg17103929+cg16265552+cg18652941+cg19741408+cg27507545+cg02374837+cg21864845+cg02002494+cg25631650+cg10918202+cg03226245+cg27556566+cg04482597+cg24332422+cg17434043+cg09488452+cg26283532+cg10483776+cg18986451+cg06124975+cg24243429+cg26571739+cg25494151+cg21515243+cg26512142+cg02664216+cg19839588+cg11502736+cg11602598+cg09311052+cg13749939+cg01090089+cg16104446+cg11226328+cg05728519+cg27566858+cg05743054+cg02760553+cg26221971+cg24420395+cg12820608+cg05485520+cg07409680+cg08688335+cg08052428+cg07251136+cg13338022+cg17298666+cg09844983+cg03734322+cg19379625+cg20623943+cg23303782+cg14534967+cg10089324+cg10740965+cg07044724+cg11490751+cg14147105+cg13030908+cg09137423+cg03569748+cg25307371+cg06590470+cg03067068+cg04969808+cg23760189+cg23337289+cg22975147+cg16109953+cg09564133+cg02396668+cg24023618+cg15008879+cg00664842+cg18584561+cg11530914+cg13412395+cg11922371+cg21961707+cg19853194+cg10397774+cg04077077+cg09455096+cg05409038+cg03270376+cg16627807+cg06408034+cg09379601+cg11522767+cg23182539+cg11342046+cg21871735+cg03945577+cg06665941+cg08655493+cg03160445+cg08664652+cg02576092+cg03311606+cg05491852+cg25947998+cg04100724+cg13038030+cg03076246+cg00734800+cg09820084+cg23171203+cg27304813+cg20538819+cg04238871+cg27539627+cg07349045+cg17524854+cg10436540+cg26519745+cg16823105+cg11050527+cg02602480+cg01302669+cg05360402+cg03672602+cg17741818+cg11399891+cg05379127+cg03984733+cg25846339+cg08862778+cg00894870+cg10083593+cg01652215+cg25578064+cg08667777+cg01337508+cg25441200+cg10869376+cg25683012+cg15309264+cg18788872+cg07758528+cg09487931+cg21813546+cg05053688+cg27570636+cg26164735+cg03157738+cg16065538+cg15459742+cg19432434+cg11355319+cg03759077+cg19699211+cg21517261+cg00426056+cg25386820+cg15042747+cg06419771+cg26125131+cg19383577+cg02098752+cg16857181+cg05104993+cg23260799+cg08099293+cg09060823+cg01166180+cg03640568+cg16139199+cg05477027+cg04122601+cg03115081, data=data1)
model.all <- glm(Status ~ ., data = data1)
step <- stepAIC(model.all, direction = "both")
step$anova ## Select the list of features based on The Akaike information criterion (AIC)
model.stepAIC <- glm(Status ~ cg15475055 + cg10455757 + cg11332822 + cg05861071 + 
                       cg08235057 + cg27261837 + cg22081084 + cg12875892 + cg00419759 + 
                       cg07368817 + cg15731816 + cg06707168 + cg09028157 + cg15358743 + 
                       cg17420288 + cg21308111 + cg24782169 + cg26112564 + cg01152906 + 
                       cg04317977 + cg01366085 + cg14875093 + cg12506373 + cg06745235 + 
                       cg22430734 + cg11087358 + cg06509239 + cg04436474 + cg05165044 + 
                       cg25642048 + cg24891034 + cg09273779 + cg16151558 + cg18125479 + 
                       cg06084907 + cg27495336 + cg08204369 + cg14022721 + cg04579211 + 
                       cg21774377 + cg01311654 + cg19872996 + cg23104823 + cg00302479 + 
                       cg14522731 + cg00717218 + cg19191984, data=data1)

se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

rocs <- performance(prediction(fitted(model.stepAIC), data1$Status), "tpr", "fpr")
aucs <- performance(prediction(fitted(model.stepAIC), data1$Status), "auc")
roc <- rocs
auc <- aucs@y.values[[1]]
#attributes(roc)$alpha.values[[1]][1] <- 1
attributes(roc)$alpha.values[[1]][1] <- 1 ## Alpha value >1 then replace it with 1. Here first five are more than 1 
#attributes(roc)$alpha.values[[1]][21:25] <- 0 ## Replace alpha vale < 1 with zero, here last five have bnegative values
ci <- ci_auc(auc, 12, 12) # Here we have 12 affcted and 12 control samples
ci_chr <- paste0(round(ci, 2), collapse = ", ") 
auc_df <- data.frame(name = c("Mulivariate"),
                     AUC = numeric(length(1)),
                     CI_lower = numeric(length(1)),
                     CI_upper = numeric(length(1)), stringsAsFactors = F)
auc_df[1, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
auc_df

save_roc_plot <- function(roc_obj) {}
## PDF image quality is better, so I will save figure in PDF
#file_name = paste0("ROC_mulivariate", "CpG", ".png")    ### change directory name here (the one you created above.
#png(filename = file_name, width = 1024, height = 768)
file_name = paste0("ROC_mulivariate_stepAIC", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()

save.image("Multiple_Linear_analysis.RData")

#### Without feature selection ########
model.all <- glm(Status ~ ., data=data1)
#step <- stepAIC(model.all, direction = "both")
#step$anova ## Select the list of features based on The Akaike information criterion (AIC)
#model.stepAIC <- glm(Status ~ cg10979903 + cg26814276 + cg18073151 + cg11125369 + cg25999722 + cg22943986 + cg12853563 + cg25729826 + cg21241839, data=data1)

se_auc <- function(auc, n1, n2) { #n1 is positive, n2 is negative
  q1 <- auc / (2 - auc)
  q2 <- (2 * auc ^ 2) / (1 + auc)
  sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc ^ 2) + (n2 - 1) * (q2 - auc ^ 2)) / (n1 * n2))
}

ci_auc <- function(auc, n1, n2, level = 0.95) {
  ci <- auc + c(-1, 1) * qnorm((1 + level) / 2) * se_auc(auc, n1, n2)
  ci[1] <- ifelse(ci[1] < 0, 0, ci[1])
  ci[2] <- ifelse(ci[2] > 1, 1, ci[2])
  ci }
save_roc_plot <- function(roc_obj) {}

rocs <- performance(prediction(fitted(model.all), data1$Status), "tpr", "fpr")
aucs <- performance(prediction(fitted(model.all), data1$Status), "auc")
roc <- rocs
auc <- aucs@y.values[[1]]
#attributes(roc)$alpha.values[[1]][1] <- 1
#attributes(roc)$alpha.values[[1]][1:5] <- 1 ## Alpha value >1 then replace it with 1. Here first five are more than 1 
#attributes(roc)$alpha.values[[1]][21:25] <- 0 ## Replace alpha vale < 1 with zero, here last five have bnegative values
ci <- ci_auc(auc, 12, 12) # Here we have 12 affcted and 12 control samples
ci_chr <- paste0(round(ci, 2), collapse = ", ") 
auc_df <- data.frame(name = c("Mulivariate"),
                     AUC = numeric(length(1)),
                     CI_lower = numeric(length(1)),
                     CI_upper = numeric(length(1)), stringsAsFactors = F)
auc_df[1, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
auc_df

save_roc_plot <- function(roc_obj) {}
## PDF image quality is better, so I will save figure in PDF
#file_name = paste0("ROC_mulivariate", "CpG", ".png")    ### change directory name here (the one you created above.
#png(filename = file_name, width = 1024, height = 768)
file_name = paste0("ROC_mulivariate All", "CpG", ".pdf")    ### change directory name here (the one you created above.
pdf(file_name)
#plot(roc, main = name, col = "red")
plot(roc, main = "Multiple linear regression", colorize=T)  
abline(0, 1, lty = 2)
### Adjust position based on figure (pdf or png)
text(0.80, 0.13, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.4)
text(0.80, 0.08, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
rect(0.65, 0.05, 0.95, 0.16)
dev.off()

##############################################
save.image("Multivariate_ROC.RData")
