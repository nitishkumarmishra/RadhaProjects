library(ROCR)


aucs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "auc")@y.values)
rocs <- lapply(dat, function(x) performance(prediction(x[, "value"], x[, "group"]), "tpr", "fpr"))

# Create directory (if it doesn't exist)
if (!dir.exists("roc_plots_18Sept2017")) dir.create("roc_plots_18Sept2017")   ### change directory name here

# Loop through ROCs to make plots.
ns <- unname(table(dat[[1]]$group))
auc_df <- data.frame(name = character(length(dat)), 
                     AUC = numeric(length(dat)), 
                     CI_lower = numeric(length(dat)), 
                     CI_upper = numeric(length(dat)), stringsAsFactors = F)

for (i in 1:length(rocs)) {
  roc <- rocs[[i]]
  auc <- aucs[[i]]
  name <- names(aucs[i])
  
  # Flip results if AUC < 0.5
  if (auc < 0.5) {
    auc <- 1 - unlist(auc)
    temp <- roc@x.values
    roc@x.values <- roc@y.values
    roc@y.values <- temp
  }
  
  ci <- ci_auc(unlist(auc), ns[1], ns[2])
  ci_chr <- paste0(round(ci, 2), collapse = ", ") 
  
  auc_df[i, "name"] <- name
  auc_df[i, c("AUC", "CI_lower", "CI_upper")] <- c(auc, ci[1], ci[2])
  
  if (auc > 0.5) {
    file_name = paste0("roc_plots_18Sept2017/roc_", name, ".png")    ### change directory name here (the one you created above.
    png(file_name)
    #plot(roc, main = name, col = "red")
    plot(roc, main = name, colorize=T)  
    abline(0, 1, lty = 2)
    text(0.65, 0.20, paste0("AUC = ", round(unlist(auc), 2)), cex = 1.5)
    text(0.65, 0.15, paste0("95% CI [", ci_chr, "]"), cex = 1.2)
    rect(0.47, 0.10, 0.82, 0.25)
    
    dev.off()
  }
}

auc_df <- arrange(auc_df, desc(AUC)) %>%
  filter(AUC > 0.5)
write.csv(auc_df, "auc_table_18Sept2017.csv", row.names = F)   ### change filename to save AUC values.