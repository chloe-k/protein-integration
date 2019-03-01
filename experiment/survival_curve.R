library("survival")
library("survminer")
require("survival")

id <- "result18_3_GM"
profile_name <- c("rna(Entrez)", "meth(Entrez)")
id <- "result18_3_GMR"
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
id <- "result28_0.6_GR"
profile_name <- c("rna(Entrez)", "rppa(Entrez)")

#-------------- Regression
desc <- c(id, "DRW", "txt")

load(file.path(datapath, paste(c("pathway_rank", id, profile_name, "DRW", "RData"), collapse = '.')))
toppath <- file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.'))
df <- read.delim(file = toppath, header = T, sep = '')


print('Getting top N pathways is done..')

# Model evaluation with top N pathway
set.seed(111)
rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]

trControl <- trainControl(method = "LOOCV", savePredictions = TRUE)
result <- train(x = X[,rankn_feats], y = survival, method="rf", metric = "RMSE", 
                trControl=trControl, importance = TRUE)
gp <- result
save(result, file=file.path('data/model/res_pa_GP_surv.RData'))

#--------------- Plot

period <- c(0:14)
#gm <- get(load('data/model/res_pa_GM_surv.RData'))
gm$pred <- gm$pred[which(gm$pred$mtry == 2),]
#gmp <- get(load('data/model/res_pa_GMP_surv.RData'))
gmp$pred <- gmp$pred[which(gmp$pred$mtry == 2),]

survival <- round(clinical$survival/365)
gm_pred <- round(gm$pred$pred/365)
gmp_pred <- round(gmp$pred$pred/365)

res_models <- list(survival, gm_pred, gmp_pred)
title <- c("Survival Curve")
xlabs <- c("Real", "iDRW(GM)", "iDRW(GMP)")

#----------------Line Plot
#perf_lineplot_multi(title, xlabs, res_models, period)

#perf_lineplot_multi <- function(title, xlabs, res_models, period) {
  
  df_list = list()
  idx <- 1
  for(i in 1:length(res_models)) {
    pred_freq <- c(1:15)*0
    temp <- table(unlist(res_models[[i]]))
    for(k in 1:length(temp)){
      pred_freq[k] <- temp[k]
    }
    pred_ratio <- c(1:15)*0
    pred_ratio[1] <- 376
    for(k in 2:length(temp)){
      pred_ratio[k] <- (pred_ratio[k-1]-pred_freq[k-1])
    }
    for(j in 1:15){
      df_list[[idx]] = data.frame(Group=c(rep(xlabs[i],15)), Ratio=pred_ratio[j]/376,period=period[j])
      idx <- idx+1
    }
  }
  
  df = Reduce(rbind, df_list)
  p <- ggplot(df, aes(x=period, y=Ratio, group=Group, color=Group)) +
    geom_line() +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(name=c("Ratio of survival")) +
    scale_x_discrete(limits = c(0:15))
  print(p + ggtitle(title))
#}
  
#----------- Line plot for two groups
  gm_model <- get(load('data/model/res_pa_GM_18_3_LOOCV.RData'))
  gmp_model <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
  
  gm_model$pred <- gm_model$pred[which(gm_model$pred$mtry == 26),]
  gm_model$pred <- cbind(gm_model$pred, survival_days=round(survival/365))
  
  gmp_model$pred <- gmp_model$pred[which(gmp_model$pred$mtry == 36),]
  gmp_model$pred <- cbind(gmp_model$pred, survival_days=round(survival/365))
  
  gmp_s <- gmp_model$pred[which(gmp_model$pred$pred == 0),]
  gmp_s <- cbind(gmp_s, survival_days=round(survival[gmp_s$rowIndex]/365))
  gmp_l <- gmp_model$pred[which(gmp_model$pred$pred == 1),]
  gmp_l <- cbind(gmp_l, survival_days=round(survival[gmp_l$rowIndex]/365))
  
  title <- c("The ratio of survival in iDRW(GMP)")
  xlabs <- c("Short-term", "Long-term")
  res_models <- list(gmp_s, gmp_l)
  
  
  gm_s <- gm_model$pred[which(gm_model$pred$pred == 0),]
  gm_s <- cbind(gm_s, survival_days=round(survival[gm_s$rowIndex]/365))
  gm_l <- gm_model$pred[which(gm_model$pred$pred == 1),]
  gm_l <- cbind(gm_l, survival_days=round(survival[gm_l$rowIndex]/365))
  
  title <- c("The ratio of survival in iDRW(GM)")
  xlabs <- c("Short-term", "Long-term")
  res_models <- list(gm_s, gm_l)
  
  
  df_list = list()
  idx <- 1
  for(i in 1:length(res_models)) {
    pred_freq <- c(1:15)*0
    temp <- table(unlist(res_models[[i]]$survival_days))
    for(k in 1:length(temp)){
      pred_freq[k] <- temp[k]
    }
    pred_ratio <- c(1:15)*0
    pred_ratio[1] <- dim(res_models[[i]])[1]
    for(k in 2:length(temp)){
      pred_ratio[k] <- (pred_ratio[k-1]-pred_freq[k-1])
    }
    for(j in 1:15){
      df_list[[idx]] = data.frame(Group=c(rep(xlabs[i],15)), Ratio=pred_ratio[j]/376, period=period[j])
      idx <- idx+1
    }
  }
  
  df = Reduce(rbind, df_list)
  
  #----For KM graph and Log-rank test
  survival <- clinical$survival
  gm_model <- get(load('data/model/res_pa_GM_18_3_LOOCV.RData'))
  gmp_model <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
  
  gm_model$pred <- gm_model$pred[which(gm_model$pred$mtry == 26),]
  gm_model$pred <- cbind(gm_model$pred, survival_days=round(survival/365))
  
  gmp_model$pred <- gmp_model$pred[which(gmp_model$pred$mtry == 36),]
  gmp_model$pred <- cbind(gmp_model$pred, survival_days=round(survival/365))
  
  total <- rbind(cbind(gm_model$pred, model="GM"), cbind(gmp_model$pred, model="GMP"))
  #### iDRW(GM)
  # create life table
  gm_lt <- gm_model$pred[order(gm_model$pred$survival_days),]
  gm_lt <- gm_lt[c("pred", "obs", "survival_days")]
  #gm_lt <- gm_lt[c("pred", "survival_days")]
  # Plot KM
  gm_rx <- survfit(Surv(survival_days, pred==1)~obs, data = gm_lt)

  title <- c("Survival curve in iDRW(GM)")
  plot(gm_rx, xlab = "Survival Period (year)", ylab="Ratio of Survival",col=c("blue", "red"), lty = 1:2, main="Survival curve in iDRW(GM)")
  legend("topright",c("Short-term survival", "Long-term survival"), lty = 1:2, col=c("blue", "red"))
  #######2019/02/28 in final revision
  # combine two survival curves (GM)
  fit_gm <- survfit(Surv(survival_days, obs==0)~pred, data=gm_lt)
  gm_plt <- ggsurvplot(fit_gm, data=gm_lt,
                        conf.int.style = "step",
                        risk.table = TRUE,
                        palette = 'jco',
                        break.x.by = 3,
                        legend.labs=c("Short-term\nsurvival", "Long-term\nsurvival"))
  gm_plt+ggtitle("Survival curve in iDRW(GM)")+
    xlab("Survival period (year)")+ylab("Ratio of Survival")

    # Log-rank test
  surv_diff_gm <- survdiff(Surv(survival_days, pred == 1)~obs, data=gm_lt)
  
  #### iDRW(GMP)
  gmp_lt <- gmp_model$pred[order(gmp_model$pred$survival_days),]
  gmp_lt <- gmp_lt[c("pred", "obs", "survival_days")]
  #gmp_lt <- as.data.frame(gmp_lt[c("pred", "survival_days")])
  
  # Plot KM
  gmp_rx <- survfit(Surv(survival_days, pred==1)~obs, data = gmp_lt)
  
  title <- c("Survival curve in iDRW(GMP)")
  ggplot(gmp_rx, xlab = "Survival Period (year)", ylab="Ratio of Survival", col=c("blue", "red") ,lty = 1:2, main="Survival curve in iDRW(GMP)")
  legend("topright",c("Short-term survival", "Long-term survival"), lty = 1:2, col=c("blue", "red"))
  ###########2019/02/28 in final revision
  # combine two survival curves (GMP)
  fit_gmp <- survfit(Surv(survival_days, obs==0)~pred, data=gmp_lt)
  gmp_plt <- ggsurvplot(fit_gmp, data=gmp_lt,
             conf.int.style = "step",
             risk.table = TRUE,
             palette = 'jco',
             break.x.by = 3,
             legend.labs=c("Short-term\nsurvival", "Long-term\nsurvival"))
  gmp_plt+ggtitle("Survival curve in iDRW(GMP)")+
    xlab("Survival period (year)")+ylab("Ratio of Survival")
  # Log-rank test
  surv_diff_gmp <- survdiff(Surv(survival_days, pred == 1)~obs, data=gmp_lt)
  
  #### Total(Short_Long term)
  # create life table
  st <- total[which(total$pred == 0),]
  lt <- total[which(total$pred == 1),]
  
  st <- st[order(st$survival_days),]
  lt <- lt[order(lt$survival_days),]
  
  ### Short-term survival
  # Plot KM
  st_rx <- survfit(Surv(survival_days, obs == 0)~model, data = st)
  
  plot(st_rx, xlab = "Period", ylab="Survival",col=c("black", "red"), lty = 1:2, main="Survival curve in Short-term survival")
  legend("topright",c("iDRW(GM)", "iDRW(GMP)"), lty = 1:2, col=c("black", "red"))
  
  # Log-rank test
  diff_st <- survdiff(Surv(survival_days, obs == 0)~model, data=st)
  
  ### Long-term survival
  # Plot KM
  lt_rx <- survfit(Surv(survival_days, obs == 1)~model, data = lt)
  
  plot(lt_rx, xlab = "Period", ylab="Survival",col=c("black", "red"), lty = 1:2, main="Survival curve in Long-term survival")
  legend("topright",c("iDRW(GM)", "iDRW(GMP)"), lty = 1:2, col=c("black", "red"))
  
  # Log-rank test
  diff_lt <- survdiff(Surv(survival_days, obs == 1)~model, data=lt)
  
  
  #-----------
  p <- ggplot(df, aes(x=period, y=Ratio, group=Group, color=Group)) +
    geom_line() +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(name=c("Ratio of survival")) +
    scale_x_discrete(limits = c(0:15))
  print(p + ggtitle(title))
  
#-----------Scatter plot
  title <- c("Survival Regression")
  
  survival <- clinical$survival
  gm_pred <- gm$pred$pred
  gmp_pred <- gmp$pred$pred
  
  res_models <- list(survival, gm_pred, gmp_pred)
  xlabs <- c("Real", "iDRW(GM)", "iDRW(GMP)")
  df_list <- list()
  
  for(i in 1:length(res_models)){
    df_list[[i]] = data.frame(model=c(rep(xlabs[i],376)), Predict=res_models[[i]], Actual=survival)
  }
  df = Reduce(rbind, df_list)
  
  p <- ggplot(df, aes(x=Actual, y=Predict, shape=model, color=model)) +
    geom_point() +
    scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) +
    geom_smooth(method="lm") +
    theme(plot.title = element_text(hjust = 0.5)) 
    #scale_x_continuous(limits = c(0,15))
  print(p + ggtitle(title))
  
#-----------Scatter plot (facet grid)
  title <- c("Survival Regression")
  
  survival <- clinical$survival
  gm_pred <- gm$pred$pred
  gmp_pred <- gmp$pred$pred
  
  res_models <- list(gm_pred, gmp_pred)
  xlabs <- c("iDRW(GM)", "iDRW(GMP)")
  df_list <- list()
  
  for(i in 1:length(res_models)){
    df_list[[i]] = data.frame(model=c(rep(xlabs[i],376)), Predict=res_models[[i]], Actual=survival)
  }
  df = Reduce(rbind, df_list)
  
  p <- ggplot(df, aes(x=Actual, y=Predict, shape=model, color=model)) +
    geom_point() +
    #scale_color_manual(values=c('#E69F00', '#56B4E9')) +
    geom_smooth(method="lm") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(0,4000), name=c("Actual survival period (days)")) +
    scale_y_continuous(limits = c(0,4000), name=c("Predicted survival period (days)"))
  print(p + ggtitle(title) + facet_grid(model ~ .))
  
#-------- Log-rank test
#GM
gm_pre <- as.numeric(gm_model$pred$pred)
gm_surv <- as.numeric(gm_model$pred$survival_days)
gm_obs <- as.numeric(gm_model$pred$obs)

survdiff(Surv(gm_surv, gm_obs)~gm_pre)

#GMP
gmp_pre <- as.numeric(gmp_model$pred$pred)
gmp_surv <- as.numeric(gmp_model$pred$survival_days)
gmp_obs <- as.numeric(gmp_model$pred$obs)

survdiff(Surv(gmp_surv, gmp_obs)~gm_pre)
