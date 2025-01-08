library(hdrcde)
nhanes <- read.delim("nhanes-aw.txt")
nhanes <- cbind(nhanes, "Female")
colnames(nhanes)[3] <- "Gender"
nhanes_m <- read.delim("nhanes-am.txt")
nhanes_m <- cbind(nhanes_m, "Male")
colnames(nhanes_m) <- colnames(nhanes)
nhanes <- rbind(nhanes, nhanes_m)
mod <- ranger::ranger(Height ~ Weight + Gender, data = nhanes)
#summary(mod)http://www.reddit.com/r/all/
resids <- mod$predictions - nhanes$Height
plot(resids ~ nhanes$Weight)
plot(resids[nhanes$Weight < 280] ~ nhanes$Weight[nhanes$Weight < 280])
score_fun <- function(y_hat, y_true){ #for upper one-sided, add for upper, subtract for lower
    (y_true - y_hat) 
}
quant_fun <- function(dat, alpha, n_size, direction){
    o <- order(dat)
    if(all.equal(alpha, 1) == TRUE) {quantile(dat, 1)}
    else if(all.equal(alpha, 0) == TRUE) {quantile(dat, 0)}
    else if(direction == "low"){dat[o][ceiling((alpha) * (n_size + 1))]}
    else if(direction != "low"){dat[o][ceiling((alpha) * (n_size + 1)) - 1]}
}

perms <- 2000
cov <- matrix(NA, nrow = perms, ncol = 1107) 
cov_m <- matrix(NA, nrow = perms, ncol = 1107) 
cov_f <- matrix(NA, nrow = perms, ncol = 1107) 
weight_vals <- matrix(NA, nrow = perms, ncol = 1107)
len <- matrix(NA, nrow = perms, ncol = 1107) 
cov_parm <- matrix(NA, nrow = perms, ncol = 1107)
len_parm <- matrix(NA, nrow = perms, ncol = 1107) 
cov_parm_m <- matrix(NA, nrow = perms, ncol = 1107) 
cov_parm_f <- matrix(NA, nrow = perms, ncol = 1107) 
set.seed(1)
for(i in seq(perms)){
    indicies <- sample(5107, 5107, FALSE)
    dat <- nhanes[indicies, ]
    dat_train <- dat[1:2000, ] #500 for mean model
    dat_cal <- dat[2001:4000, ] #1000 calibration
    dat_out <- dat[4001:5107, ] #out of sample test
    
    mod_train <- lm(Height ~ Weight + Gender, data = dat_train)
    
    pred_cal <- predict(mod_train, newdata = dat_cal)
    
    upper_scores <- score_fun(pred_cal, dat_cal$Height)
    lower_scores <- -1 * upper_scores
    
    
    dens <- density(lower_scores, kernel = "gaussian", from = range(lower_scores)[1], to = range(lower_scores)[2], n = 500)
    hpd_result <- hdr(lower_scores, prob = 90, den = dens)
    
    cdf <- ecdf(lower_scores)
    if(length(hpd_result$hdr) > 4){
        
        lower_alpha1 <- cdf(max(hpd_result$hdr)) - cdf(hpd_result$hdr[3]) + cdf(hpd_result$hdr[2]) - 0.9
        lower_alpha2 <- cdf(hpd_result$hdr[3])
        
        upper_alpha1 <- cdf(hpd_result$hdr[2]) #because HPD was computed on lower_scores, need to flip it for upper scores
        upper_alpha2 <- cdf(max(hpd_result$hdr)) #because HPD was computed on lower_scores, need to flip it for upper scores
        
        upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, length(lower_scores), "upper")
        lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, length(lower_scores), "low")
        
        upper_error2 <- quant_fun(lower_scores, 1-upper_alpha2, length(lower_scores), "upper")
        lower_error2 <- quant_fun(lower_scores, 1-lower_alpha2, length(lower_scores), "low")
        
        pred_out <- predict(mod_train, newdata = dat_out)
        
        upper_out1 <- pred_out - upper_error1 
        lower_out1 <- pred_out - lower_error1 
        upper_out2 <- pred_out - upper_error2 
        lower_out2 <- pred_out - lower_error2 
        cov[i, ] <- (dat_out$Height > lower_out1 & dat_out$Height < upper_out1 | dat_out$Height > lower_out2 & dat_out$Height < upper_out2) 
        len[i, ] <- (upper_out1 - lower_out1 + upper_out2 - lower_out2)
    } else{
        lower_alpha1 <- cdf(max(hpd_result$hdr)) - 0.9
        
        upper_alpha1 <- cdf(max(hpd_result$hdr)) #because HPD was computed on lower_scores, need to flip it for upper scores
        
        upper_error1 <- quant_fun(lower_scores, 1-upper_alpha1, length(lower_scores), "upper")
        lower_error1 <- quant_fun(lower_scores, 1-lower_alpha1, length(lower_scores), "low")
        
        pred_out <- predict(mod_train, newdata = dat_out)
        
        upper_out1 <- (pred_out - upper_error1 )
        lower_out1 <- (pred_out - lower_error1 )
        temp <- mean(upper_out1 - lower_out1)
        cov[i, ] <- (dat_out$Height > lower_out1 & dat_out$Height < upper_out1) #quantile regression and upper/lower adjustments
        len[i, ] <- upper_out1 - lower_out1
    }
    
    lower_alpha1 <- cdf(hpd_result$hdr[1])
    
    upper_alpha1 <- 0.1 - lower_alpha1 #because HPD was computed on lower_scores, need to flip it for upper scores
    
    o <- order(upper_scores)
    upper_error1 <- upper_scores[o][ceiling((2000) * (1 - upper_alpha1))]
    o <- order(lower_scores)
    lower_error1 <- lower_scores[o][ceiling((2000) * (1 - lower_alpha1))]
    
    pred_out <- predict(mod_train, newdata = dat_out)
    
    upper_out1 <- (pred_out + upper_error1)
    lower_out1 <- (pred_out - lower_error1)
    cov[i, ] <- dat_out$Height > lower_out1 & dat_out$Height < upper_out1 #quantile regression and upper/lower adjustments
    len[i] <- mean(upper_out1 - lower_out1)
    
    param_dat <- rbind(dat_train, dat_cal)
    linear_mod <- lm(Height ~ poly(Weight, 1) + Gender, data = param_dat)
    
    param_interval <- predict(linear_mod, newdata = dat_out, interval = "predict", level = 0.90)
    cov_parm[i, ] <- dat_out$Height > param_interval[, 2] & dat_out$Height < param_interval[, 3] #quantile regression and upper/lower adjustments
    len_parm[i, ] <- (param_interval[, 3] - param_interval[, 2])
    
    weight_vals[i, ] <- dat_out$Weight
    
    cov_m[i, ] <- c(cov[i, dat_out$Gender == "Male"], rep(NA, 1107 - sum(dat_out$Gender == "Male")))
    cov_f[i, ] <- c(cov[i, dat_out$Gender == "Female"], rep(NA, 1107 - sum(dat_out$Gender == "Female")))
    
    cov_parm_m[i, ] <- c(cov_parm[i, dat_out$Gender == "Male"], rep(NA, 1107 - sum(dat_out$Gender == "Male")))
    cov_parm_f[i, ] <- c(cov_parm[i, dat_out$Gender == "Female"], rep(NA, 1107 - sum(dat_out$Gender == "Female")))
    
    
}
mean(cov)
mean(len)
mean(cov_parm)
mean(len_parm)
mean(cov_m, na.rm = TRUE)
mean(cov_f, na.rm = TRUE)
mean(cov_parm_m, na.rm = TRUE)
mean(cov_parm_f, na.rm = TRUE)

sd(rowMeans(cov)) / sqrt(perms)
sd(rowMeans(len)) / sqrt(perms)
sd(rowMeans(cov_m, na.rm = TRUE), na.rm = TRUE) / sqrt(perms)
sd(rowMeans(cov_f, na.rm = TRUE), na.rm = TRUE) / sqrt(perms)

sd(rowMeans(cov_parm)) / sqrt(perms)
sd(rowMeans(len_parm)) / sqrt(perms)
sd(rowMeans(cov_parm_m, na.rm = TRUE), na.rm = TRUE) / sqrt(perms)
sd(rowMeans(cov_parm_f, na.rm = TRUE), na.rm = TRUE) / sqrt(perms)

# mean(cov[dat_out$Weight > median(dat_out$Weight)])
# mean(cov[dat_out$Weight < median(dat_out$Weight)])
# 
# 
# mean(cov_parm[dat_out$Weight < median(dat_out$Weight)])
# mean(cov_parm[dat_out$Weight > median(dat_out$Weight)])
# mean(cov_parm[dat_out$Gender == "Male"])
# mean(cov_parm[dat_out$Gender == "Female"])
## plot KDE-HPD
library(ggplot2)
x_out <- dat_out$Weight
y_plot_low1 <- lower_out1
y_plot_high1 <- upper_out1
plot_dat <- data.frame(x = nhanes$Weight, y = nhanes$Height)
plot_dat2 <- data.frame(ylow1 = y_plot_low1, yhigh1 = y_plot_high1, x = x_out, y = dat_out$Height)
ggplot(data = plot_dat2, aes(x = x)) + 
    scale_y_continuous(limits = c(50, 85)) +
    geom_ribbon(aes(ymin = ylow1, ymax = yhigh1), alpha = 1, fill = "orange") + 
    geom_point(data = plot_dat2, aes(x = x, y = y), shape = 21, alpha = 0.5, fill = "skyblue") +
    theme_classic() 

## plot parametric
x_out <- dat_out$Weight
y_plot_low1 <- param_interval[, 2]
y_plot_high1 <- param_interval[, 3]
plot_dat <- data.frame(x = nhanes$Weight, y = nhanes$Height)
plot_dat2 <- data.frame(ylow1 = y_plot_low1, yhigh1 = y_plot_high1, x = x_out, y = dat_out$Height)
ggplot(data = plot_dat2, aes(x = x)) + 
    scale_y_continuous(limits = c(50, 85)) +
    geom_ribbon(aes(ymin = ylow1, ymax = yhigh1), alpha = 1, fill = "orange") + 
    geom_point(data = plot_dat2, aes(x = x, y = y), shape = 21, alpha = 0.5, fill = "skyblue") +
    theme_classic() 