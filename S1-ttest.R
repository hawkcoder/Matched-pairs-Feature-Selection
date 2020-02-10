# ****************************************************
# Email: lsen@infervision.com
# Data: Wed Dec 20 09:58:32 2017
# Description: T-Test
#                
# *****************************************************


# *************************
# The Base t.test function
# *************************
ttest.modified <-function(x, y = NULL, paired = FALSE, method = c("Tan", "Cao", "ttest"), s0 = 1, 
                          alternative = c("two.sided", "less", "greater"),mu = 0,var.equal = FALSE, conf.level = 0.95, ...)
{
    alternative <- match.arg(alternative)
    
    if(!missing(mu) && (length(mu) != 1 || is.na(mu)))
      stop("'mu' must be a single number")
    
    if(!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || conf.level < 0 || conf.level > 1))
      stop("'conf.level' must be a single number between 0 and 1")
    
    if( !is.null(y) ) 
    {
      dname <- paste(deparse(substitute(x)),"and",deparse(substitute(y)))
      if(paired)
        xok <- yok <- complete.cases(x,y)  
      else {
        yok <- !is.na(y)
        xok <- !is.na(x)
      }
      y <- y[yok]
    } else {
      dname <- deparse(substitute(x))
      if (paired) stop("'y' is missing for paired test")
      xok <- !is.na(x)
      yok <- NULL
    }
    x <- x[xok]
    
    # 如果是paired数据， y变成了-> NULL
    if (paired) {
      
      # 修改处，加上对方法的判断
      if(method %in% c("Cao")){
        fc = as.numeric((x+0.1)/(y+0.1))
        x <- ifelse(fc >= 1, fc - 1, 1 - (1/fc))
        y <- NULL
      }else{
        x <- x-y
        y <- NULL
      }
      
    }
    
    # 计算长度，平均值， 方差
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    
    # 若是One Sample t-test, 或者 Paired t-test
    if(is.null(y)) {
      if(nx < 2) stop("not enough 'x' observations")
      df <- nx-1             # 自由度
      
      # 修改处，加上对方法的判断
      if(method %in% c("Tan", "Cao")){
        stderr <- sqrt(vx/nx) + s0
      }else{
        stderr <- sqrt(vx/nx)
      }
      
      if(stderr < 10 *.Machine$double.eps * abs(mx))
        stop("data are essentially constant")
      
      tstat <- (mx-mu)/stderr 
      
      
      
      method <- ifelse(paired, "Paired t-test", "One Sample t-test")
      estimate <- setNames(mx, ifelse(paired, "mean of the differences", "mean of x"))
      
    # 若是Wlech,或者 Two Sample t-test
    } else {
      ny <- length(y)
      if(nx < 1 || (!var.equal && nx < 2)) stop("not enough 'x' observations")
      if(ny < 1 || (!var.equal && ny < 2)) stop("not enough 'y' observations")
      if(var.equal && nx+ny < 3) stop("not enough observations")
      
      my <- mean(y)
      vy <- var(y)
      method <- paste(ifelse(!var.equal, "Welch", "Two Sample t-test"))
      estimate <- c(mx,my)
      names(estimate) <- c("mean of x","mean of y")
      
      
      if(var.equal) {
        df <- nx+ny-2
        v <- 0
        if(nx > 1) v <- v + (nx-1)*vx
        if(ny > 1) v <- v + (ny-1)*vy
        v <- v/df
        stderr <- sqrt(v*(1/nx+1/ny))
      } else {
        stderrx <- sqrt(vx/nx)
        stderry <- sqrt(vy/ny)
        stderr <- sqrt(stderrx^2 + stderry^2)
        df <- stderr^4/(stderrx^4/(nx-1) + stderry^4/(ny-1))
      }
      if(stderr < 10 *.Machine$double.eps * max(abs(mx), abs(my)))
        stop("data are essentially constant")
      tstat <- (mx - my - mu)/stderr
    }
    
    # 检测是单边(less, greater)还是双边(alternative)
    if (alternative == "less") {
      pval <- pt(tstat, df)
      cint <- c(-Inf, tstat + qt(conf.level, df) )
    } else if (alternative == "greater") {
      pval <- pt(tstat, df, lower.tail = FALSE)
      cint <- c(tstat - qt(conf.level, df), Inf)
    }else {
      pval <- 2 * pt(-abs(tstat), df)
      alpha <- 1 - conf.level
      cint <- qt(1 - alpha/2, df)
      cint <- tstat + c(-cint, cint)
    }
    
    
    cint <- mu + cint * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if(paired || !is.null(y)) "difference in means" else "mean"
    attr(cint,"conf.level") <- conf.level
    
    
    rval <- list(statistic = tstat, parameter = df, p.value = pval,
                 conf.int = cint, estimate = estimate, null.value = mu,
                 alternative = alternative,
                 method = method, data.name = dname)
    class(rval) <- "htest"
    return(rval)
  }




# *************************
# Test Statistic: (1) paired t-test
# *************************
test_pttest = function(data_case, data_ctrl)
{
  if(any(dim(data_case) != dim(data_ctrl))) {
    stop("Error: dim of data_case and data_ctrl are not equal!")
  }
  
  # build the resutl dataframe
  df_result = data.frame(features= rownames(data_case), pvalue = NA)
  
  for (i in 1:nrow(data_case)) {
    df_result[i,"pvalue"] = ttest.modified(as.numeric(data_case[i,]), as.numeric(data_ctrl[i,]), paired = T, method = "ttest")$p.value
  }
  
  df_result = df_result[order(df_result$pvalue),]
  
  return(df_result)
  
}

# *************************
# Test Statistic: (2) modified paired t-test
# *************************

test_mpttest = function(data_case, data_ctrl)
{
  if(any(dim(data_case) != dim(data_ctrl))) {
    stop("Error: dim of data_case and data_ctrl are not equal!")
  }
  
  
  # Step1: get S0, we set to the median of Sk
  data_case = as.matrix(data_case)
  data_ctrl = as.matrix(data_ctrl)
  data_diff = data_case - data_ctrl          # the differial of case and contrl
  
  d_hat = rowMeans(data_diff)                # the mean difference of `data_diff` 
  n = ncol(data_case)                        # the number of cases 
  k = nrow(data_case)                        # the number of features
  
  s = as.numeric(sqrt(rowSums((data_diff - d_hat)^2) / (n-1)))    # the variacnce of each feature
  s0 = mean(s)                               # we need this positive constant 
  
  # Step2: do this modified t.test function
  df_result = data.frame(features= rownames(data_case), pvalue = NA)
  for (i in 1:nrow(data_case)) {
    df_result[i,"pvalue"] = ttest.modified(as.numeric(data_case[i,]), as.numeric(data_ctrl[i,]), paired = T, method = "Tan", s0 = s0)$p.value
  }
  
  df_result = df_result[order(df_result$pvalue),]
  
  return(df_result)
}


# *************************
# Test Statistic: (3) fold-change paired t-test
# *************************
test_fcpttest = function(data_case, data_ctrl)
{
  
  if(any(dim(data_case) != dim(data_ctrl))) {
    stop("Error: dim of data_case and data_ctrl are not equal!")
  }
  
  # Step1: get S0, we set to the median of Sk
  data_case = as.matrix(data_case)
  data_ctrl = as.matrix(data_ctrl)
  data_diff = data_case - data_ctrl          # the differial of case and contrl
  
  d_hat = rowMeans(data_diff)                # the mean difference of `data_diff` 
  n = ncol(data_case)                        # the number of cases 
  k = nrow(data_case)                        # the number of features
  
  s = as.numeric(sqrt(rowSums((data_diff - d_hat)^2) / (n-1)))    # the variacnce of each feature
  s0 = mean(s)                               # we need this positive constant 
  
  # Step2: do this modified t.test function
  df_result = data.frame(features= rownames(data_case), pvalue = NA)
  for (i in 1:nrow(data_case)) {
    df_result[i,"pvalue"] = ttest.modified(as.numeric(data_case[i,]), as.numeric(data_ctrl[i,]), paired = T, method = "Cao", s0 = s0)$p.value
  }
  
  df_result = df_result[order(df_result$pvalue),]
  
  return(df_result)
}





