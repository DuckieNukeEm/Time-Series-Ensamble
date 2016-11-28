
ensablem_ts_forecast = function (x,
                                 start = c(1900, 1), 
                                 frequency = 26, #frequency of the data set
                                 h = 26, #how many time steps forward we want to use
                                 models = c("HoltWinters", "STL", "ARIMA", "TBATS","Naive"), #the models
                                 zero_na = T, #replace NA'/ Inf with 0's
                                 rm_na = F, #should we remove NA's
                                 CI_type = 'mean', #type used to calculate the CI
                                 use = c("Y","Y"), #Do we want to keep the historical or future parts
                                 scale = 'F', #do we want to scale the data?
                                 log_transform = NULL, #if you want to log transform the data first,
                                 models_to_calc_with = models, #what models do you want to use in the calculations?
                                 growth_limit = 0.01, #what is the avearge growth limit of a model before we cut it from consideraionts?
                                 blend_threshhold = NA, #what is our percentage change we want to use for the blend (NA menas do not blend)
                                 to_index = NA,
                                 ma_adj = NA,
                                 ...
)
{
  #checking the inout data
  if(!is.vector(x) | !is.numeric(x)){stop("x must be a vector")}
  if(NROW(x) < frequency * 2) {stop("x must be at least twice as long as the frequency")}
  
  
  #This df will hold the data for historical records
  df = data.frame(x)
  #if log transform the data is true, then log transform. note I'm not forcing replacing NA values with a number...so yea
  if(!is.null(log_transform))
  {
    df[,1] = log(df[,1])
  }
  #removing NA's or replacing the NA's with zero, if removing NA's, then replace is truned off  
  if(rm_na & zero_na){ zero_na = F}
  
  if (zero_na) {
    df[is.na(df[,1]) | is.infinite(df[,1]) | is.nan(df[,1]),1] = 0
  }
  if(rm_na){
    
    df_rows_rm = which(is.na(df$x))
    df = data.frame( x = df[!is.na(df$x),])
    if(length(df[,1]) < (frequency * 2)) 
    {
      frequency = as.integer(nrow(df$x)/2)
    }
    
  }
  
  # Creating a list of Data frames, each column of the dataframe will hold the forecasted values
  # and chec dataframe is unique to EACH forecasting method
    #now building the proper dataframes
     ft_mean = lapply(models, function(x) data.frame(matrix(0,nrow = h, ncol = length(to_index))))
       ######################
       # WORKING HERE YA BASTARD!!!!!!
       ######################
    
    names(ft_mean) = models
    ft_upper = ft_mean
    ft_lower = ft_mean
  #this will hold each models forecasted values
  ft = data.frame(rep(0, h))
  #error flag incase any of the methods fail
  error_flag = FALSE
  
  names(ft) = names(df)
  #running the models
  for (m in models)
  {
    for (i in to_index){
      h_adj =  (nrow(df) - to_index[i]) + h
      ret = models_runs(method = m, 
                        start= start, 
                        frequency = frequency, 
                        df = df[,1], 
                        h = h_adj, 
                        ma_adj = ma_adj)
      
    
        frcst = ret$forecast
        
        if(NROW(ret$df) == 1){
          models_to_calc_with = setdiff(models_to_calc_with, m)
          error_flag = TRUE
        } else { #Not an error
          if(to_index[i] == nrow(df)){
              df[,m] =ret$df
              }
        }
        
        
        #updating dataframes
        
        if(!error_flag){ #if there wasn't an error
          
            
            frcst_mean =  cbind(ft_mean, as.double(temp_f[(h_adj - h + 1):h,]$mean))
            frcst_upper = cbind(ft_upper, as.double(temp_f[(h_adj - h + 1):h,]$upper))
            frcst_lower = cbind(ft_lower, as.double(temp_f[(h_adj - h + 1):h,]$lower))
            
            
          } else {
            
            frcst_mean = as.double(frcst$mean)
            frcst_upper = as.numeric(frcst$upper[,2])
            frcst_lower = as.numeric(frcst$lower[,2])
          }
          
          
        # ft[,m] = frcst_mean
        # ft[,paste(c(m,"upper"), collapse = "_" )] =  frcst_upper
        # ft[,paste(c(m,"lower"), collapse = "_" )] =  frcst_lower
        # df[,paste(c(m,"upper"), collapse = "_" )] =  df[,m]
        # df[,paste(c(m,"lower"), collapse = "_" )] =  df[,m]
          
        } else { #IF THERE WAS AN ERROR
          ft[,m] = rep(0,h)
          ft[,paste(c(m,"upper"), collapse = "_" )] =  rep(NA,h)
          ft[,paste(c(m,"lower"), collapse = "_" )] =  rep(NA,h)
          df[,m] = rep(0,nrow(df))
          df[,paste(c(m,"upper"), collapse = "_" )] =  rep(NA,nrow(df))
          df[,paste(c(m,"lower"), collapse = "_" )] =  rep(NA,nrow(df))
        }
        
        #checking to make sure that the forecast didn't exploded to the moon or to the earths center
        if(growth_limit != 0.0 & !is.na(growth_limit)){
          explosion_check = try(abs(mean(c(diff(ft[,m]),0)/ft[,m])))
          if(inherits(explosion_check, 'try-error')) {
            models_to_calc_with = setdiff(models_to_calc_with, m)
            print(paste(c("Model ",m,"exceeded % change threeshold, removing model from consideration"), collapse = " "))
          } else if (explosion_check > growth_limit & !is.na(explosion_check)) {
            models_to_calc_with = setdiff(models_to_calc_with, m)
            print(paste(c("Model ",m,"exceeded % change threeshold, removing model from consideration"), collapse = " "))}
          
        }
        #if the model passed, we will add it to the forecast
        if(m %in% models_to_calc_with){
          ft[,1] = ft[,1] + ft[,m]
        }
        #reseting the error flag
        error_flag = FALSE
      }
  
      ft_mean = ft_mean[,2:ncol(ft_mean)]
      ft_upper = ft_mean[,2:ncol(ft_upper)]
      ft_lower = ft_mean[,2:ncol(ft_lower)]
  
      for (i in ncol(ft_mean):2 {
        if(ft_mean[2*frequency, i]/ft_mean[2*frequency, i-1] > blend_threshhold){
          ft_mean[, i] = 0.5*ft_mean[, i] + 0.5*ft_mean[, i-1]
          ft_upper[, i] = 0.5*ft_upper[, i] + 0.5*ft_upper[, i-1]
          ft_lower[, i] = 0.5*ft_lower[, i] + 0.5*ft_lower[, i-1]
        }
        
      }
  
  
      ft[,m] = frcst_mean
      ft[,paste(c(m,"upper"), collapse = "_" )] =  frcst_upper
      ft[,paste(c(m,"lower"), collapse = "_" )] =  frcst_lower
      df[,paste(c(m,"upper"), collapse = "_" )] =  df[,m]
      df[,paste(c(m,"lower"), collapse = "_" )] =  df[,m]
      
  
  
  }
  
  if(length(models_to_calc_with) == 0)
  {
    #just in case ALL the models fail, I can still generate a naive linear forecast
    LM = lm(df[,1]~index(df[,1]))
    print(LM)
    ft[,"HoltWinters"] = coef(LM)[2]*(max(index(df[,1]))+1):(max(index(df[,1]))+h) + coef(LM)[1]
    ft[,"HoltWinters_upper"] =  ft[,"HoltWinters"] + sd(LM$residuals)*1.645
    ft[,"HoltWinters_lower"] =  ft[,"HoltWinters"] - sd(LM$residuals)*1.645
    
    df[,"HoltWinters"] = as.numeric(LM$fitted.values)
    df[,"HoltWinters_upper"] =  df[,"HoltWinters"] 
    df[,"HoltWinters_lower"] =  df[,"HoltWinters"] 
    models_to_calc_with = "HoltWinters"
  }  
  
  #Averaging the values together
  ft[,1] = ft[,1]/length(models_to_calc_with)
  
  #backfilling the NA (if we opeted to remove them) with NA's
  if(rm_na){
    m = data.frame(matrix(0, nrow = length(df_rows_rm), ncol =ncol(df) ))
    names(m) = names(df)
    df = rbind(m,df)
  }
  
  names(df) = names(ft)
  
  #if we added a use vectory, appending it to the end of the data  
  if(length(use) == 2) {
    df$use = use[1]
    ft$use = use[2]
  } else if (length(use) == 1 & !is.na(use)){
    df$use = use
    ft$use = use
  } 
  
  if(scale) #if scale is turned on
  {
    ft$x = df_scale(as.double(ft$x), method = 'Hard')
  }
  
  df = rbind(df,ft)
  
  #gotta revert our numbers if we log transofrmed them first
  if(!is.null(log_transform)){
    df[,n_v(df)] = exp(df[,n_v(df)])
  }
  
  
  #adding in the CI now
  #we need to check if certain methods have been dropped
  t_to_use = c( TRUE , #x
                ifelse('HoltWinters' %in% models_to_calc_with, c(TRUE, TRUE, TRUE), c(FALSE,FALSE,FALSE)), #HoltWinter placeholder
                ifelse('STL' %in% models_to_calc_with, c(TRUE, TRUE, TRUE), c(FALSE,FALSE,FALSE)),
                ifelse('ARIMA' %in% models_to_calc_with, c(TRUE, TRUE, TRUE), c(FALSE,FALSE,FALSE)),
                ifelse('TBATS' %in% models_to_calc_with, c(TRUE, TRUE, TRUE), c(FALSE,FALSE,FALSE)),
                ifelse('Naive' %in% models_to_calc_with, c(TRUE, TRUE, TRUE), c(FALSE,FALSE,FALSE)),
                TRUE
  )
  df = ensamble_ts_CI(df, CI_type)
  
  
  return(df)
}



