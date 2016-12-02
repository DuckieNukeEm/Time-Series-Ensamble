









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
                                 to_index = c(0),
                                 ma_adj = 0,
                                 CI = 95,
                                 de_season = F,
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
    #inisilizing a null list
      ft_mean = vector()
      ft_upper = vector()
      ft_lower = vector()
    #now building the proper dataframes
    ft_mean = lapply(models, function(x) data.frame(matrix(0,nrow = (h + NROW(df)), ncol = length(to_index))))
    names(ft_mean) = models
    ft_upper = ft_mean
    ft_lower = ft_mean

  #error flag incase any of the methods fail
  error_flag = FALSE
  
  names(ft) = names(df)
  #running the models
  for (m in models)
  {
    
    temp_list = lapply(to_index, function(x) models_runs(method =  m, 
                                                         start = start,
                                                         frequency = frequency,
                                                         df = df[1:(NROW(df) - abs(x))],
                                                         h = h + abs(x),
                                                         surpress.error = F,
                                                         season = de_season,
                                                         CI = CI))
                                                         
      for(i in 1:length(temp_list)){
        ft_mean[[m]][,i] = temp_list[[i]][,1]
        ft_upper[[m]][,i] = temp_list[[i]][,2]
        ft_lower[[m]][,i] = temp_list[[i]][,3]
      }
   
    ###
    #Check to see if the forecast exploded (ie grew or shrank WAY to fast)
    ###
    
    if(growth_limit != 0.0 & !is.na(growth_limit)){
      check = apply(ft_mean[[m]], 2, explosion_check, threshhold = growth_limit)
      if(sum(check) < ncol(ft_mean[[m]])){
        print(paste(c("We've had some explosions, for model ", m," at position", to_index[check], "going to remove them"),collapse = " "))
        #zeroing out all fields that 'exploded'
        ft_mean[[m]][,check] = rep(NA,nrow(ft_mean[[m]]))
        ft_upper[[m]][,check] = rep(NA,nrow(ft_upper[[m]]))
        ft_lower[[m]][,check] = rep(NA,nrow(ft_lower[[m]]))
        }
      
      } 
    
    
    ###
    #now we are going to smooth the data (if we want to smooth it)
    ###
      if(!is.na(blend_threshhold)){
          ft_mean[[m]] = df_smooth(ft_mean[[m]], blend_threshhold = blend_threshhold, start_index = (NROW(df)+1))
          ft_lower[[m]] = df_smooth(ft_lower[[m]], blend_threshhold = blend_threshhold, start_index = (NROW(df)+1))
          ft_upper[[m]] = df_smooth(ft_upper[[m]], blend_threshhold = blend_threshhold, start_index = (NROW(df)+1))
      }
   
        
    ###
    # if we log transforemd the variables, we want to transform them back.
    ###
    if(log_transform){
      zero_check = apply(ft_mean[[m]], 2, mean)
      zero_check = !is.na(zero_check)
      
      ft_mean[[m]][,zero_check] = exp(ft_mean[[m]][,zero_check])
      ft_upper[[m]][,zero_check] = exp(ft_upper[[m]][,zero_check])
      ft_lower[[m]][,zero_check] = exp(ft_lower[[m]][,zero_check])
    }
    
    ###
    # Now aggregating to the 'final' column of each
    ###
      if(ncol(ft_mean[[m]]) > 1){ #there is more than one column
        ft_mean[[m]]$Final = rowMeans(ft_mean[[m]], na.rm = TRUE)
        ft_upper[[m]]$Final = rowMeans(ft_upper[[m]], na.rm = TRUE)
        ft_lower[[m]]$Final = rowMeans(ft_lower[[m]], na.rm = TRUE)
      } else { #there is only one columns (rowMeans doesn't work)
        ft_mean[[m]]$Final = ft_mean[[m]][,1]
        ft_upper[[m]]$Final = ft_upper[[m]][,1]
        ft_lower[[m]]$Final = ft_lower[[m]][,1]
        
      }
    ###
    # Finally, testing to see if $Final is viable (usable) 
    ###
      if(is.na(mean(ft_mean[[m]]$Final))){
          models_to_calc_with = setdiff(models_to_calc_with, m)
      }
  } #ending the model for loop
  
  ###now we have a ft_mean list that contains multiple dataframes, each for a seperae modeling techinuq
  ###that has been cleaned and smoothed (if need be) and aggregated together we just gotta paste the bits together :D
  
    ft = data.frame(lapply(models, function(x) data.frame(ft_mean[[x]]$Final,
                                                     ft_upper[[x]]$Final,
                                                      ft_lower[[x]]$Final
                                                     )))
    names(ft) = #adding name vector
      as.vector( 
              unlist(
                    lapply(models, function(x) c(as.character(x), 
                                                 paste(c(x,"upper"), collapse ="_"),
                                                 paste(c(x,"lower"), collapse ="_")
                                                 )
                           )
                  )
            )
    
    ###
    # Adding in the average value of all the models (that survived) to acutally produce our forecast
    ###
    
    if(length(models_to_calc_with) > 1) {
        ft$x = rowMeans(ft[,models_to_calc_with], na.rm = T)
        ft$x_upper = rowMeans(ft[,match(models_to_calc_with, names(ft))+1], na.rm = T) 
        ft$x_lower = rowMeans(ft[,match(models_to_calc_with, names(ft))+2], na.rm = T)
    } else if (length(models_to_calc_with) == 1){
        ft$x = ft[,models_to_calc_with]
        ft$x_upper = ft[,match(models_to_calc_with, names(ft))+1] 
        ft$x_lower = ft[,match(models_to_calc_with, names(ft))+2]
    } else {stop("No models were found to work, yeeesh")}
    
    ft = ft[,c((ncol(ft)-2):ncol(ft), 1:(ncol(ft)-3))]  #reordering the columns
          
    ###
    # making the actual historical record correct
    ###
        if(!is.null(log_transform)){
          df = exp(df)
        }
    
      ft[1:NROW(df), c("x","x_upper","x_lower")] = df
      
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



