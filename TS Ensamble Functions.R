
library(dplyr)

########################################################################
# 
# this function holds the methods to run the TS analysis
#
########################################################################

models_runs = function(method, 
                       start, 
                       frequency, 
                       df, 
                       h, surpress.error = F, 
                       season = F,
                       ma_adj = 0, 
                       CI = 95, 
                       return_method = c('stack','list')) {
  if(ma_adj > 0){
    df = ma(df, n = ma_adj)
  }
  
  if(season & method != 'STL'){
    seas = stl(ts(df, start = start, frequency = frequency), s.window = h, robust = T)
    
    df = seas$time.series[,2] #for the trend compponent
    
    seas = list( seasonal = as.double(rep(seas$time.series[(nrow(seas$time.series) - frequency + 1):nrow(seas$time.series),1],
                     ceiling(h/frequency))[1:h]),
                 error = as.numeric(seas$time.series[nrow(seas$time.series),3])
                 )
  } else {
    seas = list(seasonl = 0, error = 0)
    
  }
  if (method == 'HoltWinters') {
    
    M = try(HoltWinters(ts(df, start = start, frequency = frequency)), silent = T)
    
  } else if (method == 'ARIMA'){
    for(d in 0:5)
    {
      print('try')
      M = try(arima(as.numeric(df), order = c(1,d,1), seasonal = list(order = c(1,d,1), period = frequency)),silent = T)
      if(!inherits(M,"try-error")) {break}
    }
    if(inherits(M,"try-error")) { #ERROR
      print("Couldn't find a non stationary means for ARIMA; switching to auto arima")
      M = auto.arima(df, parallel =  T, stepwise = F)
    } 
  }else if (method == 'STL'){
    M = try(stl(ts(df, start = start, frequency = frequency), s.window = "periodic", robust = T, s.degree = 1), silent = F)
  } else if  (method == 'TBATS') {
    M = try(tbats(ts(df, start = start, frequency = frequency), seasonal.periods = frequency), silent = F)
  } else if (method == 'lm') {
    sub_df = data.frame(y = df[(NROW(df) - 51):NROW(df)], t = 1:52, wt = 1:52/52)
    M = lm(y ~ t, data = sub_df, weight = sub_df$wt) 
    # M = lm(y ~ t, data = sub_df) 
    
  } else {
    
    M == 1
  }
  
  
  if(inherits(M,"try-error")) {
    if(!surpress.error){
      print(paste(c(method," Failed, will not use"), sep = ""))}
    df = 0
    frcst = 0
  } else if (method == 'HoltWinters'){
    
    df[(NROW(df) -  nrow(data.frame(M[1])) + 1):NROW(df)] = data.frame(M[1])[,1]
    df[is.na(df)]=df[is.na(df)]
    frcst = forecast(M,h=h, level = CI)
    
  } else if(method == 'ARIMA')
  {
    frcst = forecast(M, h=h, level = CI)
    df = df - as.numeric(M$residuals) 
    
  } else if(method == 'STL'){
    df = rowSums(data.frame(M[1])[,1:2])
    frcst = forecast(M,h=h,method = 'arima', level = CI)
    
  } else if(method == 'TBATS'){
    df = as.numeric(M$fitted.values)
    frcst = forecast(M, h=h, level = CI)
  } else if(method == 'lm'){
    lm_predict = as.double(predict(M, newdata = data.frame(t = 52:(52 + h))))
    # lm_predict = lm_predict + (df[NROW(df)] - lm_predict[1]) #adjusting it so that the last point is right on target
    z_score = qnorm(float(CI)/100)
    sd_df = sd(df)
    frcst = structure(list(
      mean = lm_predict[2:NROW(lm_predict)],
      lower = data.frame(lm_predict[2:NROW(lm_predict)] - z_score*sd_df
      ),
      upper = data.frame(lm_predict[2:NROW(lm_predict)] + z_score*sd_df
      )
    )
    )               
    
  } else {
    df
    z_score = qnorm(float(CI)/100)
    frcst =  structure(list(
      mean = rep(df[NROW(df)], h),
      lower = data.frame(rep(df[NROW(df)] - z_score*df(df), h)
      ),
      upper = data.frame(rep(df[NROW(df)] + z_score*df(df), h)
      )
    )
    )
    
  }
  if(return_method[1] == 'list') {
    output = structure(list(df = df, forecast = frcst))
  }
  else{
    
    print(str(as.numeric(frcst$mean +  seas$season + seas$error)))
    
    output = data.frame(mean = as.numeric(frcst$mean +  seas$season + seas$error),
                       upper = as.numeric(frcst$upper + seas$season + seas$error), 
                       lower = as.numeric(frcst$lower + seas$season + seas$error))

    df = data.frame(mean = as.numeric(df), 
                    upper = as.numeric(df), 
                    lower = as.numeric(df))
 
    output = rbind(df, output)
  }
  return(output)
  
}


########################################################################
# 
# this function applies a moving average to the data
#
########################################################################



ma = function(arr, n = 3){
	res = arr
	for(i in n:length(arr)){
		res[i] = mean(arr[(i-n+1):i])
	}
	return(res)
}

########################################################################
# 
# This function takes a vector, does some testing to it, and then 
# forecast out 4 diffrent TS series forecast
# averages the result together, and returns a data frame with the result
# and CI of the 4 individual forecasts
#
########################################################################





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
                                 Date_V, #just a date vector, used to calcuate which rows we need to use for blending
                                 num_of_elements = 1, 
                                 jumps = 1,
                                 datetype = 'months', 
                                 location =  0,
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
  #this will hold the forecasted values
    ft = data.frame(rep(0, h))
  #error flag incase any of the methods fail
    error_flag = FALSE
    
  names(ft) = names(df)
  #running the models
    for (m in models)
        {
          ret = models_runs(method = m, start= start, frequency = frequency, df = df[,1], h = h, ma_adj = ma_adj)
          
          
          frcst = ret$forecast
          
          if(NROW(ret$df) == 1){
                  models_to_calc_with = setdiff(models_to_calc_with, m)
                  error_flag = TRUE
            } else { #Not an error
                df[,m] =ret$df
                
            }
        
          
        #updating dataframes
       
        if(!error_flag){ #if there wasn't an error
         
            #if we want to do a blend check
              if(!is.na(blend_threshhold)) {
                temp_f =blend(Date_V = Date_V, 
                        df =  df[,1],
                        start = start, 
                        frequency = frequency,
                        method = m,
                        blend_threshhold = 0.05,
                        h = h,
                        datetype = datetype,
                        jumps = jumps,
                        location =location,
                        ma_adj = ma_adj,
                        frcst = data.frame(mean = frcst$mean,
                                           upper = frcst$upper[,2],
                                           lower = frcst$lower[,2]))
                
                          frcst_mean =  as.double(temp_f$mean)
                          frcst_upper = as.double(temp_f$upper)
                          frcst_lower = as.double(temp_f$lower)
                          
                        
              } else {
                
                frcst_mean = as.double(frcst$mean)
                frcst_upper = as.numeric(frcst$upper[,2])
                frcst_lower = as.numeric(frcst$lower[,2])
              }

            
            ft[,m] = frcst_mean
            ft[,paste(c(m,"upper"), collapse = "_" )] =  frcst_upper
            ft[,paste(c(m,"lower"), collapse = "_" )] =  frcst_lower
            df[,paste(c(m,"upper"), collapse = "_" )] =  df[,m]
            df[,paste(c(m,"lower"), collapse = "_" )] =  df[,m]
            
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




########################################################################
# 
# this function holds the methods to run the TS analysis
#
########################################################################


blend = function(Date_V, df, start = c(2001,1), frequency = 26, ma_adj = NA, method = 'HoltWinters', blend_threshhold, h = 26, datetype = "months", jumps = 1, location = 0, frcst)
{ 
      
      
      prev_mo = create_date_vector(EndDt =  max(Date_V),
                                       datetype = datetype,
                                       jumps = jumps, 
                                       location =  location,
                                       num_of_elements = 2)[2]
    
      period_diff = sum(Date_V > prev_mo)
      ret = models_runs(method = method, start= start, frequency = frequency, df = df[Date_V < prev_mo], h = (h + period_diff), surpress.error = T, ma_adj = ma_adj)  
      
      frcst2 = ret$forecast

      if(NROW(ret$df) == 1){ #ERROR
        return(frcst)
      } else {#IF NOT ERROR, which means we got it to work
    
      mean1 = frcst$mean[frequency*2]
      mean2 = frcst2$mean[frequency * 2 + period_diff]
  
      if(is.na(mean1) | is.na(mean2) | mean2 == 0) {return(frcst) # if we have issues with the means, just don't modify the new standings
 
      } else if (abs(mean1/mean2 - 1) < blend_threshhold) {
        return(frcst)
      } else {
        print("blending")
        print(abs(mean1/mean2 - 1) )
   
        adj_fcst_range = (NROW(frcst2$mean)- NROW(frcst$mean) + 1): NROW(frcst2$mean)

        frcst$mean = as.double((frcst$mean + frcst2$mean[adj_fcst_range])/2)
        frcst$upper = as.double((frcst$upper + as.numeric(frcst2$upper[adj_fcst_range,2]))/2)
        frcst$lower = as.double((frcst$lower + as.numeric(frcst2$lower[adj_fcst_range,2]))/2)
        

      }

      return(frcst)
      
     }
}
########################################################################
# 
# This function takes a length of period, then performance a ensablem_ts_forecast
# for that many pay periods and averages them all together to create a moving  'average'
#
########################################################################

ts_moving_average = function(x,
                             date_v, #you need to pass in a date vector, foo!
                             num_of_elements = 3, 
                            jumps = 1,
                             datetype = 'months', 
                             location =  0,
                             #The below elements are all controll factors for ensamble_ts_forecast
                                 start = c(1900, 1), 
                                 frequency = 26, #frequency of the data set
                                 h = 26, #how many time steps forward we want to use
                                 models = c("HoltWinters", "STL", "ARIMA", "TBATS","Naive"), #the models
                                 zero_na = T, #replace NA'/ Inf with 0's
                                 rm_na = F,
                                 CI_type = 'mean', #type used to calculate the CI
                                 use = c("Y","Y"), #Do we want to keep the historical or future parts
                                 scale = 'F', #do we want to scale the data?
                                 log_transform = NULL, #if you want to log transform the data first,
                                  models_to_calc_with = models,
                                  growth_limit = 0.01,
                                blend_threshhold = 0.07,
                            ma_adj = NA,
                                 ...){

  date_vector = create_date_vector(EndDt =  max(date_v),
                                   datetype = datetype ,
                                   jumps = jumps, 
                                   location =  location,
                                   num_of_elements = num_of_elements)
  Num_of_success = 0
  #time to itterate
  for (i in 1:num_of_elements){
    #subsetting the data
      print(paste(c("working on iteration",i,"of",num_of_elements),collapse = " "))
      x1 = x[date_v <= date_vector[i]]
      length_of_elements = NROW(x) - NROW(x1)
 
    #ensablme time  
      df_temp =  try(ensablem_ts_forecast( x = x1,
                              start = start,
                              frequency = frequency,
                              h = h + length_of_elements,
                              models =  models,
                              zero_na = zero_na,
                              rm_na = rm_na,
                              CI_type = CI_type,
                              use = use,
                              scale = scale,
                              log_transform = log_transform,
                              models_to_calc_with = models_to_calc_with,
                              growth_limit = growth_limit,
                              Date_V = date_v[date_v <= date_vector[i]],
                              num_of_elements = 3, 
                              jumps = 1,
                              datetype = 'months', 
                              location =  0,
                              blend_threshhold = 0.05,
                              ma_adj = ma_adj
                              ))
      if(inherits(df_temp,"try-error")) {next}
      Num_of_success= Num_of_success + 1
      #creating the DF that weill either hold all the data or using that df
      #to append growing values to it
      if(i == 1){
          master_df = df_temp
      } else {
        master_df[(NROW(x) +1):nrow(master_df), n_v(master_df)] = 
                              master_df[(NROW(x) +1):nrow(master_df), n_v(master_df)] +
                              df_temp[(NROW(x) +1):nrow(df_temp),n_v(df_temp)]
      }
      #testing the waters to see if we can do the next itteration (taknig the average change and scalring it up to the next level)      
        if(NROW(x1) - length_of_elements/i * (i+1)/i < frequency * 2 & i != num_of_elements ) {break}
        
        
      } #fin
    #now, averaging all those values together
    
    if(Num_of_success == 1){return(master_df) #if i = 1, then we are already done
    } else if(Num_of_success == 0) {return(x)
    } else{#if i > 1, then we need to average the non char rows together
      master_df[(NROW(x) +1):nrow(master_df), n_v(master_df)] = 
                                master_df[(NROW(x) +1):nrow(master_df), n_v(master_df)]/Num_of_success
      return(master_df)
    }
  
  }

########################################################################
# 
# This function Creates a vector of dates for whatever reason
#
########################################################################

create_date_vector = function(EndDt, 
                              datetype = 'months', 
                              jumps = 1, 
                              location = 1, # 1is the start of th month, 0 is end of the month
                              num_of_elements){
  seq(as.Date(eom(as.Date(EndDt)) + location),
      length = num_of_elements,
      by = paste(c('-',jumps, ' ',datetype), collapse = "")
      )
  
}
  



eom <- function(date_v) {
  next_date = seq(as.Date(date_v), length = 2, by = '1 months')[2]
  SOM = as.Date(paste(c(substr(next_date, 1, 8), '01'), collapse = ""))
  return(SOM - 1)
}

########################################################################
# 
# This function takes a data frame from the ensablem_ts_forecast
# and creates a confidence interabal for the final the ensamble forecast
#
########################################################################
ensamble_ts_CI = function(x, style = c("max-min","mean","median")){
  if(!is.data.frame(x)) {stop("x must be a data frame")}
  

  CI_upper_loc = sapply(names(x), function(x) substrRight(x, 5)) %in% 'upper'
  CI_lower_loc = sapply(names(x), function(x) substrRight(x, 5)) %in% 'lower'
  
  
  if(style == 'max-min'){
    
    if(sum(CI_upper_loc) + sum(CI_lower_loc) == 0)  {stop("You need to give me some CI to work with here, buddy")}
    upper = as.vector(apply(x[,CI_upper_loc], 1, max, na.rm = T))
    lower = as.vector(apply(x[,CI_lower_loc], 1, min, na.rm = T))
    
  } else if (style == "mean"){
    
    if(sum(CI_upper_loc) + sum(CI_lower_loc) == 0)  {stop("You need to give me some CI to work with here, buddy")}
    upper = as.vector(apply(x[,CI_upper_loc], 1, mean, na.rm = T))
    lower = as.vector(apply(x[,CI_lower_loc], 1, mean, na.rm = T))
    
  } else if (style == "median") {
    
    if(sum(CI_upper_loc) + sum(CI_lower_loc) == 0) {stop("You need to give me some CI to work with here, buddy")}
    upper = as.vector(apply(x[,CI_upper_loc], 1, median, na.rm = T))
    lower = as.vector(apply(x[,CI_lower_loc], 1, median, na.rm = T))
    
  } else {
    
    upper = x[,1] + sd(x[,1])*1.645
    lower = x[,1] - sd(x[,1])*1.645
    
  }
  
  
  x[,paste(c(names(x)[1],"upper"), collapse = "_" )] = upper
  x[,paste(c(names(x)[1],"lower"), collapse = "_" )] = lower  
  
  
  return(x)
}


########################################################################
# 
# These functions takes a data frame and returns a vector of either
# numeric variables or catecorical variables
#
########################################################################

n_v = function(df){ #n_v = Numerice_Vector
 test_if_df(df)
  
 vars = names(df)
 n_vars = vars[sapply(df[,vars], class) %in% c('numeric','integer')]
 
 return(n_vars)
}

c_v = function(df){ #c_v = Categorical_Vector
  test_if_df(df)
  
  vars = names(df)
  c_vars = vars[sapply(df[,vars], class) %in% c('factor','character')]
  
  return(c_vars)
}


########################################################################
# 
# This Function takes a data fraeme and two vectors, one of new fields
# and the other of thevalue to populate those fields with
#
########################################################################
add_df_values = function(df, field_names, field_values, add_date_vector = NULL){
  #testing values
  if(!is.vector(field_names)){ stop("field_names must be a vector")}
  if(!is.vector(field_values) & !is.data.frame(field_values)){ stop("field_values must be a values")}
  if(length(field_names) != length(field_values)) {stop("Field_names and field_values must be the same length")}
  if(length(field_names) == 0 | length(field_values) == 0) {stop('field_names and field_values MUST have length greater than zero....foo!')}
  test_if_df(df)
  
  for(i in 1:length(field_names))
  {
    df[,field_names[i]] = field_values[i]
    
  }
  
  if(!is.null(add_date_vector)){
    df[,names(add_date_vector)] =  add_date_vector
      }
  
  return(df)
}


########################################################################
# 
# This Function takes a data fraeme and and another vector, it then
# Multiples each numerical input of the DF by that vector
#
########################################################################
m_x_v = function(df, v, only_numeric = T){ #m_x_v = matrix times vector
  #testing values
    test_if_df(df)
    if(!is.vector(v)){stop("v needs to be a vector")}
    if(nrow(df) != length(v) & length(v) != 1) {stop('df and v need to be the same length, or v needs to be a scalar')}
  if(only_numeric)
  {
    df[,n_v(df)] =df[,n_v(df)] * v
  } else {
    df = df * v
  }
  
  
  return(df)
}

########################################################################
# 
# Some substring formulas to make my life just a little bit eiaser
#
########################################################################


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n){
  substr(x, 1, n)
}

########################################################################
# 
# a simple function to test if a passed data frame is, infact, a dataframe
#
########################################################################

test_if_df = function(df){
  if('data.table' %in% loadedNamespaces()) {
    if(!is.data.frame(df) & !is.data.table(df)) {stop("df isn't a data frame or data table")}
  } else if(!is.data.frame(df)) {stop("df isn't a data frame")}
  return(TRUE)
  
  
}


########################################################################
# 
# A Function that takes a vector and scales the interval into another range, either hard or soft scale
#
########################################################################

df_scale = function(vt, range = c(-Inf, Inf), vt_range = c(-Inf,Inf), method = c('Hard','Solid Shift','Scale', 'compress')){
      #HarD: Anything over or under the given values is hardcoded to the given values
      #Solid Shift: line is forced up or down so that the max/min point is right at the limit of the range (if the range is exceded)
     #scale  (b - a)(x - min)/(max - min) + a
      #factor by which the value exced the threshhold, then it's multipled the inverser of that value
    
  
    #Users ID10T check
        if(!is.vector(vt)) {stop("vt isn't a vector")}
        if(!is.numeric(vt)) {stop("vt must be numeric")}
        if(length(range) +  length(vt_range) != 4 ) {stop("range or vt_range must have two elements each")}
        if(length(method) > 1) {stop("please pick one method!!!!")}
    method = tolower(method)   
  
    if(sum(is.infinite(vt_range)) == 2 & method == 'Scale')
      {vt_range[1] = min(vt)
        vt_range[2] = max(vt)
      }
  
    
      if(sum(is.infinite(range)) == 2 ){
        #nothing to scale, ya dipshit
        
      } else if (method == 'hard'){
              if(range[1] != -Inf) {vt = ifelse( vt < range[1], range[1],  vt)}
              if(range[2] !=  Inf) {vt = ifelse( vt > range[2], range[2],  vt)}
      } else if (method == 'solid shift') {
          if(range[1] != - Inf){
              delta = range[1] - min(vt)  
              if(delta > 0){vt = vt + delta}
          } else if (range[2] != Inf) {
            delta = max(vt)  - range[2]
            if(delta > 0){vt = vt - delta}
          }
      
      } else if (method == 'scale') {
        if(is.infinite(range[1]) | is.infinite(range[2])) {stop("for the scale method, there needs to be hard limits")}
        vt = (range[2] - range[1]) * (vt - vt_range[1])/(vt_range[2] - vt_range[1]) + range[1]
      } else if(method == 'factor') {
          if(range[1] != -Inf) {
              if(range[1] == 0 ) {vt = df_scale(vt, range = c(0,Inf), method = 'Solid Shift')
              print('rescale')
              } else  {
              delta = abs(min(vt)/range[1])
              if( min(vt) < range[1] ) {vt = vt*(1/delta)}
              }
          }
          if(range[2] != Inf){
              delta = max(vt)/range[2]
              
              if(max(vt) > range[2]) {vt = vt*(1/delta)}
          }
           
      } 
  
  return(vt)
  
  
  
}



########################################################################
# 
# A function that will take a ts vector and decompose it into the 3 ranges, basically a wrapper for stl
#
########################################################################
stl_w = function(df, type = c("trend","seasonal","error"), s.window = "periodic", ...){
  t = stl(df, s.window = s.window, ...)$time.series
  if(length(type) > 1){type = type[1]}
  
  if(type == "seasonal") { ret = t[,1]}
  if(type == "trend") {ret = t[,2]}
  if(type == "error") {ret = t[,3]}
  
  return(ret)
}

########################################################################
# 
# A function that will take a vector of names and length and create an empty df around it
#
########################################################################






