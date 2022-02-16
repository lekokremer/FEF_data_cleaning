#Functions

#open_concat
#1. find files based on location and site
#2. concatonate all files in the folder
#3. add cap-rod no to each df
#4. eliminate duplicated dates and rename columns as needed

# interfiles is a subsection of the filepath that links the working directory
# and the location folder
# location is the catchment name and a folder that contains site folders

opn_concat <- function(interfiles, location, site) {
   file_path <- paste(getwd(),interfiles,location,site, sep='/')
   path_list <- paste(file_path, list.files(file_path), sep= '/')
   data <- lapply(path_list, function(x) {
      dat <- read.table(x, skip = 0, header = TRUE, sep = ",", row.names = NULL, as.is = TRUE)
      # for each item in path list, grab the cap_rod number
      dat$rod_no <- unlist(strsplit(x, "_"))[11]
      return(dat)
   })
   combined.data <- do.call(rbind, data)
   combined.data <- combined.data %>%
      mutate(datetime = lubridate::mdy_hm(formatted_datetime))%>%
      distinct()%>%
      arrange(datetime) 
   drops <- c("formatted_datetime")
   combined.data <- combined.data[ , !(names(combined.data) %in% drops)]
   rename(combined.data, water_temp = wtemp_a_1)
   return(combined.data)
}


#checks TimeSteps
checkTimeSteps<- function(df=stage_raw$datetime){
      checkts<- c(1,2+which(diff(diff(df))!=0))
      return(checkts)
}

###Dygraph for raw temperature data
DyTemp<- function(threshold='.2',airtemp='y'){
      tempRaw<- mutate(stage_raw, flag = ifelse(c(0,abs(diff(water_temp)))/water_temp > threshold, 48,1))
      if(airtemp == 'y'){
            tsTemp<- xts(dplyr::select(tempRaw, datetime,ID,water_temp,logger_temp), order.by=tempRaw$datetime)
      }else if(airtemp=='n'){
            tsTemp<- xts(dplyr::select(tempRaw, datetime,ID,water_temp), order.by=tempRaw$datetime)
      }
      #launch dygraph using xts object
      dygraph(tsTemp) %>% 
            #adds time series launcher
            dyRangeSelector() %>% 
            #adds highlight/fade (controled by alpha) of series and formats circle.
            #dyHighlight(highlightCircleSize = 4, 
            #            highlightSeriesBackgroundAlpha = 0.2,
            #            hideOnMouseOut = TRUE)%>% 
            dyAxis('y',label='Degree C',valueRange = c(-10,50 ))%>%
            dyAxis('y2',label='ID',independentTicks=T)%>%
            dySeries('ID',axis='y2')%>%
            #Assigns legend to follow cursor, also can choose 'always' to always see it
            dyLegend(show = "always")
}         

#Dygraph for raw stage - option to flag deviations by slope change 
DyRawStage<- function(df=stage_raw_prep,threshold = 0.2, flag='TRUE',max=1200){
      if (flag){
            stage_df<- mutate(df, flag = ifelse(c(0,abs(diff(wtr_ht_avg)))/wtr_ht_avg > threshold, 1000,100))
            tsStage<- xts(dplyr::select(stage_df, datetime,ID,flag,wtr_ht_pt,wtr_ht_avg), order.by=stage_df$datetime)
            dygraph(tsStage) %>% 
                  dyAxis('y',label='mm',valueRange = c(-100, max))%>%
                  dyAxis('y2',label='ID',independentTicks=T)%>%
                  dySeries('ID',axis='y2')%>%
                  dyRangeSelector() %>%
                  #dyHighlight(highlightCircleSize = 4, 
                  #        highlightSeriesBackgroundAlpha = 0.2,
                  #        hideOnMouseOut = TRUE)%>%
                  #dyOptions(drawPoints = TRUE, pointSize = 2)%>%
                  dyLegend(show = "always")
      }else{
            tsStage<- xts(dplyr::select(stage_df, datetime,ID,water_ht1,water_avg), order.by=stage_df$datetime)
            dygraph(tsStage) %>% 
                  dyAxis('y',label='mm',valueRange = c(-100, max))%>%
                  dyAxis('y2',label='ID',independentTicks=T)%>%
                  dySeries('ID',axis='y2')%>%
                  dyRangeSelector() %>%
                  #dyHighlight(highlightCircleSize = 4, 
                  #        highlightSeriesBackgroundAlpha = 0.2,
                  #        hideOnMouseOut = TRUE)%>%
                  #dyOptions(drawPoints = TRUE, pointSize = 2)%>%
                  dyLegend(show = "always")
      }
}

DyBatt<- function(){
      tsBatt<- xts(dplyr::select(stage_raw, datetime,batt_pw), order.by=stage_raw$datetime)
      dygraph(tsBatt) %>% 
            dyRangeSelector()%>%
            dyHighlight(highlightCircleSize = 4, 
                        highlightSeriesBackgroundAlpha = 0.2,
                        hideOnMouseOut = TRUE)%>%
            dyOptions(drawPoints = TRUE, pointSize = 2)%>%
            dyLegend(show = "always")
}

AdjStage <- function(df=stage_raw,maxgap=8){
      #create Stage Adj dataframe and adj wt ht column
      stageAdj<-df%>%
            mutate(adj_wtr_ht = wtr_ht_avg)
      #for loop that assigns corrected offset within ranges of IDs define in Cor
      for(i in 1:length(vert_correction$ID)){
            if(i<length(vert_correction$ID)){
                  #assigns final ID for each range
                  x<- vert_correction$ID[i+1]-1
                  #adds cumulative offset
                  stageAdj$adj_wtr_ht[vert_correction$ID[i]:x] <-stageAdj$wtr_ht_avg[vert_correction$ID[i]:x]+vert_correction$cumOffset[i]
            }
            #for last set of IDs, assigns the window from the final ID in vert_correction until end of TS
            if(i==length(vert_correction$ID)){
                  stageAdj$adj_wtr_ht[vert_correction$ID[i]:length(stageAdj$wtr_ht_avg)]<-
                        stageAdj$wtr_ht_avg[vert_correction$ID[i]:length(stageAdj$wtr_ht_avg)]+vert_correction$cumOffset[i]
                  
            }
      }
      #deletes Bad IDs and interpolates between then if gap is less than 8 in a row
      stageAdj <- mutate(stageAdj, adj_wtr_ht = na.approx(ifelse(ID %in% bad_id,NA, adj_wtr_ht),maxgap=maxgap,na.rm=F))
      return(stageAdj)    
}

dyStageAdj<- function(df= stageAdj,max=1200){
      tsStageAdj<- xts(dplyr::select(df, datetime,wtr_ht_avg,adj_wtr_ht,ID), order.by=df$datetime)
      dygraph(tsStageAdj) %>% 
            dyAxis('y',label='mm',valueRange = c(-150, max))%>%
            dyAxis('y2',label='ID',independentTicks=T)%>%
            dySeries('ID',axis='y2')%>%
            dyRangeSelector() %>%
            #dyHighlight(highlightCircleSize = 4, 
            #        highlightSeriesBackgroundAlpha = 0.2,
            #        hideOnMouseOut = TRUE)%>%
            dyOptions(drawPoints = FALSE, pointSize = 2)%>%
            dyLegend(show = "always")
}

#join to main stage database and interpolate offsets and calculate final stage (relative to reference gage position)
interpStage <- function(){
      
      stageCor<- left_join(stageAdj,stageOffset)%>%
            arrange(datetime)
      stageCor$offset[1] <- 0
      stageCor<-mutate(stageCor,interp_offset= ifelse(datetime <= lastdate,
                                                      na.approx(offset,na.rm=F),last_offset))%>%
            mutate(interp_offset=ifelse(is.na(interp_offset),0,interp_offset))%>%
            mutate(final_stage = adj_wt_ht+interp_offset)
      return(stageCor)
}
