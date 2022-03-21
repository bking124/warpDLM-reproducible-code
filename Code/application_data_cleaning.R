#Load needed packages
library(magrittr)
library(lubridate)
library(tidyverse)

#Import datasets
full_df <- read_csv("./Data/Cincinnati_Fire_Incidents.zip")
incidentTypes <- unique(full_df$INCIDENT_TYPE_ID)
full_df$DATE <-  full_df$CREATE_TIME_INCIDENT %>% 
  as.POSIXct(., format="%m/%d/%Y %I:%M:%S %p", tz="EST") %>%
  as.Date(., tz="EST")
full_dates <- full_df$DATE %>% unique() %>% sort()

#Take only incidents between 2017-01-01 and 2021-02-01
full_dates <- subset(full_dates, full_dates < '2021-02-01' & full_dates >= '2017-01-01')
full_df <- subset(full_df, DATE < '2021-02-01' & DATE >= '2017-01-01')

#Useful functions for filling in zeroes on missing dates
checkDateFunc <- function(x, nonZeros){
  sub <- subset(nonZeros, dates==x)
  result <- ifelse(nrow(sub)==0, 0, sub$Freq)
  return(result)
}

countDFProcessing <- function(df){
  df$CREATE_TIME_INCIDENT %<>% as.POSIXct(format="%m/%d/%Y %I:%M:%S %p", tz="EST")
  df$Date <- as.Date(df$CREATE_TIME_INCIDENT, tz="EST")
  dates <- sort(df$Date)
  fullDateRange <- seq(from=dates[1], to=dates[length(dates)], by="day")
  nonZeroCounts  <- as.data.frame(table(dates))
  nonZeroCounts$dates %<>% as.Date()
  
  dfCounts = data.frame(Date=fullDateRange)
  dfCounts$Count <- fullDateRange %>% sapply(., FUN= checkDateFunc, nonZeros=nonZeroCounts)
  dfCounts$Month <- month(dfCounts$Date)
  return(dfCounts)
}

#Get heroin and overdose incidents
overdoseIDs <- incidentTypes[grep('23', incidentTypes)]
heroin_responses <- c("HEROIF", "HEROIN-COMBINED", "HEROINF - FIRE ONLY", 
                      "HEROINF-FIRE ONLY", "HERON F", "HERONF")
heroin <- subset(full_df, INCIDENT_TYPE_ID %in% heroin_responses)
otherOD <- subset(full_df, INCIDENT_TYPE_ID %in% overdoseIDs)

heroinCounts <- countDFProcessing(heroin)
otherCounts <- countDFProcessing(otherOD)


#Do some basic data viz
countsByMonth <- group_by(heroinCounts, Month) %>% summarise(n=sum(Count))
barplot(countsByMonth$n, names.arg=countsByMonth$Month)

ggplot(data=heroinCounts,aes(x=Date, y=Count)) + geom_col()
ggplot(data=otherCounts,aes(x=Date, y=Count)) + geom_col()

ggplot(data=heroinCounts,aes(x=Date, y=Count)) + geom_line()
ggplot(data=otherCounts,aes(x=Date, y=Count)) + geom_line()

mean(heroinCounts$Count)
var(heroinCounts$Count)
length(unique(heroinCounts$Count))/length(heroinCounts$Count)
table(heroinCounts$Count)[1]/length(heroinCounts$Count)

mean(otherCounts$Count)
var(otherCounts$Count)
length(unique(otherCounts$Count))/length(otherCounts$Count)
table(otherCounts$Count)[1]/length(otherCounts$Count)

#Save datasets
save(heroinCounts,file="./Data/heroinCallsCounts")
save(otherCounts,file="./Data/drugCallsCounts")