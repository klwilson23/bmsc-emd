library(dplyr)
library(ggplot2)
lice <- read.csv("Data/sealice_fish_data.csv")
sites <- read.csv("Data/sealice_site_data.csv")

lice <- lice[!is.na(lice$length),] %>%
  group_by(site_id,year,day,month,location,species) %>%
  mutate(n=n()) %>%
  ungroup()

lice_data <- merge(lice,sites,by=c("site_id","year","day","month","location"),all.x=TRUE)
lice_sub <- lice_data %>%
  filter(species=="sockeye")
plot(chalA~length,data=lice_sub)
temp_salt <- sites[!is.na(sites$salt) & !is.na(sites$temp) & !is.na(sites$P_ratio) & !is.na(sites$C_ratio),]
plot(salt~temp,temp_salt)
