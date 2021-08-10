
rm(list=ls()); gc()

library(tidyverse)
library(lubridate)
library(readstata13)


make_life_table<-function(dat){
  ### create proportions of pop with each outcome by age/year
  dat<-dat%>%
    mutate(m = var / pop)
  
  ### convert to probability 
  ### age_period (n) = 1 for all cases
  ### a = 0.5 (avg age of TPR for ppl in year, within-period survival)
  dat<-dat%>%
    mutate(q = 1 * m / (1 + (1 - 0.5) * m),
           p = 1 - q)
  ### make cumulative survival
  dat<-dat%>%
    mutate(lx = 1e5 * cumprod(c(1, p))[-nrow(dat)])
  ### deaths
  dat<-dat%>%
    mutate(d = -c(diff(lx),0))
  ## person-years in each group
  dat<-dat%>%
    mutate(L = (lx - d) * 1 + d * 0.5,
           t = sum(L)- cumsum(L) + L)
  ## life expectancy (time to TPR)
  dat<-dat%>%
    mutate(e = t/lx)
  ### cum prevalence
  dat<-dat%>%
    mutate(c = 1-lx/1e5)
  
  dat<-dat %>% 
    mutate(s = 1 - c) %>% 
    mutate(var = as.numeric(var),
           pop = as.numeric(pop)) %>% 
    mutate(se = s * 
             sqrt(
               sum(var / 
                     (pop * (pop - var)))
             ))
  
  dat<-dat %>% 
    mutate(se = ifelse(age==0, 0, se)) # set se to 0 when age == 0 for variance on s, c=1-s
  
  return(dat)
}


## population data are available here: https://seer.cancer.gov/popdata/download.html
pop<-read_fwf("./data/us.1990_2019.singleages.adjusted.txt",
              fwf_widths(c(4, 2, 2, 3, 2, 1, 1, 1, 2, 8),
                         c("year", "state", "st_fips",
                           "cnty_fips", "reg", "race",
                           "hisp", "sex", "age", "pop"))) 


## filter population for 2014-2018
pop<-pop%>%
  mutate(year = as.integer(year)) %>% 
  filter(year>=2014) %>% 
  filter(year<=2018)%>%
  mutate(pop = as.integer(pop))%>%
  mutate(race_ethn =
           case_when(
             race==1 & hisp ==0 ~ "White",
             race==2 ~ "Black",
             race==3 ~ "AI/AN",
             race==4 ~ "Asian/PI",
             hisp==1 ~ "Hispanic")) 


## create fipscode
pop_fips<-pop %>%
  mutate(fipscode = as.numeric(
    paste(st_fips, cnty_fips,
          sep = "")))%>% 
  mutate(fipscode =
           case_when(fipscode==36047 ~ 36061,
                     fipscode==36005 ~ 36061,
                     fipscode==36081 ~ 36061,
                     fipscode==36085 ~ 36061,
                     T ~ fipscode)) ## mega new york to match afcars


## filter those 20 most populated counties
top_pops<-pop_fips %>%
  ungroup() %>% 
  group_by(fipscode) %>%
  summarise(pop = sum(pop)) %>%
  arrange(desc(pop)) %>%
  ungroup() %>%
  mutate(rank = 1:n()) %>%
  filter(rank<=20) 

write_csv(top_pops, "./data/top_pops.csv")


## of the top 20 counties, we focus only those in california
pop_fips <- pop_fips %>% 
  ungroup() %>% 
  filter(state=="CA") %>%
  filter(fipscode%in%top_pops$fipscode) %>% 
  filter(age<=18) %>% 
  mutate(age = as.integer(age)) %>% 
  group_by(year, fipscode, age, race_ethn) %>% 
  summarise(pop = sum(pop)) %>% 
  ungroup()


## load first CPS contact data

inv<-read.dta13("./data/first_inv_cws.dta") %>% 
  rename(year = inv_year, 
         fipscode = rptfips,
         race_ethn = race1,
         age = inv_age) %>% 
  filter(age<19) %>%
  group_by(year, fipscode, race_ethn, age) %>% 
  summarise(first_inv = sum(inv)) %>% 
  ungroup()

sub<-read.dta13("./data/first_sub_cws.dta") %>% 
  rename(year = sub_year, 
         fipscode = rptfips,
         race_ethn = race1,
         age = sub_age) %>% 
  filter(age<19) %>%
  group_by(year, fipscode, race_ethn, age) %>% 
  summarise(first_sub = sum(sub)) %>% 
  ungroup()

fc<-read.dta13("./data/first_fc_cws.dta") %>% 
  rename(year = fc_year, 
         fipscode = rptfips,
         race_ethn = race1,
         age = fc_age) %>% 
  filter(age<19) %>%
  group_by(year, fipscode, race_ethn, age) %>% 
  summarise(first_fc = sum(fc)) %>% 
  ungroup()

tpr<-read.dta13("./data/first_tpr_cws.dta") %>% 
  rename(year = tpr_year, 
         fipscode = rptfips,
         race_ethn = race1,
         age = tpr_age) %>% 
  filter(age<19) %>%
  group_by(year, fipscode, race_ethn, age) %>% 
  summarise(first_tpr = sum(tpr)) %>% 
  ungroup()


## full join because california cws source data were not imputed for race/ethnicity.
dat <-pop_fips %>%
  full_join(inv) %>%
  full_join(sub) %>%
  full_join(fc) %>%
  full_join(tpr)

dat[is.na(dat)]<-0

check<-dat %>% 
  group_by(year, fipscode, race_ethn) %>% 
  summarise(pop = sum(pop),
            first_inv = sum(first_inv),
            first_sub = sum(first_sub),
            first_fc = sum(first_fc),
            first_tpr = sum(first_tpr))


## set up for lifetable loop
tab_dat<-dat %>% 
  pivot_longer(names_to = "varname",
               values_to = "var",
               cols = first_inv:first_tpr)


## set up as 5 year total life table
tab_dat<-tab_dat %>% 
  ungroup() %>% 
  group_by(fipscode, age, race_ethn, varname) %>% 
  summarise(pop = sum(pop),
            var = sum(var)) %>% 
  ungroup()

tot_dat<-tab_dat %>% 
  group_by(fipscode, age, varname) %>% 
  summarise(pop = sum(pop),
            var = sum(var)) %>% 
  mutate(race_ethn = "Total")

tab_dat<-tab_dat %>% 
  bind_rows(tot_dat)

tab_dat <- tab_dat %>%
  filter(pop != 0)


### run life tables by race, sex
vars<-unique(tab_dat$varname)
race<-unique(tab_dat$race_ethn)
fips<-unique(tab_dat$fipscode)
tables_out<-list()


### descriptive table out for appendix
appx<-tab_dat %>% 
  group_by(fipscode, race_ethn, varname) %>% 
  summarise(pop = sum(pop), var = sum(var)) %>% 
  pivot_wider(names_from = varname, values_from = var) %>% 
  mutate(county =
           case_when(
             fipscode==6001 ~ "CA: Alameda",
             fipscode==6037 ~ "CA: Los Angeles",
             fipscode==6059 ~ "CA: Orange",
             fipscode==6065 ~ "CA: Riverside",
             fipscode==6071 ~ "CA: San Bernadino",
             fipscode==6073 ~ "CA: San Diego",
             fipscode==6085 ~ "CA: Santa Clara"
           )) %>% 
  filter(!is.na(county)) %>% ### keep only california data
  ungroup() %>% 
  select(-fipscode) %>% 
  select(county, race_ethn, pop, first_inv, first_sub, first_fc, first_tpr) %>% 
  write_csv("./data/appx_pop_var_table.csv")

library(xtable)
appx<-appx %>% 
  select(county, race_ethn, pop, first_inv, first_sub, first_fc, first_tpr) %>% 
  rename(County = county,
         `Race/ethnicity` = race_ethn,
         `5 year population total` = pop,
         `First investigation` = first_inv,
         `First substantiation` = first_sub,
         `First placement` = first_fc,
         `First TPR` = first_tpr)

appx_out<-xtable(appx,
                 digits = 0,
                 caption = "First events by county and race, 5-year pooled event counts and population")

print(appx_out,
      include.rownames = F,
      file = "appx_tab.tex"
)

counter<-0
for(h in 1:length(vars)){
  for(r in 1:length(fips)){
    for(y in 1:length(race)){
      counter<-counter + 1
      
      temp<-tab_dat %>%
        filter(varname == vars[h],
               fipscode == fips[r],
               race_ethn == race[y])
      
      tables_out[[counter]]<-make_life_table(temp)
    }
  }
}

tables<-bind_rows(tables_out)


## combine tables

tables_within<-tables %>%
  group_by(race_ethn, fipscode, age, varname) %>%
  summarise(c_mn = mean(c),
            v_within = mean(se^2))

tables_between<-tables %>%
  left_join(tables_within) %>%
  group_by(race_ethn, fipscode, age, varname) %>%
  summarise(v_between = mean((c - c_mn)^2))

tables_comb<-tables_within %>%
  left_join(tables_between) %>%
  mutate(se_tot = sqrt(v_within + (1 + 1/5)*v_between)) %>%
  select(fipscode, varname, race_ethn, age, c_mn, se_tot)

tables_comb<-tables_comb %>%
  mutate(c_upr = c_mn + 1.96 * se_tot,
         c_lwr = c_mn - 1.96 * se_tot) %>%
  filter(age == 17) %>%
  mutate(c_upr = ifelse(c_upr>1, 1, c_upr),
         c_lwr = ifelse(c_lwr<0, 0, c_lwr))


### format names
tables_comb<-tables_comb %>%
  mutate(county =
           case_when(
             fipscode==6001 ~ "CA: Alameda",
             fipscode==6037 ~ "CA: Los Angeles",
             fipscode==6059 ~ "CA: Orange",
             fipscode==6065 ~ "CA: Riverside",
             fipscode==6071 ~ "CA: San Bernadino",
             fipscode==6073 ~ "CA: San Diego",
             fipscode==6085 ~ "CA: Santa Clara"
           )) %>%
  filter(is.na(county) == FALSE) %>%
  filter(is.na(race_ethn) == FALSE)


write_csv(tables_comb, "county_tables.csv")

