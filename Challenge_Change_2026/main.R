#main analysis code for PCAN2

rm(list=ls())
cat("\014")

library(readxl)
library(lme4)
library(lmerTest)
library(yhat)
#SEM related libraries
library(lavaan)
library(lavaanPlot)

source('- functions/z_score.R')

## read in PCA outputs, demographics, other data. combine.
ind = read.table('VWP/PCAindices.txt', header = T, sep = ',')
q_subs = unique( ind[ind$cond == 'Quiet','sub'] )
n_subs = unique( ind[ind$cond == 'Noise','sub'] )
v_subs = unique( ind[ind$cond == 'Vocoding','sub'] )

temp_shared_1 = n_subs[(n_subs %in% q_subs)]
temp_shared_2 = v_subs[(v_subs %in% q_subs)]
# this set of subjects has data for quiet, day1, and day2.
# same set of subjects for each challenge condition, so can just use one
data = ind[ind$sub %in% temp_shared_1,]
#sub 34 doesn't have noise data. drop
data = data[data$sub != 34,]


#read in demographic / other data
demos = read_excel('demos.xlsx')

ieee = read.table('IEEE/avg_scores.csv',header=T)
colnames(ieee) = c('ID','PBSB_ID','noise','ss','voc')
ieee$noise_z = z_score(ieee$noise)
ieee$ss_z = z_score(ieee$ss)
ieee$voc_z = z_score(ieee$voc)

stroop = read.table('Stroop/avg_scores.csv', header = T)
stroop$cc_z = z_score(stroop$score)
colnames(stroop)[1] = 'Stroop_ID'

wj = read_excel('WJ Oral Comp/wj_scores.xlsx')
wj$lang_z = z_score(wj$avg)
colnames(wj)[1] = 'PBSB_ID'


## more set up -- get all data in one df 
d1_data = data[data$day==1,]
d2_data = data[data$day==2,]
d2_data = d2_data[d2_data$sub!=11,]#sub 11 doesn't have both noise and vocoding on day2, drop

d1_subs = unique(d1_data$sub)
d2_subs = unique(d2_data$sub)
temp_d1 = d1_data[d1_data$sub %in% d2_subs,]
dall_data = rbind(temp_d1, d2_data)

qall_data = dall_data[dall_data$cond=='Quiet',]
colnames(qall_data) = c('PC1_q','PC2_q','PC3_q','sub','day','cond')
nall_data = dall_data[dall_data$cond=='Noise',]
colnames(nall_data) = c('PC1_n','PC2_n','PC3_n','sub','day','cond')
vall_data = dall_data[dall_data$cond=='Vocoding',]
colnames(vall_data) = c('PC1_v','PC2_v','PC3_v','sub','day','cond')

nall_d1 = nall_data[nall_data$day==1,]
nall_d2 = nall_data[nall_data$day==2,]
vall_d1 = vall_data[vall_data$day==1,]
vall_d2 = vall_data[vall_data$day==2,]
colnames(nall_d1) = c('PC1_d1_N', 'PC2_d1_N','PC3_d1_N','sub','day','cond')
colnames(nall_d2) = c('PC1_d2_N', 'PC2_d2_N','PC3_d2_N','sub','day','cond')
colnames(vall_d1) = c('PC1_d1_V', 'PC2_d1_V','PC3_d1_V','sub','day','cond')
colnames(vall_d2) = c('PC1_d2_V', 'PC2_d2_V','PC3_d2_V','sub','day','cond')

sem_pc = merge(nall_d1[,c('PC1_d1_N', 'PC2_d1_N','sub')],
               nall_d2[,c('PC1_d2_N', 'PC2_d2_N','sub')],
               by='sub')
sem_pc = merge(sem_pc,
               vall_d1[,c('PC1_d1_V', 'PC2_d1_V','sub')],
               by='sub')
sem_pc = merge(sem_pc,
               vall_d2[,c('PC1_d2_V', 'PC2_d2_V','sub')],
               by='sub')
sem_pc = merge(sem_pc,
               qall_data[,c('PC1_q', 'PC2_q','sub')],
               by='sub')


#drop language and activation rate from nlvwp -- dont have any relevance
new_outs = merge(demos[,c('PBSB_ID','Stroop_ID','age','PTA')], ieee, by = 'PBSB_ID')
new_outs = merge(new_outs[,c('PBSB_ID','Stroop_ID','age','PTA','noise_z','ss_z','voc_z')], 
             stroop[,c('Stroop_ID','cc_z')], by = 'Stroop_ID')
new_outs = new_outs[,c(2:ncol(new_outs))]
colnames(new_outs)[1] = 'sub'
new_outs$pta_z = z_score(new_outs$PTA)
new_outs$age_z = z_score(new_outs$age)

sem_dat = merge(new_outs, sem_pc, by='sub')


### find out whether challenge is sig different from quiet or not ###
d1_dat = dall_data[dall_data$day==1,]
d1_dat$sub = as.factor(d1_dat$sub)
d1_dat$cond = as.factor(d1_dat$cond)

contrasts(d1_dat$cond) = matrix(c(1/3,-2/3,1/3,
                                  -1,0,1), 
                                ncol=2)

chall_m_1 = lmer(PC1 ~ cond + (1|sub), d1_dat)
summary(chall_m_1)

chall_m_2 = lmer(PC2 ~ cond + (1|sub), d1_dat)
summary(chall_m_2)


#### Question 1 ####
#does quiet act as anchor to subregion of space, accounting for some covariates?
q_data = merge(dall_data[dall_data$cond != 'Quiet',], qall_data[,c('PC1_q','PC2_q','sub')], by='sub')

q_data$sub = as.factor(q_data$sub)
q_data$cond = as.factor(q_data$cond)
q_data$day = as.factor(q_data$day)
contrasts(q_data$day) = c(-1,1)
contrasts(q_data$cond) = c(-1,1)

q_pc1_mod1 = lmer( PC1 ~ (cond * day) + (PC1_q + PC2_q) * (cond + day) + (cond+day|sub),  q_data )
summary(q_pc1_mod1)

q_pc2_mod1 = lmer( PC2 ~ (cond * day) + (PC1_q + PC2_q) * (cond + day) + (cond+day|sub),  q_data )
summary(q_pc2_mod1)


#post-hoc analyses of interaction
noise_dat = q_data[q_data$cond=="Noise",]
voc_dat = q_data[q_data$cond=="Vocoding",]

post_hoc_1 = lmer( PC1 ~ day * (PC1_q + PC2_q) + (1|sub), noise_dat)
summary(post_hoc_1)
post_hoc_2 = lmer( PC1 ~ day * (PC1_q + PC2_q) + (1|sub), voc_dat)
summary(post_hoc_2)


#further analysis, including interactions between cc/demos and factors
q1_sub_data = merge(q_data, new_outs, by='sub')
q1_sub_data$sub = as.factor(q1_sub_data$sub)
q1_sub_data$cond = as.factor(q1_sub_data$cond)
q1_sub_data$day = as.factor(q1_sub_data$day)
contrasts(q1_sub_data$day) = c(-1,1)
contrasts(q1_sub_data$cond) = c(-1,1)

pc1_demos = lmer( PC1 ~ (cc_z + pta_z + age_z) * (cond + day) + cond * day + (cond + day | sub),  q1_sub_data )
summary(pc1_demos)

pc1_demos_quiet = lmer( PC1 ~ cc_z * (cond + day + PC1_q + PC2_q) + (cond + day | sub),  q1_sub_data )
summary(pc1_demos_quiet)

pc2_demos = lmer( PC2~ (cc_z + pta_z + age_z)* (cond + day) + cond * day + (cond + day | sub),  q1_sub_data )
summary(pc2_demos)

pc2_demos_quiet = lmer( PC2 ~ cc_z * (cond + day + PC1_q + PC2_q) + (cond + day | sub),  q1_sub_data )
summary(pc2_demos_quiet)


#cognitive control predicting quiet location directly? 
cc_q1 = lm(PC1_q ~ cc_z + age_z + pta_z, sem_dat)
summary(cc_q1)

cc_q2 = lm(PC2_q ~ cc_z + age_z + pta_z, sem_dat)
summary(cc_q2)

#### Question 2 ####
#quiet accounts for sig variance. 
#Beyond this though, is there some consistency in location when we just look at noise / day?
#so outside of whether quiet acts as an anchor, what is the extent of the local region in which people shift around?
#do people do the same sort of thing across noise types, basically...


q2_n_pc1 = aps(sem_pc,"PC1_d1_N",list("PC1_d1_V","PC1_d2_N","PC1_d2_V"))
q2_n_pc1 = commonality(q2_n_pc1)

q2_v_pc1 = aps(sem_pc,"PC1_d1_V",list("PC1_d1_N","PC1_d2_V","PC1_d2_N"))
q2_v_pc1 = commonality(q2_v_pc1)

q2_n_pc2 = aps(sem_pc,"PC2_d1_N",list("PC2_d1_V","PC2_d2_N","PC2_d2_V"))
q2_n_pc2 = commonality(q2_n_pc2)

q2_v_pc2 = aps(sem_pc,"PC2_d1_V",list("PC2_d1_N","PC2_d2_V","PC2_d2_N"))
q2_v_pc2 = commonality(q2_v_pc2)

#extend to Question 1
#regression by quiet alone
q2_n_p1_reg = lm(PC1_d1_N~PC1_q+PC2_q,sem_pc)
q2_v_p1_reg = lm(PC1_d1_V~PC1_q+PC2_q,sem_pc)
q2_n_p2_reg = lm(PC2_d1_N~PC1_q+PC2_q,sem_pc)
q2_v_p2_reg = lm(PC2_d1_V~PC1_q+PC2_q,sem_pc)

#include all other factors
q2_n_p1_reg_q = lm(PC1_d1_N~PC1_d1_V+PC1_d2_N+PC1_d2_V+PC1_q+PC2_q,sem_pc)
q2_v_p1_reg_q = lm(PC1_d1_V~PC1_d1_N+PC1_d2_V+PC1_d2_N+PC1_q+PC2_q,sem_pc)
q2_n_p2_reg_q = lm(PC2_d1_N~PC2_d1_V+PC2_d2_N+PC2_d2_V+PC1_q+PC2_q,sem_pc)
q2_v_p2_reg_q = lm(PC2_d1_V~PC2_d1_N+PC2_d2_V+PC2_d2_N+PC1_q+PC2_q,sem_pc)


#### Q3: outcomes ####

## first we want to predict outcomes by general location
modelraw <- '
                PC1 =~ PC1_d1_N + PC1_d2_N + PC1_d1_V + PC1_d2_V
                PC2 =~ PC2_d1_N + PC2_d2_N + PC2_d1_V + PC2_d2_V
                PC1 ~~ PC2

                #within challenge measurement covariance
                PC1_d1_N ~~ PC1_d2_N
                #PC1_d1_V ~~ PC1_d2_V  
                PC2_d1_N ~~ PC2_d2_N
                PC2_d1_V ~~ PC2_d2_V  

                #within day measurement covariance
                PC1_d1_N ~~ PC1_d1_V
                PC1_d2_N ~~ PC1_d2_V
                PC2_d1_N ~~ PC2_d1_V
                PC2_d2_N ~~ PC2_d2_V

                PC2 ~ cc_z

                noise_z + voc_z + ss_z ~ age_z + pta_z + PC1 + PC2

                age_z ~~ pta_z
'
#fitting 
fitraw <- sem(modelraw, data = sem_dat, effect.coding="loadings",missing="FIML")
#summary of the fits
summary(fitraw,fit.measures=TRUE,  standardized=T, rsquare=T)


# SEM Model of change in PCs on outcomes
modelchg <- '   PC1 =~ PC1_d1_N + PC1_d2_N + PC1_d1_V + PC1_d2_V
                PC2 =~ PC2_d1_N + PC2_d2_N + PC2_d1_V + PC2_d2_V
                PC1 ~~ PC2

                #within challenge measurement covariance
                PC1_d1_N ~~ PC1_d2_N
                #PC1_d1_V ~~ PC1_d2_V  
                PC2_d1_N ~~ PC2_d2_N
                PC2_d1_V ~~ PC2_d2_V  

                #within day measurement covariance
                PC1_d1_N ~~ PC1_d1_V
                PC1_d2_N ~~ PC1_d2_V
                PC2_d1_N ~~ PC2_d1_V
                PC2_d2_N ~~ PC2_d2_V

                #include additional quiet dimensions
                PC1 ~ PC1_q + PC2_q
                PC2 ~ PC2_q + PC1_q + cc_z

                PC1_q ~~ PC2_q

                #regress outcomes by dimensions and quiet!
                noise_z + voc_z + ss_z ~ age_z + pta_z + PC1 + PC2 + PC1_q + PC2_q

                age_z ~~ pta_z

'
#fitting 
fit_chg_all <- sem(modelchg, data = sem_dat, effect.coding="loadings",missing="FIML")
#summary of the fits
summary(fit_chg_all,fit.measures=TRUE,  standardized=T, rsquare=T)


## we show a negative effect of PC2 on noise outcomes
## this is not what we found in previous work, which saw a positive effect of these profiles in interaction with speeded sentences
## can't run an interaction of latent factor within SEM model, but can pull out latent factor values and do an interaction there to predict noise / vocoding outcomes

#extract latent factor values
lf_mod = 
  '
                PC1 =~ PC1_d1_N + PC1_d2_N + PC1_d1_V + PC1_d2_V
                PC2 =~ PC2_d1_N + PC2_d2_N + PC2_d1_V + PC2_d2_V
                PC1 ~~ PC2

                PC1 ~ PC1_q + PC2_q
                PC2 ~ PC2_q + PC1_q + cc_z

                PC1_q ~~ PC2_q

  '
fit_lf = sem(lf_mod, data=sem_dat) 
summary(fit_lf,fit.measures=TRUE,  standardized=T, rsquare=T)

sem_dat$PC1_lf = predict(fit_lf)[,1]
sem_dat$PC2_lf = predict(fit_lf)[,2]

exp_mod_noise = lm( noise_z ~ (PC1_lf + PC2_lf) * ss_z + pta_z , sem_dat  )
summary(exp_mod_noise)

exp_mod_voc = lm( voc_z ~ (PC1_lf + PC2_lf) * ss_z + pta_z , sem_dat  )
summary(exp_mod_voc)
