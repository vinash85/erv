## generate a single cell samples Expression(ij) = Treatment(i) + Patient(j) + eij.

patient = rep(sample(10),each=1000)
patient.exp = rnorm(n = 10000, mean=patient, sd=1)
treatment = rep(c(0,1), each=5000)
treatment.exp = rnorm(n = 10000, mean=treatment, sd=1)
noise.exp = rnorm(10000)
exp = treatment.exp + 1/5* patient.exp + 1/5* noise.exp

library(nlme)
extract_nmle_table <- function (m1){
    mod = summary(m1)
    beta <- m1$coefficients$fixed[2] #$fixed is not needed
    se <- m1$varFix[2]
    t <- beta/se
    p<- anova(m1)[2,4]
    table=data.frame(cbind(beta,se,t,p))
    return(table)
}
eval.nlme = function(gene.expression, data.dt){
    # data.dt = data.table with colums  Response1 and dataset (factor to control)
    tryCatch(
        {
            data.dt$col = gene.expression
            data.dt = data.dt[!is.na(col)]
            data.dt = data.dt[dataset %in% names(which(table(data.dt$dataset) > 10))]  # only choose those with at least 10 obervations 
            m1 <- lme(col~ Response1, random=~1|dataset, data=data.dt)
            list(m1, extract_nmle_table(m1))
        },error=  function(e) rep(NA,4))
}

data.dt= dt1 = data.table(Response1=treatment, dataset=factor(patient))
aa = eval.nlme(gene.expression=exp, data.dt =dt1)


data.dt = dt1
library(lme4) 
data.dt$col = exp
cntrl <- lmer(col~ 1 + 1|dataset, data=data.dt[Response1==0])
treated1 <- lmer((col)~ 1 + 1|dataset, data=data.dt[Response1==1])
treated <- lmer((col-1)~ 1 + 1|dataset, data=data.dt[Response1==1])

combined1 = lm(exp~treatment+patient)
combined = lm(exp~treatment)
combined = lm(exp~patient)

order.inx = rep(1:10, each=1000)
treatment.mcmc = rep(NA,100)
patient.effect.mcmc = matrix(NA, nrow=100, ncol=10)
treatment.curr = 0
for (ii in seq(100)) {
	data.dt[,exp.res.treatment:=col - Response1*treatment.curr]
	cntrl <- lmer(exp.res.treatment~ 1 + 1|dataset, data=data.dt[Response1==0])
	treated <- lmer(exp.res.treatment ~ 1 + 1|dataset, data=data.dt[Response1==1])

	patient.effect.mcmc[ii,] = coef.order = unlist(c(coef(cntrl)$dataset, coef(treated)$dataset))
	patient.effect = coef.order[order.inx]
	data.dt[,exp.res.patient:=col - patient.effect]

	combined = lm(exp.res.patient~Response1, data=data.dt)
	treatment.curr  = treatment.mcmc[ii] = coef(combined)["Response1"]
	print(ii)
}

combined = lm(col~Response1, data=data.dt)
patient.order  =unique(patient)
data.dt= dt1 = data.table(Response1=treatment, dataset=factor(patient), col = exp)
order.inx = rep(1:10, each=1000)
treatment.mcmc = rep(NA,100)
patient.effect.mcmc = matrix(NA, nrow=100, ncol=10)
treatment.curr = 0
data.dt[,exp.res.patient:=col]
for (ii in seq(100)) {
	combined = lm(exp.res.patient~Response1, data=data.dt)
	treatment.curr  = treatment.mcmc[ii] = coef(combined)["Response1"]
	
	data.dt$exp.res.treatment= data.dt$col - data.dt$Response1*treatment.curr
	# cntrl <- lmer(exp.res.treatment~ 1 + 1|dataset, data=data.dt[Response1==0])
	# treated <- lmer(exp.res.treatment ~ 1 + 1|dataset, data=data.dt[Response1==1])
	treated <- lmer(exp.res.treatment ~ 1|dataset, data=data.dt)

	# patient.effect.mcmc[ii,] = coef.order = unlist(c(coef(cntrl)$dataset, coef(treated)$dataset))
	patient.effect.mcmc[ii,] = coef.order = unlist( coef(treated)$dataset[patient.order,])
	patient.effect = coef.order[order.inx]
	data.dt$exp.res.patient= data.dt$col - patient.effect
	# data.dt[,exp.res.patient:=col - patient.effect]
	print(ii)
}