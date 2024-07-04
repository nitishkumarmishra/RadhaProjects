## R code for t-test and two factor ANNOVA test
dat <- read.csv("Annova_File.txt", header=TRUE, sep="\t") # This is the same order as he give me

t.test(dat$Male[1:4], dat$Male[5:8])
t.test(dat$Female[1:4], dat$Female[5:8])

dat <- read.csv("ANNOVA_Modified_IL1A.txt", header=TRUE, sep="\t") #Here I have change the order
summary(aov(ILB1~as.factor(Gender),data=dat))
summary(aov(ILB1~as.factor(PF_Group),data=dat))
summary(aov(ILB1~as.factor(Gender)*as.factor(PF_Group),data=dat))
leveneTest(ILB1~as.factor(Group)*as.factor(Gender),data=dat)

