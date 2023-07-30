library(data.table)
source("MMM_code.R")
dat = fread("D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/data.csv")
missing = rep(F,nrow(dat))
for (i in 1:(ncol(dat)-1)) {
  missing.tmp = dat[,..i]==-1
  missing = missing | missing.tmp
}
question_names = fread("./WPI/question_names.csv")
ques_names = question_names$ques_names

missing = as.vector(missing)
dat.new = dat[!missing]
# clean age misspecification >150
dat.new = dat.new[dat.new$age<100]
A = dat.new[,1:116]
A[A==2] = 0
X = dat.new[,117:118]
age.scale = scale(X[["age"]])
X[,1] = age.scale
X[X$gender==3,"gender"] = 0
X[X$gender==2,"gender"] = -1
X = as.matrix(X)
A = as.matrix(A)

K = 3
R = ncol(X)

# jml
res.jml = gom.jml.self(A,K)
fwrite(res.jml.sirt$g,"./WPI/Pi_jml.csv")
fwrite(res.jml.sirt$lambda[,3:5],"./WPI/Theta_jml.csv")
sum(apply(res.jml$Theta, 1, var))
res.jml.sirt = gom.jml(data.frame(A),K)
rowsums(res.jml.sirt$g)
sum(apply(res.jml.sirt$lambda[,3:5],1,var))


transpose = Rfast::transpose
AA = mat.mult(A,transpose(A))
AA = Diag.fill(AA,0)
XX = mat.mult(X,transpose(X))
XX = Diag.fill(XX,0)
AA.evd = eigs(AA,K+1)
XX.evd = eigs(XX,K+1)
alpha_min = (AA.evd$values[K]-AA.evd$values[K+1])/XX.evd$values[1]
if (R<K){
  # alpha_max = (AA.evd$values[1])/(XX.evd$values[R])
  alpha_max = AA.evd$values[K]/XX.evd$values[R]
} else {
  alpha_max = (AA.evd$values[K])/(XX.evd$values[K]-XX.evd$values[K+1])
}

# alpha_max = (AA.evd$values[K]-AA.evd$values[K+1])/XX.evd$values[R]
alpha_seq = seq(from=alpha_min,to=alpha_max,length=100)# fail sometimes
alpha_seq1 = seq(from=0,to=alpha_max,length=100)

err = future_sapply(alpha_seq1,find_best_alpha,A,X,K,future.seed=T)
ggplot()+
  geom_line(aes(alpha_seq1,err))+
  geom_smooth(aes(alpha_seq1,err))+
  xlab("alpha")
# alpha = 0.23,0.4,0.6
res1 = cov_cluster(A,X,K,.23)
res.null = null_cluster(A,K)
res2 = cov_cluster(A,X,K,0.6)
res3 = cov_cluster(A,X,K,0.4)
mean(abs(res1$Pi-res2$Pi))# about 0.13 QUITE much
mean(abs(res1$Pi-res3$Pi))# 0.14 quite much
mean(abs(res2$Pi-res3$Pi))# 0.044 quite little
head(res1$Pi)
head(res.null$Pi)
mean(abs(res2$Pi-res.null$Pi)) # 0.11 not quite much
fwrite(res2$Theta,"D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/cov_theta.csv")
fwrite(res1$Theta,"D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/cov_theta1.csv")
fwrite(res.null$Theta,"D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/null_theta.csv")

Theta = fread("D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/cov_theta.csv")
Theta.null = fread("./WPI/null_theta.csv")
sum(apply(Theta.null, 1, var))
Theta.cov1 = fread("./WPI/cov_theta.csv")
sum(apply(Theta.cov1, 1, var))
Theta.cov2 = fread("./WPI/cov_theta1.csv")
sum(apply(Theta.cov2, 1, var))
idx1 = order(Theta$V1,decreasing = T)[1:10]
idx2 = order(Theta$V2,decreasing = T)[1:10]
idx3 = order(Theta$V3,decreasing = T)[1:10]
idx = c(idx1,idx2,idx3)
idx = unique(idx)
idx = sort(idx)
# idx = order(apply(Theta,1,var),decreasing = T)[1:30]
idx = 31:60
Theta.subset = Theta[idx,]
rownames(Theta.subset) = ques_names[idx]
colnames(Theta.subset) = c("Profile 1", "Profile 2", "Profile 3")
library(tidyverse)
Theta.subset2 = Theta.subset %>% 
  # as_tibble() %>% 
  rownames_to_column("Questions") %>%
  pivot_longer(-Questions, names_to = "Profile", values_to = "value") %>% 
  mutate(
    Profile = factor(Profile),
    Questions = factor(Questions,levels=unique(Questions))
  )
ggplot(Theta.subset2,aes(Profile,Questions))+
  geom_raster(aes(fill=value))+ 
  geom_text(aes(label=round(value,3))) + 
  theme(axis.ticks = element_blank(), panel.background = element_blank(),legend.position = "none") +
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette = "Reds", direction = 1)


library(ggtern)
pi.dat = data.frame(res1$Pi)
pi.dat$Age = dat.new$age
pi.dat$Gender = dat.new$gender
pi.dat$Gender[dat.new$gender==1] = "Male"
pi.dat$Gender[dat.new$gender==2] = "Female"
pi.dat$Gender[dat.new$gender==3] = "Other"
pi.dat$Gender = as.factor(pi.dat$Gender)
ggtern(data = pi.dat, aes(X1,X2,X3)) + 
  geom_point(aes(col=Age))+
  theme_hidegrid()+
  theme_hideticks()+
  theme_hidelabels()+
  xlab("Profile 1")+
  ylab("Profile 2")+
  zlab("Profile 3")+
  scale_color_distiller(palette = "YlOrRd")
ggtern(data = pi.dat, aes(X1,X2,X3)) + 
  geom_point(aes(col=Gender))+
  theme_hidegrid()+
  theme_hideticks()+
  theme_hidelabels()+
  xlab("Profile 1")+
  ylab("Profile 2")+
  zlab("Profile 3")+
  scale_color_brewer(palette = "Set2")

questions = c(35, 68, 94, 91, 37, 96, 19, 97, 34, 76,
              116, 69, 33, 8, 39, 95, 93, 46, 113, 25,
              114, 40, 73, 115, 45, 105, 106, 10, 90, 85)
names_ques = c("Were you shy with other boys [or girls]?",
               "Do you worry too much about little things?",
               "Do you get tired of people quickly?",
               "Is it easy to make you laugh?",
               "Did you ever have a strong desire to run away from home?",
               "Do you get tired of work?",
               "Have you ever had fits of dizziness?",
               "Do your interests change frequently?",
               "Did the other children let you play with them?",
               "Can you do the little chores of the day without worrying over them?",
               "Do you like outdoor life?",
               "Do you think you worry too much when you have an unfinished job on your hands?",
               "As a child did you like to play alone better than to play with other children?",
               "Do you have the sensation of falling when going to sleep?",
               "Did the teachers in school generally treat you right?",
               "Do you get tired of amusements quickly?",
               "Is it easy to get you cross or grouchy?",
               "Do you find your way about easily?",
               "Can you stand the sight of blood?",
               "Have you ever fainted away?",
               "Can you stand pain quietly?",
               "Have your employers generally treated you right?",
               "Can you sit still without fidgeting?",
               "Can you stand disgusting smells?",
               "Do you get used to new places quickly?",
               "Did you ever have dyspepsia [indigestion]?",
               "Did you ever have asthma or hay fever [allergies]?",
               "Do ideas run through your head so that you cannot sleep?",
               "Have you a good appetite?",
               "Did you ever have the habit of biting your finger nails?")
questions = rev(questions)
names_ques = rev(names_ques)
Theta = fread("D:/Programmes/Covariate-Assisted-Spectral-GoM/WPI/cov_theta1.csv")
Theta.subset = Theta[questions,]
rownames(Theta.subset) = names_ques
colnames(Theta.subset) = c("Profile 1", "Profile 2", "Profile 3")
library(tidyverse)
Theta.subset2 = Theta.subset %>% 
  # as_tibble() %>% 
  rownames_to_column("Questions") %>%
  pivot_longer(-Questions, names_to = "Profile", values_to = "value") %>% 
  mutate(
    Profile = factor(Profile),
    Questions = factor(Questions,levels=unique(Questions))
  )
# find_row_names = function(i){
#   return(names_ques[i])
# }
# Theta.subset2$Questions = sapply(Theta.subset2$Questions,find_row_names)
  
ggplot(Theta.subset2,aes(Profile,Questions))+
  geom_raster(aes(fill=value))+ 
  geom_text(aes(label=round(value,3))) + 
  theme(axis.ticks = element_blank(), panel.background = element_blank(),legend.position = "none") +
  xlab("")+
  ylab("")+
  scale_fill_distiller(palette = "Reds", direction = 1)

# 2 seems more interpretable to socially passive, unhealthy and normal
