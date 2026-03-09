#Install packages
install.packages("vegan") 
library(vegan) 
install.packages("car") 
library(car) 
install.packages("lme4") 
library(lme4) 
install.packages("lmerTest") 
library(lmerTest) 
install.packages("emmeans") 
library(emmeans) 
install.packages("nlme") 
library(nlme) 
install.packages("dunn.test") 
library(dunn.test) 
install.packages("ggplot2") 
library(ggplot2) 
install.packages("reshape2") 
library(reshape2) 
install.packages("gridExtra") 
library(gridExtra)
install.packages("devtools") 
library(devtools) 
install.packages("adonis2") 
library(adonis2) 
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") 
library(pairwiseAdonis) 
library(multcomp)
library(multcompView)
install.packages("ARTool")
library(ARTool)
library(multcomp)
library(multcompView)

AIC(bacillota.m1, bacillota.m2, bacillota.m3)

bacillota.m1 <- lmer(Bacillota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
bacillota.m2 <- lmer(Bacillota ~ Se*Ye + (1|Room), vegan)
bacillota.m3 <- lmer(Bacillota ~ Se*Ye + (1|Tank), vegan)



vegan <- read.csv('VeganProportionsSeYe.csv', sep=',', header=TRUE, check.names=F)
vegan
vegan$Tank <- as.factor(vegan$Tank)
vegan$Diet <- as.factor(vegan$Diet)
vegan$Room <- as.factor(vegan$Room)
vegan$Se <- as.factor(vegan$Se)
vegan$Ye <- as.factor(vegan$Ye)

class(vegan)

#Phylum
bacillota.m <- lmer(Bacillota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.bacillota <- art(Bacillota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.bacillota
anova(bacillota.m)
bacillota.e <- emmeans(bacillota.m, pairwise ~ Se*Ye, type = "response", adjust = "dunn")
bacillota.e

pseudomonadota.m <- lmer(Pseudomonadota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.pseudomonadota <- art(Pseudomonadota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.pseudomonadota
anova(art.pseudomonadota)
pseudomonadota.e <- emmeans(pseudomonadota.m, pairwise ~ Se*Ye, type = "response", adjust = "dunn")
pseudomonadota.e

actinomycetota.m <- lmer(Actinomycetota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.actinomycetota <- art(Actinomycetota ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.actinomycetota
anova(art.actinomycetota)
actinomycetota.e <- emmeans(actinomycetota.m, pairwise ~ Se*Ye, type = "response", adjust = "dunn")
actinomycetota.e

#Genera
aeromonas.m <- lmer(Aeromonas ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.aeromonas <- art(Aeromonas ~ Se*Ye + (1|Tank) + (1|Room), vegan)
art.aeromonas
anova(art.aeromonas)
aeromonas.e <- emmeans(aeromonas.m, pairwise ~ Se|Ye, type = "response", adjust = "dunn")
aeromonas.e

streptococcus.m <- lmer(Streptococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.streptococcus <- art(Streptococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.streptococcus
anova(art.streptococcus)
streptococcus.e <- emmeans(streptococcus.m, pairwise ~ Se|Ye, type = "response", adjust = "dunn")
streptococcus.e

peptostreptococcus.m <- lmer(Peptostreptococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.peptostreptococcus <- art(Peptostreptococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.peptostreptococcus
anova(art.peptostreptococcus)
peptostreptococcus.e <- emmeans(peptostreptococcus.m, pairwise ~ Ye|Se, type = "response", adjust = "dunn")
peptostreptococcus.e


romboutsia.m <- lmer(Romboutsia ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.romboutsia <- art(Romboutsia ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.romboutsia
anova(art.romboutsia)
romboutsia.e <- emmeans(romboutsia.m, pairwise ~ Se|Ye, type = "response", adjust = "dunn")
romboutsia.e
romboutsia.cld <- cld(object = romboutsia.e, adjust = "dunn", Letters = letters, alpha = 0.05)
romboutsia.cld

mammaliicoccus.m <- lmer(Mammaliicoccus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.mammaliicoccus <- art(Mammaliicoccus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.mammaliicoccus
anova(art.mammaliicoccus)
mammaliicoccus.e <- emmeans(mammaliicoccus.m, pairwise ~ Se|Ye, type = "response", adjust = "dunn")
mammaliicoccus.e
mammaliicoccus.cld <- cld(object = mammaliicoccus.e, adjust = "dunn", Letters = letters, alpha = 0.05)
mammaliicoccus.cld

carnobacterium.m <- lmer(Carnobacterium ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.carnobacterium <- art(Carnobacterium ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.carnobacterium
anova(art.carnobacterium)
  carnobacterium.e <- emmeans(carnobacterium.m, pairwise ~ Se|Ye, type = "response", adjust = "dunn")
carnobacterium.e
carnobacterium.cld <- cld(object = carnobacterium.e, adjust = "dunn", Letters = letters, alpha = 0.05)
carnobacterium.cld

vagococcus.m <- lmer(Vagococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.vagococcus <- art(Vagococcus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.vagococcus
anova(art.vagococcus)
vagococcus.e <- emmeans(vagococcus.m, pairwise ~ Se*Ye, type = "response", adjust = "dunn")
vagococcus.e
vagococcus.cld <- cld(object = vagococcus.e, adjust = "dunn", Letters = letters, alpha = 0.05)
vagococcus.cld

peptoniphilus.m <- lmer(Peptoniphilus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.peptoniphilus <- art(Peptoniphilus ~ Se*Ye + (1|Room) + (1|Tank), vegan)
art.peptoniphilus
anova(art.peptoniphilus)
peptoniphilus.e <- emmeans(peptoniphilus.m, pairwise ~ Ye|Se, type = "response", adjust = "dunn")
peptoniphilus.e
peptoniphilus.cld <- cld(object = peptoniphilus.e, adjust = "dunn", Letters = letters, alpha = 0.05)
peptoniphilus.cld
 