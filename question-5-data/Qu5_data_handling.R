library(ggplot2)
install.packages(SemiPar)
library(SemiPar)


#Question 2
ggplot(data = virus_clean, aes(x = virion_volume))+
  geom_histogram()

ggplot(data = virus_clean, aes(x = genome_length))+
  geom_histogram()

ggplot(data = virus_clean, aes(x = virion_volume, y = genome_length))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab

ggplot(data = virus_clean, aes(x = log_genome_length, y = log_virion_volume))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("log(Genome Length), kb") +
  ylab("log(Virion volume), nm^3")
  

virus_model<- lm(formula = log_virion_volume ~ log_genome_length, data = virus_clean)
coef(virus_model)
summary(virus_model)

exp(7.0748)
exp(1.5152)









