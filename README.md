---
output:
  pdf_document: default
  html_document: default
---

# Reproducible research: version control and R

The answers for Questions 1, 2 and 3 can be found in my `logistic_growth` repo, found using this link:

<https://github.com/anonymous-student1/logistic_growth>

### Question 4

```{r load packages required for the following analysis}
    
install.packages(c("ggplot2", "tidyverse", "SemiPar", "gridExtra", "dplyr"))
    
library(c(ggplot2, tidyverse, SemiPar, gridExtra, dplyr))
```

1.  

The image address of this simulated walk is: <https://f692516e24f04defb196d930e13d7f70.app.posit.cloud/graphics/a55500bc-0980-45f8-b58e-ed3b9549a33e.png>

The output of the code produces two random walks that are very different. These walks both have a total of 500 steps, which are colour coded using the same colour gradient to indicate the start and finish. This is produced by the function colour = time, and so the same numbered steps in both graphs are represented by the same shade in the blue gradient. The origin coordinate of these are both at (0.0) at time = 0, and progress from there to produce different walks. The angle at which the next step can be taken can be any angle from the previous point, with the largest change that can occur in both the x or y axes being a coordinate change of 0.25, causing an step to be taken either vertically or horizontally (with no change in the other axis due to the interaction between sin and cos in the n_steps function).

The walk in the left hand simulation does show are more compact walk, spanning across 4 X coordinates and 6 Y coordinates, while the right hand walk spans 7 X coordinates and 7 Y coordinates. The right hand walk covered a larger areas of space on the axis, but this is due to random sampling of the angles between 0 and 2\*pi, as opposed to anything else here.

A long straight line appears when the same angle is randomly sampled twice in a row, as is seen twice in the last 20 steps of the left hand walk, but this occurs infrequently.

There is no clear order or pattern observed in the direction or order of steps taken in either simulated walk, however both walks here do begin in with positive x values and negative y values, with the left hand walk later moving into positive values of y and the right hand walk moving further into the negative values of y.

There are lots of time points at which the random steps are bunched together in small areas of the axis, and a few where there is some direction to steps all leading the same way, however there are only a few of these areas, that appear as 'lines' away from the clusters, likely due to the decreased probability of randomly sampling a similar angle enough times consecutively in order to see this effect.

2.  

A random seed is described as a number or vector that is used to start up a pseudorandom number generator. This mimics the properties of the formation of a uniform distribution and can be used to draw numbers at random to use in a model. To generate random numbers, a recursive algorithm (a Random Number Generator) is required and the seed is the starting state of this algorithm.

3.  

```{r edit script}

random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
  
  set.seed(1)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}

data1 <- random_walk(500)

plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

data2 <- random_walk(500)

plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

grid.arrange(plot1, plot2, ncol=2)

```

4.  This has been committed to the forked repo, with this image showing the change made
   <img width="580" alt="Screenshot 2023-12-05 at 15 12 17" src="https://github.com/anonymous-student1/reproducible-research_homework/assets/150151047/d390dc2e-99f1-47f8-9229-c6c57e7e25d3">


### Question 5

1.  

```{r}

    file.choose() 
        
    virus_data\<-read.csv("/cloud/project/question-5-data/Cui_etal2014.csv") 
        
    summary(virus_data)
```

The summary of the virus data shows that there are 33 rows and 13 columns of data in the dsDNA viruses dataset.

2.  

```{r plots to visualise}

library(ggplot2)

ggplot(data = virus_clean, aes(x = virion_volume))+
  geom_histogram()

ggplot(data = virus_clean, aes(x = genome_length))+
  geom_histogram()
  
```

Plotting these measurements on a histogram shows that these both have a significant right skew, with data spanning over several orders of magnitude, so the most appropriate transformation to use to fit to a linear model is a logarithm transformation.

```{r apply transformation}
virus_clean$Virion.volume..nm.nm.nm.<- as.numeric(virus_clean$Virion.volume..nm.nm.nm.)
  
colnames(virus_clean)[10] = "virion_volume"
colnames(virus_clean)[12] = "genome_length"

virus_clean <- mutate(virus_clean, log_virion_volume = log(virion_volume))
virus_clean <- mutate(virus_clean, log_genome_length = log(genome_length))
```

This has applied the transformation to the data for virion volume and genome length, and has stored these in new columns added to the dataset being worked on (this has been created in order not to add columns to the raw data set)

3\.

Applying the log transformation to the equation V = $\beta$L$\alpha$ produces the equation: logV = log$\beta$ + $\alpha$logL.

This takes the form of a straight line with the formula y = ax + b, and so a linear model can be used to deduce the values of $\alpha$ and $\beta$.

```{r finding values of alpha and beta}

virus_model<- lm(formula = log_virion_volume ~ log_genome_length, data = virus_clean)
coef(virus_model)
summary(virus_model)
```

The output of the linear models produces a value of $\alpha$, the exponent, to be 1.5152, with a p value of 6.44e10, a statistically signficant p value. When this is rounded to 3 significant figures, 1.52, this matches the allometric exponent value by Cui et al., 2014 exactly.

The value output for the scaling value is 7.0748, with a p value of 2.28e-10, which is also statistically signficant, but this represents the value of log$\beta$, and so the value of $\beta$ can be found by exp(7.0748) = 1181.807. Cui et al., 2014 state that the value of the scaling component, $\beta$, is 1182, which again matches the paper once rounded.

4.  

```{r figure}

ggplot(data = virus_clean, aes(x = log_genome_length, y = log_virion_volume))+
  geom_point()+
  geom_smooth(method = "lm")+
  xlab("log [Genome Length (kb)]") +
  ylab("log [Virion volume (nm3)]") +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold"))+
  theme(axis.title.y = element_text(face = "bold"))
```

5.  

```{r estimate volume when genome = 300kb}
b<-1181.807.   
L<-300000.   
a<-1.5152.   
    
V <- b*(L^a)
V
```

V is estimated to be 235223486693nm^3^

sink(file = "package-versions.txt")
sessionInfo()
sink()
