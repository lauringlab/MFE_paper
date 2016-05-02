# MFE distributions Maximum likelihood
JT McCrone  
April 22, 2016  




In this document we'll try to fit the distribution of the MFE using maximum likelihood methods and the r package bbmle. We'll also rely on Rs built in distributions and for the time being we'll be using the data in Ashley _et. al_ until the influenza data is ready.


```r
data.df <- read.csv("../data/flu.csv", stringsAsFactors = F)

ggplot(data.df, aes(x = Total.Hist)) + geom_histogram()
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](mll_fits_files/figure-html/the_data-1.png)

To start we'll use the models from Sanjuan 2004 and to make things even easier we'll just look at the simple noncompounded models first and save the compounded ones for later. If I understand the paper correctly they just model the distribution of the negative effects (<1.0 , not necessarily significant) that are not lethal. We could model the lethal effects with a zero-inflated model, but we'll save that for later.

I'll subset the data above to include the non-lethal negative fitnesses

```r
model.df <- mutate(data.df, Fitness = Total.CDF)  # change this to other columns if you'd like to model those as well.
model.df <- subset(model.df, Fitness > 0 & Fitness < 1)
ggplot(model.df, aes(x = Fitness)) + geom_histogram()
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](mll_fits_files/figure-html/unnamed-chunk-1-1.png)

### functions



```r
gammaNLL <- function(shape, scale, data) {
    # gamma negative log likelihood
    -sum(dgamma(data, shape = shape, scale = scale, log = T))  # the negative of the sum of the prob of seeing the data give the shape and scale parameters. Log transform the probabilities
}

gamma_fit <- function(fitness) {
    m = mean(fitness)
    vm = var(fitness)/mean(fitness)
    
    gammamodel <- mle2(gammaNLL, start = list(shape = m/vm, scale = vm), data = list(data = fitness))
    gammamodel
    AIC(gammamodel, k = 2)
}


betaNLL <- function(a, b, data) {
    # print(data)
    -sum(dbeta(data, shape1 = a, shape2 = b, log = T))
}

beta_fit <- function(fitness) {
    data <- fitness
    # print(data)
    betamodel <- mle2(betaNLL, start = list(a = 5, b = 5), data = list(data = fitness))
    betamodel
    AIC(betamodel, k = 2)
}
# plot(profile(betamodel))


WeibullNLL <- function(a, b, data) {
    -sum(dweibull(data, shape = a, scale = b, log = T))
    
}

weibull_fit <- function(fitness) {
    weibullmodel <- mle2(WeibullNLL, start = list(a = 5, b = 5), data = list(data = fitness))
    weibullmodel
    AIC(weibullmodel, k = 2)
}


lognormalNLL <- function(a, b, data) {
    -sum(dlnorm(data, meanlog = a, sdlog = b, log = T))
    
}

lnorm_fit <- function(fitness) {
    lnormmodel <- mle2(lognormalNLL, start = list(a = mean(log(fitness)), b = sd(log(fitness))), 
        data = list(data = fitness))
    lnormmodel
    AIC(lnormmodel, k = 2)
}


expNLL <- function(a, data) {
    -sum(dexp(data, rate = a, log = T))
}

exp_fit <- function(fitness) {
    expmodel <- mle2(expNLL, start = list(a = 1/mean(fitness)), data = list(data = fitness))
    expmodel
    AIC(expmodel, k = 1)
}
```


## Table


```r
make_table <- function(data) {
    fitness <- data[which(data < 1)]
    x <- data.frame(Distribution = c("Exponential", "Gamma", "Beta", "Weibull", 
        "Lognormal"), AIC = c(exp_fit(fitness), gamma_fit(fitness), beta_fit(fitness), 
        weibull_fit(fitness), lnorm_fit(fitness)))
}

total <- make_table(data.df$Total.CDF)

random <- make_table(data.df$Random.CDF)
surface <- make_table(data.df$Surface.CDF)
internal <- make_table(data.df$Internal.CDF)


all_table <- data.frame(Distribution = c("Exponential", "Gamma", "Beta", "Weibull", 
    "Lognormal"), Total = total$AIC, Random = random$AIC, Surface = surface$AIC, 
    Internal = internal$AIC)

knitr::kable(all_table)
```



Distribution         Total       Random     Surface     Internal
-------------  -----------  -----------  ----------  -----------
Exponential     109.850233    84.840199    49.11733    61.563198
Gamma            -2.575103     5.428507   -23.98510     9.812420
Beta            -68.684491   -48.754863   -42.91293   -28.769167
Weibull         -33.842623   -17.001904   -35.77080    -4.475457
Lognormal        12.485825    17.124983   -21.11320    18.134522




## Adding uniform distributions

In this analyis we add a uniform distribution. We ARE NOT letting the parameters vary but rather are letting a proportion of the variants fall on a uniform distribution between 0 and some parameter b.


###Exponential + uniform


```r
expUniNLL <- function(a, p, m) {
    expll <- p * dexp(model.df$Fitness, rate = a)
    unill <- (1 - p) * dunif(model.df$Fitness, 0, m)
    likeli <- sum(expll, unill, na.rm = T)
    LL = -sum(log(likeli))
    
    # print(c(a,p,m)) print(LL)
    
    if (is.finite(LL)) 
        return(LL) else return(1000)
}

expUnimodel <- mle2(expUniNLL, start = list(a = 1/mean(model.df$Fitness), p = 0.5, 
    m = 0.5), method = "L-BFGS-B", lower = c(0, 0, 0), upper = c(100, 1, 1))
expUnimodel
```

```
## 
## Call:
## mle2(minuslogl = expUniNLL, start = list(a = 1/mean(model.df$Fitness), 
##     p = 0.5, m = 0.5), method = "L-BFGS-B", lower = c(0, 0, 0), 
##     upper = c(100, 1, 1))
## 
## Coefficients:
##        a        p        m 
## 1.408008 1.000000 0.000000 
## 
## Log-likelihood: 3.58
```

```r
AIC(expUnimodel, k = 3)
```

```
## [1] 1.831134
```





### Gamma + uniform


```r
gammaUniNLL <- function(a, b, p, m) {
    # gamma negative log likelihood
    gammall <- p * dgamma(model.df$Fitness, shape = a, scale = b)
    unill <- (1 - p) * dunif(model.df$Fitness, 0, m)
    likeli <- sum(gammall, unill, na.rm = T)
    LL = -sum(log(likeli))
    # LL<-sum(log(p*dgamma(model.df$Fitness,shape=a,scale=b)+(1-p)*dunif(model.df$Fitness,0,m)))
    # # the negative of the sum of the prob of seeing the data give the shape
    # and scale parameters. Log transform the probabilities print(c(a,b,p,m))
    # print(LL)
    
    if (is.finite(LL)) 
        return(LL) else return(1000)
}

m = mean(model.df$Fitness)
vm = var(model.df$Fitness)/mean(model.df$Fitness)

gammaUniModel <- mle2(gammaUniNLL, start = list(a = m/vm, b = vm, p = 0.5, m = 0.5), 
    method = "L-BFGS-B", lower = list(p = 0, m = 0), upper = list(p = 1, m = 1))
gammaUniModel
```

```
## 
## Call:
## mle2(minuslogl = gammaUniNLL, start = list(a = m/vm, b = vm, 
##     p = 0.5, m = 0.5), method = "L-BFGS-B", lower = list(p = 0, 
##     m = 0), upper = list(p = 1, m = 1))
## 
## Coefficients:
##         a         b         p         m 
## 1.0000000 0.7102244 1.0000000 0.0000000 
## 
## Log-likelihood: 3.58
```

```r
AIC(gammaUniModel, k = 4)
```

```
## [1] 8.831134
```

 The extra distributions aren't really adding anything.






