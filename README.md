### CCJS 710 Lesson 9 - Thursday 10/27/22

* We now turn our attention to the problem of analyzing counted outcome data. Counted data can generally be grouped into 3 categories: (1) bounded counts; (2) unbounded counts; and (3) censored counts. We will discuss each of these categories.
* For today's class, we will consider the problem of bounded count data. Here is a reading that addresses the issue ([link](https://link.springer.com/article/10.1007/s10940-017-9346-9)).
* Let's begin by reading in a dataset (emailed to you) with 1000 cases and 2 variables (x and y):

```R
df <- read.csv(file="df.csv",header=T,sep=",")
head(df)
tail(df)
```

* Here is the output:

```Rout
> df <- read.csv(file="df.csv",header=T,sep=",")
> head(df)
  X          x y
1 1 -1.1155665 3
2 2 -0.8795643 5
3 3 -0.7568539 2
4 4 -0.2033934 4
5 5 -1.1694296 0
6 6  0.2639083 8
> tail(df)
        X           x y
995   995 -1.73083611 2
996   996  0.03160026 9
997   997  0.18312608 6
998   998 -0.51101843 7
999   999  0.53383795 7
1000 1000 -0.44075654 5
> 
```

* Next, we create a table showing the distribution of *df$y*:

```R
table(df$y,exclude=NULL)
```

which gives us:


```Rout
> table(df$y,exclude=NULL)

  0   1   2   3   4   5   6   7   8   9  10  11  12 
 17  40  54 109 103 124 100 121 108  88  74  43  19 
> 
```

* Notice that these counts are bounded (range of 0-12).
* What kinds of variables might take this form? (discussion)
* Here is a scatterplot of *df$x* and *df$y*:

```R
plot(x=df$x,y=df$y)
```

which gives us this plot:

<p align="left">
<img src="/gfiles/xy-scatterplot.png" width="600px">
</p>

* Now, let's construct a statistical model of this relationship.

```
m1 <- glm(cbind(y,12-y)~1+x,data=df,family=binomial(link="logit"))
summary(m1)
logLik(m1)
```

* Here is the output:


```Rout
> m1 <- glm(cbind(y,12-y)~1+x,data=df,family=binomial(link="logit"))
> summary(m1)

Call:
glm(formula = cbind(y, 12 - y) ~ 1 + x, family = binomial(link = "logit"), 
    data = df)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-3.0461  -0.6559  -0.0394   0.6634   3.9342  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept) 0.002481   0.019914   0.125    0.901    
x           0.998452   0.024688  40.443   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 3139.9  on 999  degrees of freedom
Residual deviance: 1097.8  on 998  degrees of freedom
AIC: 3731.8

Number of Fisher Scoring iterations: 4

> logLik(m1)
'log Lik.' -1863.885 (df=2)
> 
```

* So, we can see that we have an intercept estimate of approximately zero and a regression coefficient estimate of approximately 1. How do we interpret these estimates?

* Let's further consider this question by writing down the likelihood function:

```R
library(maxLik)

ll2 <- function(parms)
  {
   a <- parms[1]
   b <- parms[2]
   p <- exp(a+b*df$x)/(1+exp(a+b*df$x))
   pt1 <- factorial(12)
   pt2 <- factorial(df$y)*factorial(12-df$y)
   pt3 <- p^df$y
   pt4 <- (1-p)^(12-df$y)
   pmf <- pt1/pt2*pt3*pt4
   lpmf <- log(pmf)
   return(lpmf)
  }

m2 <- maxLik(ll2,start=c(0.02384042,1.02384023),
             method="BHHH",finalHessian="BHHH")
summary(m2)
```

which yields the following results:

```rout
> library(maxLik)
> 
> ll2 <- function(parms)
+   {
+    a <- parms[1]
+    b <- parms[2]
+    p <- exp(a+b*df$x)/(1+exp(a+b*df$x))
+    pt1 <- factorial(12)
+    pt2 <- factorial(df$y)*factorial(12-df$y)
+    pt3 <- p^df$y
+    pt4 <- (1-p)^(12-df$y)
+    pmf <- pt1/pt2*pt3*pt4
+    lpmf <- log(pmf)
+    return(lpmf)
+   }
> 
> m2 <- maxLik(ll2,start=c(0.02384042,1.02384023),
+              method="BHHH",finalHessian="BHHH")
> summary(m2)
--------------------------------------------
Maximum Likelihood estimation
BHHH maximisation, 3 iterations
Return code 8: successive function values within relative tolerance limit (reltol)
Log-Likelihood: -1863.885 
2  free parameters
Estimates:
     Estimate Std. error t value Pr(> t)    
[1,] 0.002481   0.019715   0.126     0.9    
[2,] 0.998451   0.025243  39.553  <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
--------------------------------------------
> 
```

* Next, we extract the parameter estimates form this model:

```R
a <- coef(m2)[1]
a
b <- coef(m2)[2]
b
```

which we then use to plot the functional form of the relationship between *df$x* and *df$y*:

```R
p <- exp(a+b*df$x)/(1+exp(a+b*df$x))
plot(x=df$x,y=p)
```

which yields:

<p align="left">
<img src="/gfiles/fform.png" width="600px">
</p>

* Next, we check on the model's fit. First, we calculate the expected frequencies:

```R
# calculate expected cell frequencies

e0  <- sum(choose(12,0)*p^0*(1-p)^(12-0))
e1  <- sum(choose(12,1)*p^1*(1-p)^(12-1))
e2  <- sum(choose(12,2)*p^2*(1-p)^(12-2))
e3  <- sum(choose(12,3)*p^3*(1-p)^(12-3))
e4  <- sum(choose(12,4)*p^4*(1-p)^(12-4))
e5  <- sum(choose(12,5)*p^5*(1-p)^(12-5))
e6  <- sum(choose(12,6)*p^6*(1-p)^(12-6))
e7  <- sum(choose(12,7)*p^7*(1-p)^(12-7))
e8  <- sum(choose(12,8)*p^8*(1-p)^(12-8))
e9  <- sum(choose(12,9)*p^9*(1-p)^(12-9))
e10 <- sum(choose(12,10)*p^10*(1-p)^(12-10))
e11 <- sum(choose(12,11)*p^11*(1-p)^(12-11))
e12 <- sum(choose(12,12)*p^12*(1-p)^(12-12))

# collect expected frequencies into a vector

e <- c(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12)
e
```

which yields:

```Rout
> e <- c(e0,e1,e2,e3,e4,e5,e6,e7,e8,e9,e10,e11,e12)
> e
 [1]  15.42094  40.75232  67.83955  90.60545 106.18798 114.69966
 [7] 117.21441 114.38691 106.40160  92.89133  72.51787  44.88954
[13]  16.19243
> 
```

* Then, we collect the observed frequencies:

```R
# collect observed frequencies

table(df$y)

o0  <- table(df$y)[1]
o1  <- table(df$y)[2]
o2  <- table(df$y)[3]
o3  <- table(df$y)[4]
o4  <- table(df$y)[5]
o5  <- table(df$y)[6]
o6  <- table(df$y)[7]
o7  <- table(df$y)[8]
o8  <- table(df$y)[9]
o9  <- table(df$y)[10]
o10 <- table(df$y)[11]
o11 <- table(df$y)[12]
o12 <- table(df$y)[13]

# put the observed frequencies into a vector

o <- c(o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12)
o

```

which yields:

```Rout
> table(df$y)

  0   1   2   3   4   5   6   7   8   9  10  11  12 
 17  40  54 109 103 124 100 121 108  88  74  43  19 
> 
> o0  <- table(df$y)[1]
> o1  <- table(df$y)[2]
> o2  <- table(df$y)[3]
> o3  <- table(df$y)[4]
> o4  <- table(df$y)[5]
> o5  <- table(df$y)[6]
> o6  <- table(df$y)[7]
> o7  <- table(df$y)[8]
> o8  <- table(df$y)[9]
> o9  <- table(df$y)[10]
> o10 <- table(df$y)[11]
> o11 <- table(df$y)[12]
> o12 <- table(df$y)[13]
> 
> # put the observed frequencies into a vector
> 
> o <- c(o0,o1,o2,o3,o4,o5,o6,o7,o8,o9,o10,o11,o12)
> o
  0   1   2   3   4   5   6   7   8   9  10  11  12 
 17  40  54 109 103 124 100 121 108  88  74  43  19 
> 
```

* So, the test for model fit is:

```R
data.frame(o,e,o-e)
chi.sq <- (o-e)^2/e
sum(chi.sq)
dof <- 13-2
dof
pval <- 1-pchisq(sum(chi.sq),dof)
pval
```

which gives:

```Rout
> data.frame(o,e,o-e)
     o         e       o...e
0   17  15.42094   1.5790594
1   40  40.75232  -0.7523185
2   54  67.83955 -13.8395493
3  109  90.60545  18.3945451
4  103 106.18798  -3.1879848
5  124 114.69966   9.3003380
6  100 117.21441 -17.2144075
7  121 114.38691   6.6130861
8  108 106.40160   1.5984039
9   88  92.89133  -4.8913307
10  74  72.51787   1.4821315
11  43  44.88954  -1.8895428
12  19  16.19243   2.8075695
> chi.sq <- (o-e)^2/e
> sum(chi.sq)
[1] 11.37183
> dof <- 13-2
> dof
[1] 11
> pval <- 1-pchisq(sum(chi.sq),dof)
> pval
[1] 0.4126589
> 
```

* Now, we replot the joint distribution of *df$x* and *df$y*:

```R
plot(x=df$x,y=df$y)
```

and we now add the expected distribution of *y* conditional on different values of *x* to the plotspace:

```R
# calculate E(y|x=xstar) for range of xstar values

eyx <- vector()
xstar <- vector()

for (i in seq(from=1,to=1000,by=1))
  {
    # rescale 1 to 1000 --> -4 to 4

    lower <- -4
    upper <- 4
    min.i <- 1
    max.i <- 1000
    xstar[i] <- ((upper-lower)*(i-min.i)/(max.i-min.i))+lower

    pstar <- exp(a+b*xstar[i])/(1+exp(a+b*xstar[i]))
    p0  <- choose(12,0)*pstar^0*(1-pstar)^(12-0)
    p1  <- choose(12,1)*pstar^1*(1-pstar)^(12-1)
    p2  <- choose(12,2)*pstar^2*(1-pstar)^(12-2)
    p3  <- choose(12,3)*pstar^3*(1-pstar)^(12-3)
    p4  <- choose(12,4)*pstar^4*(1-pstar)^(12-4)
    p5  <- choose(12,5)*pstar^5*(1-pstar)^(12-5)
    p6  <- choose(12,6)*pstar^6*(1-pstar)^(12-6)
    p7  <- choose(12,7)*pstar^7*(1-pstar)^(12-7)
    p8  <- choose(12,8)*pstar^8*(1-pstar)^(12-8)
    p9  <- choose(12,9)*pstar^9*(1-pstar)^(12-9)
    p10 <- choose(12,10)*pstar^10*(1-pstar)^(12-10)
    p11 <- choose(12,11)*pstar^11*(1-pstar)^(12-11)
    p12 <- choose(12,12)*pstar^12*(1-pstar)^(12-12)
    eyx[i] <- p0*0+p1*1+p2*2+p3*3+p4*4+p5*5+p6*6+
              p7*7+p8*8+p9*9+p10*10+p11*11+p12*12
 }

# annotate the plot with a curve showing the relationship
# between xstar and E(y|x=xstar) implied by the statistical model

lines(x=xstar,y=eyx,pch=19,lty=1,lwd=3,col="blue")
```

and the resulting plot is:

<p align="left">
<img src="/gfiles/eyx-plot.png" width="600px">
</p>
