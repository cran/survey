# Bell-McAffrey standard errors run and are greater than linearized

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey"))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey")
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541 0.000175 ***
      api99         0.90329    0.02734  33.039 2.33e-12 ***
      stypeH      -19.38726    5.43114  -3.570 0.004398 ** 
      stypeM      -18.15821    6.07011  -2.991 0.012267 *  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

# Bell-McAffrey degrees of freedom go down

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey",
      degf = TRUE))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey", 
          degf = TRUE)
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541  0.00218 ** 
      api99         0.90329    0.02734  33.039 2.39e-07 ***
      stypeH      -19.38726    5.43114  -3.570  0.01452 *  
      stypeM      -18.15821    6.07011  -2.991  0.02824 *  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

---

    Code
      summary(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey-2",
      degf = TRUE))
    Output
      
      Call:
      svyglm(formula = api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey-2", 
          degf = TRUE)
      
      Survey design:
      dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
      
      Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
      (Intercept)  99.85645   18.02030   5.541  0.00989 ** 
      api99         0.90329    0.02734  33.039  3.8e-05 ***
      stypeH      -19.38726    5.43114  -3.570  0.03412 *  
      stypeM      -18.15821    6.07011  -2.991  0.05386 .  
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
      (Dispersion parameter for gaussian family taken to be 710.3237)
      
      Number of Fisher Scoring iterations: 2
      

# Bell-McAffrey degf work with confint.svyglm()

    Code
      confint(svyglm(api00 ~ api99 + stype, design = dclus1, std.errors = "Bell-McAffrey",
      degf = TRUE))
    Output
                        2.5 %      97.5 %
      (Intercept)  56.2792308 143.4336784
      api99         0.8387537   0.9678245
      stypeH      -33.1110597  -5.6634606
      stypeM      -31.8944852  -4.4219434
      attr(,"degf")
      (Intercept)       api99      stypeH      stypeM 
         6.308131    7.061340    5.304156    8.979698 

---

    Code
      confint(svyglm(as.factor(sch.wide) ~ api99 + stype, design = dclus1, family = "quasibinomial",
      std.errors = "Bell-McAffrey", degf = TRUE))
    Output
                        2.5 %       97.5 %
      (Intercept) -0.50252632  6.077377146
      api99       -0.00549704  0.004222495
      stypeH      -3.04477253  0.830203522
      stypeM      -3.27515434 -0.014736530
      attr(,"degf")
      (Intercept)       api99      stypeH      stypeM 
         6.308131    7.061340    5.304156    8.979698 

