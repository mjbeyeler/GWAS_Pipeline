phenotype.name <- PHENOTYPE.NAME

colnames(Dgrp2_Inversions)[1] <- 'line.id'
Dgrp2_Inversions$line.id <- as.factor(paste('line_',
                                            substr(Dgrp2_Inversions$line.id, 6, 100),
                                            sep=''))

colnames(Dgrp2_Infection) <- c('line.id', 'infection.status')
Dgrp2_Infection$line.id <- as.factor(paste('line_',
                                           substr(Dgrp2_Infection$line.id, 7, 100),
                                           sep=''))
# colnames(Dgrp_Flystock_Ids)[2] <- 'line.id'
# Dgrp_Flystock_Ids$line.id <- as.factor(paste('line_',
#                                              substr(Dgrp_Flystock_Ids$line.id, 6, 100),
#                                              sep=''))


# Giving standardized column names, dropping NA values, and DROPPING NON-DGRP2 lines
if(SEXUAL.DIMORPHISM == T) {
  colnames(Phenotype_Raw) <- c('line.id', 'phenotype', 'sex')
  Phenotype_Raw <- Phenotype_Raw[!is.na(Phenotype_Raw$phenotype), ]
  levels.to.drop <- levels(Phenotype_Raw$line.id)[!levels(Phenotype_Raw$line.id)
                                                  %in% levels(Dgrp2_Infection$line.id)]
  Phenotype_Raw <- Phenotype_Raw[!Phenotype_Raw$line.id %in% levels.to.drop, ]
  Female <- Phenotype_Raw[Phenotype_Raw$sex=='f', -3]
  Male <- Phenotype_Raw[Phenotype_Raw$sex=='m', -3]
  
  
  # Lastly, discarding all lines that do not occurr in BOTH sexes.
  print(paste('Unique to Female: ',
              unique(Female$line.id[which(!Female$line.id %in% Male$line.id)])))
  Female <- Female[which(Female$line.id %in% Male$line.id), ]
  print(paste('Unique to Male: ',
              unique(Male$line.id[which(!Male$line.id %in% Female$line.id)])))
  Male <- Male[which(Male$line.id %in% Female$line.id), ]
  
} else {
  colnames(Phenotype_Raw) <- c('line.id', 'phenotype')
  Phenotype_Raw <- Phenotype_Raw[!is.na(Phenotype_Raw$phenotype), 1:2]
  levels.to.drop <- levels(Phenotype_Raw$line.id)[!levels(Phenotype_Raw$line.id)
                                                  %in% levels(Dgrp2_Infection$line.id)]
  Phenotype_Raw <- Phenotype_Raw[!Phenotype_Raw$line.id %in% levels.to.drop, ]
}

# This part is specific to mass data, because they were stored in Flystock ID format.
# Female$line.id <- Dgrp_Flystock_Ids$line.id[match(Phenotype_Raw$line.id, Dgrp_Flystock_Ids$Stock.No)][Phenotype_Raw$sex=='f']
# Male$line.id <- Dgrp_Flystock_Ids$line.id[match(Phenotype_Raw$line.id, Dgrp_Flystock_Ids$Stock.No)][Phenotype_Raw$sex=='m']
# # Here, I assume there was a typo. Cf. protocol
# Female$line.id[is.na(Female$line.id)] = 'line_355'
# Male$line.id[is.na(Male$line.id)] = 'line_355'
# #End of mass-specific part.


if(SEXUAL.DIMORPHISM == T) {
  occurring.levels <- levels(droplevels(Female$line.id))
} else {
  occurring.levels <- levels(droplevels(Phenotype_Raw$line.id))
}



phenotype.adjusted = F

# CONTINUATION SEXUALLY DIMORPHIC -----------------------------------------
if(SEXUAL.DIMORPHISM == T) {
  
  
  # Assessing Normality -----------------------------------------------------
  
  alpha.normality <- NORMALITY.SIGNIFICANCE.LEVEL
  
  # Using visuals
  NormalityHistogram(c(Female$phenotype, Male$phenotype))
  NormalityHistogram(Female$phenotype)
  NormalityHistogram(Male$phenotype)
  
  qqnorm(Female$phenotype)
  qqline(Female$phenotype,col=2)
  
  qqnorm(Male$phenotype)
  qqline(Male$phenotype,col=2)
  
  # Using a goodness of fit chi-squared test
  # p.total <- ChisqForNormality(c(Female$phenotype, Male$phenotype))
  # print(p.total)
  # p.female <- ChisqForNormality(Female$phenotype)
  # p.male <- ChisqForNormality(Male$phenotype)
  # print(p.female)
  # print(p.male)
  
  
  # Shapiro test has a sample size limit
  if(length(c(Female$phenotype, Male$phenotype)) < 5000) {
    p.total <- shapiro.test(c(Female$phenotype, Male$phenotype))$p.value
  } else {
    p.total <- ChisqForNormality(c(Female$phenotype, Male$phenotype))
  }
  print(paste('p-value total', p.total))
  
  
  p.female <- shapiro.test(Female$phenotype)$p.value
  p.male <- shapiro.test(Male$phenotype)$p.value
  print(p.female)
  print(p.male)
  mean.female.raw <- sapply(occurring.levels, function(x) mean(Female$phenotype[Female$line.id==x]))
  mean.male.raw <- sapply(occurring.levels, function(x) mean(Male$phenotype[Male$line.id==x]))
  # Log-transforming the phenotype if male or female is not distributed normally.
  if(p.female < alpha.normality | p.male < alpha.normality) {
    print('At least one of the phenotypes was not normally distributed. Log-transformation performed.')
    Female$phenotype <- log(Female$phenotype)
    Male$phenotype <- log(Male$phenotype)
    phenotype.adjusted = T
  }
  
  
  # Writing a table that can be used for Online GWAS
  mean.female <- sapply(occurring.levels, function(x) mean(Female$phenotype[Female$line.id==x]))
  mean.male <- sapply(occurring.levels, function(x) mean(Male$phenotype[Male$line.id==x]))
  write.table(data.frame(substr(occurring.levels, 6, 100), mean.female), 
              file=paste('Outputs/Gwas-Online-Input-', phenotype.name,
                         'Female.csv', sep=''),
              sep=',', row.names=F, col.names=F, quote=F)
  
  write.table(data.frame(substr(occurring.levels, 6, 100), mean.male), 
              file=paste('Outputs/Gwas-Online-Input-', phenotype.name,
                         'Male.csv', sep=''),
              sep=',', row.names=F, col.names=F, quote=F)
  
  # Phenotype difference
  phenotype.diff <- (mean.female.raw - mean.male.raw) / ((mean.male.raw + mean.female.raw) / 2)
  
  write.table(data.frame(substr(occurring.levels, 6, 100), phenotype.diff), 
              file=paste('Outputs/Gwas-Online-Input-', phenotype.name,
                         'Dimorphism.csv', sep=''),
              sep=',', row.names=F, col.names=F, quote=F)
  
  
  # Phenotype Adjustment ----------------------------------------------------
  
  used.inversions <- INVERSIONS.CONSIDERED
  
  Adjustment_Data_Female <- data.frame(phenotype=Female$phenotype,
                                       Dgrp2_Infection[match(Female$line.id,
                                                             Dgrp2_Infection$line.id), ],
                                       Dgrp2_Inversions[match(Female$line.id,
                                                              Dgrp2_Inversions$line.id),
                                                        used.inversions]
  )
  Adjustment_Data_Male <- data.frame(phenotype=Male$phenotype,
                                     Dgrp2_Infection[match(Male$line.id,
                                                           Dgrp2_Infection$line.id), ],
                                     Dgrp2_Inversions[match(Male$line.id,
                                                            Dgrp2_Inversions$line.id),
                                                      used.inversions]
  )
  
  # Writing the whole model to a table.
  write.table(Adjustment_Data_Female, file=paste('Outputs/',
                                                 phenotype.name,
                                                 '-Adjustment-Data-Female.txt', sep=''),
              sep='\t', row.names=F)
  write.table(Adjustment_Data_Female, file=paste('Outputs/',
                                                 phenotype.name,
                                                 '-Adjustment-Data-Male.txt', sep=''),
              sep='\t', row.names=F)
  
  female.model <- with(Adjustment_Data_Female,
                       lmer(phenotype ~ infection.status + In.2L.t + In.2R.NS + 
                              In.3R.P + In.3R.K + In.3R.Mo + (1|line.id)))
  male.model <- with(Adjustment_Data_Male,
                     lmer(phenotype ~ infection.status + In.2L.t + In.2R.NS + 
                            In.3R.P + In.3R.K + In.3R.Mo + (1|line.id)))
  
  Adjustment_Prediction_Data <- expand.grid(line.id=occurring.levels,
                                            infection.status=factor('n', levels=c('n', 'y')),
                                            In.2L.t=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.2R.NS=factor("ST",
                                                            levels=c("ST", "INV/ST", "INV")),
                                            In.3R.P=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.3R.K=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.3R.Mo=factor("ST",
                                                            levels=c("ST", "INV/ST", "INV")))
  
  adjusted.female <- predict(female.model, Adjustment_Prediction_Data)
  adjusted.male <- predict(male.model, Adjustment_Prediction_Data)
  plot(mean.female, adjusted.female)
  
  phenotype.diff.adj <- (adjusted.male - adjusted.female) / ((adjusted.male + adjusted.female) / 2)
  plot(phenotype.diff, phenotype.diff.adj)
  
  
  # Writing the Adjusted Phenotypes to Be Used for FaST-LMM -----------------
  
  
  if(phenotype.adjusted == T) {
    write.table('# Phenotype was log-transformed!',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Female-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
    write.table('# Phenotype was log-transformed!',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Male-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
    write.table('# Phenotype was log-transformed!',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Dimorphism-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
  }
  else {
    write.table('# Phenotype was not log-transformed.',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Female-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
    write.table('# Phenotype was not log-transformed.',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Male-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
    write.table('# Phenotype was not log-transformed.',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Dimorphism-Commented.txt', sep=''),
                col.names=F, row.names=F, quote=F)
  }
  write.table(data.frame(occurring.levels, occurring.levels, adjusted.female),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Female-Commented.txt', sep=''),
              col.names=F, row.names=F, quote=F, append=T)
  
  write.table(data.frame(occurring.levels, occurring.levels, adjusted.male),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Male-Commented.txt', sep=''),
              col.names=F, row.names=F, quote=F, append=T)
  
  # Writing the difference one, for now using raw difference
  write.table(data.frame(occurring.levels, occurring.levels, phenotype.diff),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Dimorphism-Commented.txt', sep=''),
              col.names=F, row.names=F, quote=F, append=T)
  
  
  # For compatibility-reasons with FaST-LMM, I also have to write an uncommented version!
  
  write.table(data.frame(occurring.levels, occurring.levels, adjusted.female),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Female.txt', sep=''),
              col.names=F, row.names=F, quote=F)
  
  write.table(data.frame(occurring.levels, occurring.levels, adjusted.male),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Male.txt', sep=''),
              col.names=F, row.names=F, quote=F)
  
  # Writing the difference one, for now using raw difference
  write.table(data.frame(occurring.levels, occurring.levels, phenotype.diff),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Dimorphism.txt', sep=''),
              col.names=F, row.names=F, quote=F)
  
  
  
  
  # CONTINUATION NON-DIMORPHIC PHENOTYPE --------------------------------------
} else {
  
  # Assessing Normality -----------------------------------------------------
  
  alpha.normality <- NORMALITY.SIGNIFICANCE.LEVEL
  
  # Using visuals
  NormalityHistogram(c(Phenotype_Raw$phenotype))
  
  qqnorm(Phenotype_Raw$phenotype)
  qqline(Phenotype_Raw$phenotype,col=2)
  
  # Using a goodness of fit chi-squared test
  if(length(c(Female$phenotype, Male$phenotype)) < 5000) {
    p.total <- shapiro.test(c(Female$phenotype, Male$phenotype))$p.value
  } else {
    p.total <- ChisqForNormality(c(Female$phenotype, Male$phenotype))
  }
  print(p.total)
  
  # Log-transforming the phenotype if male or female is not distributed normally.
  if(p.total < alpha.normality) {
    phenotype.adjusted = T
    print('The phenotype was not normally distributed. Log-transformation performed.')
    Phenotype_Raw$phenotype <- log(Phenotype_Raw$phenotype)
  }
  
  
  # Writing a table that can be used for Online GWAS
  mean.total <- sapply(occurring.levels, function(x) mean(Phenotype_Raw$phenotype[Phenotype_Raw$line.id==x]))
  write.table(data.frame(substr(occurring.levels, 6, 100), mean.total), 
              file=paste('Outputs/Gwas-Online-Input-', phenotype.name,
                         '.csv', sep=''),
              sep=',', row.names=F, col.names=F, quote=F)
  
  
  
  # Phenotype Adjustment ----------------------------------------------------
  
  used.inversions <- INVERSIONS.CONSIDERED
  Adjustment_Data_Total <- data.frame(phenotype=Phenotype_Raw$phenotype,
                                      Dgrp2_Infection[match(Phenotype_Raw$line.id,
                                                            Dgrp2_Infection$line.id), ],
                                      Dgrp2_Inversions[match(Phenotype_Raw$line.id,
                                                             Dgrp2_Inversions$line.id),
                                                       used.inversions]
  )
  
  # Writing the whole model to a table.
  write.table(Adjustment_Data_Total, file=paste('Outputs/',
                                                phenotype.name,
                                                '-Adjustment-Data-Total.txt', sep=''),
              sep='\t', row.names=F)
  total.model <- with(Adjustment_Data_Total,
                      lmer(phenotype ~ infection.status + In.2L.t + In.2R.NS + 
                             In.3R.P + In.3R.K + In.3R.Mo + (1|line.id)))
  Adjustment_Prediction_Data <- expand.grid(line.id=occurring.levels,
                                            infection.status=factor('n', levels=c('n', 'y')),
                                            In.2L.t=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.2R.NS=factor("ST",
                                                            levels=c("ST", "INV/ST", "INV")),
                                            In.3R.P=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.3R.K=factor("ST",
                                                           levels=c("ST", "INV/ST", "INV")),
                                            In.3R.Mo=factor("ST",
                                                            levels=c("ST", "INV/ST", "INV")))
  adjusted.total <- predict(total.model, Adjustment_Prediction_Data)
  plot(mean.total, adjusted.total)
  
  # Writing the Adjusted Phenotypes to Be Used for FaST-LMM -----------------
  
  if(phenotype.adjusted == T) {
    write.table('# Phenotype was log-transformed!',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Total.txt', sep=''),
                col.names=F, row.names=F, quote=F)
  }
  else {
    write.table('# Phenotype was not log-transformed.',
                file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                           '-Total.txt', sep=''),
                col.names=F, row.names=F, quote=F)
  }
  
  write.table(data.frame(occurring.levels, occurring.levels, adjusted.total),
              file=paste('Outputs/Fast-Lmm-Input-', phenotype.name,
                         '-Total.txt', sep=''),
              col.names=F, row.names=F, quote=F, append=T)
  
}


# Lastly, this script also writes all the ocurring lines in the raw phenotype to a text file that can be read by Plink. For convenience
# WriteBare(data.frame(droplevels(Phenotype_Raw$line.id), droplevels(Phenotype_Raw$line.id)), paste('../Outputs/Plink-Phenotypes-', PHENOTYPE.NAME, '.txt', sep=''))

# mass-specific:
WriteBare(data.frame(levels(droplevels(Female$line.id)), levels(droplevels(Female$line.id))), paste('Outputs/Plink-Lines-', PHENOTYPE.NAME, '.txt', sep=''))

