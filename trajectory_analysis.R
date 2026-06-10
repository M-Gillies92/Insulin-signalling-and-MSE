
# Load packages and dataset #

rm(list=ls())
library(nlme)
library(reshape2)
library(ggplot2)
library(haven)

load("data.rdata")
library(dplyr)
library(mice)

attach(data_JUL24)
data_JUL24$ID <- 1:nrow(data_JUL24)
data_JUL24 <- data_JUL24[, c("ID", setdiff(names(data_JUL24), "ID"))]
data_JUL24 <- as.data.frame(data_JUL24)


data_JUL24$HOMA_IR_BBS <- (Insulin_BBS * Glucose_BBS)/22.5
data_JUL24$HOMA_IR_TF3 <- (insulin_TF3 * glucose_TF3)/22.5
data_JUL24$HOMA_IR_TF4 <- (insulin_TF4 * glucose_TF4)/22.5
data_JUL24$HOMA_IR_F24 <- (Insulin_F24 * Glucose_F24)/22.5

# Create refractive spherical equivalent variables for 12 and 15 y/o #

data_JUL24$SER_TF1 <- ff1900 + (ff1901 / 2)
data_JUL24$SEL_TF1 <- ff1910 + (ff1911 / 2)
data_JUL24$SER_TF3 <- fh7310 + (fh7311 / 2)
data_JUL24$SEL_TF3 <- fh7320 + (fh7321 / 2)

# Take the mean of left and right eye to give a mean spherical equivalent #

data_JUL24$AutorefF7 <- (f7vs208 + f7vs213)/2
data_JUL24$AutorefF10 <- (fdms155 + fdms165)/2
data_JUL24$AutorefF11 <- (fevs208 + fevs213)/2
data_JUL24$AutorefTF1 <- (data_JUL24$SER_TF1 + data_JUL24$SEL_TF1)/2
data_JUL24$AutorefTF3 <- (data_JUL24$SER_TF3 + data_JUL24$SEL_TF3)/2

detach(data_JUL24)

attach(data_JUL24)

data_JUL24$age7 <- data_JUL24$f7003c/12
data_JUL24$age10 <- data_JUL24$fd003c/12
data_JUL24$age11 <- data_JUL24$fe003c/12
data_JUL24$age12 <- data_JUL24$ff0011a/12
data_JUL24$age15 <- data_JUL24$fh0011a/12

data_JUL24$Zsc_RminL_sph_7 <- scale(f7vs205 - f7vs210)
data_JUL24$Zsc_RminL_sph_10 <- scale(fdms151 - fdms161)
data_JUL24$Zsc_RminL_sph_11 <- scale(fevs205 - fevs210)
data_JUL24$Zsc_RminL_sph_12 <- scale(ff1900 - ff1910)
data_JUL24$Zsc_RminL_sph_15 <- scale(fh7310 - fh7320)
data_JUL24$Zsc_RminL_cyl_7 <- scale(f7vs206 - f7vs211)
data_JUL24$Zsc_RminL_cyl_10 <- scale(fdms152 - fdms162)
data_JUL24$Zsc_RminL_cyl_11 <- scale(fevs206 - fevs211)
data_JUL24$Zsc_RminL_cyl_12 <- scale(ff1901 - ff1911)
data_JUL24$Zsc_RminL_cyl_15 <- scale(fh7311 - fh7321)

data_JUL24$outlier_age7 <- ifelse(abs(data_JUL24$Zsc_RminL_sph_7) > 4 |
                                    abs(data_JUL24$Zsc_RminL_cyl_7) > 4 |
                                    f7vs206                         >= 4 |
                                    f7vs211                         >= 4, 1, 0)

data_JUL24$outlier_age10 <- ifelse(abs(data_JUL24$Zsc_RminL_sph_10) > 4 |
                                     abs(data_JUL24$Zsc_RminL_cyl_10) > 4 |
                                     fdms152                         >= 4 |
                                     fdms162                         >= 4, 1, 0)

data_JUL24$outlier_age11 <- ifelse(abs(data_JUL24$Zsc_RminL_sph_11) > 4 |
                                     abs(data_JUL24$Zsc_RminL_cyl_11) > 4 |
                                     fevs206                         >= 4 |
                                     fevs211                         >= 4, 1, 0)

data_JUL24$outlier_age12 <- ifelse(abs(data_JUL24$Zsc_RminL_sph_12) > 4 |
                                     abs(data_JUL24$Zsc_RminL_cyl_12) > 4 |
                                     ff1901                         >= 4 |
                                     ff1911                         >= 4, 1, 0)

data_JUL24$outlier_age15 <- ifelse(abs(data_JUL24$Zsc_RminL_sph_15) > 4 |
                                     abs(data_JUL24$Zsc_RminL_cyl_15) > 4 |
                                     fh7311                         >= 4 |
                                     fh7321                         >= 4, 1, 0)
data_JUL24$mse_R_7 <- f7vs205 + (0.5 * f7vs206)
data_JUL24$mse_L_7 <- f7vs210 + (0.5 * f7vs211)
data_JUL24$mse_R_10 <- fdms151 + (0.5 * fdms152)
data_JUL24$mse_L_10 <- fdms161 + (0.5 * fdms162)
data_JUL24$mse_R_11 <- fevs205 + (0.5 * fevs206)
data_JUL24$mse_L_11 <- fevs210 + (0.5 * fevs211)
data_JUL24$mse_R_12 <- ff1900 + (0.5 * ff1901)
data_JUL24$mse_L_12 <- ff1910 + (0.5 * ff1911)
data_JUL24$mse_R_15 <- fh7310 + (0.5 * fh7311)
data_JUL24$mse_L_15 <- fh7320 + (0.5 * fh7321)

data_JUL24$avMSE_7 <- ifelse(data_JUL24$outlier_age7==0, 
                             0.5 * (data_JUL24$mse_R_7 +
                                      data_JUL24$mse_L_7),
                             NA)

data_JUL24$avMSE_10 <- ifelse(data_JUL24$outlier_age10==0, 
                              0.5 * (data_JUL24$mse_R_10 +
                                       data_JUL24$mse_L_10),
                              NA)

data_JUL24$avMSE_11 <- ifelse(data_JUL24$outlier_age11==0, 
                              0.5 * (data_JUL24$mse_R_11 +
                                       data_JUL24$mse_L_11),
                              NA)

data_JUL24$avMSE_12 <- ifelse(data_JUL24$outlier_age12==0, 
                              0.5 * (data_JUL24$mse_R_12 +
                                       data_JUL24$mse_L_12),
                              NA)

data_JUL24$avMSE_15 <- ifelse(data_JUL24$outlier_age15==0, 
                              0.5 * (data_JUL24$mse_R_15 +
                                       data_JUL24$mse_L_15),
                              NA)

# Create combined insulin, refractive error, glucose and HOMA-IR variables #

attach(data_JUL24)


detach(data_JUL24)

data_JUL24$visits <- 5 - is.na(data_JUL24$avMSE_7) - is.na(data_JUL24$avMSE_10) - is.na(data_JUL24$avMSE_11) - is.na(data_JUL24$avMSE_12) - is.na(data_JUL24$avMSE_15)


data_JUL24$avMSE_7[data_JUL24$avMSE_7 < -500] <- NA # https://www.ocl-online.de/en/extreme-axial-length-high-myopia?utm_source=chatgpt.com #

data_JUL24 <- data_JUL24[data_JUL24$c804 == 1, ]
data_JUL24 <- data_JUL24[
  (data_JUL24$age7   > 7  & data_JUL24$age7   <= 8 | is.na(data_JUL24$age7)) &
    (data_JUL24$age10 > 10 & data_JUL24$age10 <= 11 | is.na(data_JUL24$age10)) &
    (data_JUL24$age11 > 11 & data_JUL24$age11 <= 12 | is.na(data_JUL24$age11)) &
    (data_JUL24$age12 > 12 & data_JUL24$age12 <= 13 | is.na(data_JUL24$age12)) &
    (data_JUL24$age15 > 15 & data_JUL24$age15 <= 16 | is.na(data_JUL24$age15)),
]

data_JUL24a <- data_JUL24[,c("ID", "ccc900", "Glucose_BBS", "Insulin_BBS", "HOMA_IR_BBS", "Adiponectin_f9", "f9ms026a", "fddd303", "HBA1C_f9", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

data_JUL24b <- data_JUL24a
data_JUL24b$visits <- as.numeric(data_JUL24b$visits)
data_JUL24b$avMSE_7 <- as.numeric(data_JUL24b$avMSE_7)
data_JUL24b$avMSE_10 <- as.numeric(data_JUL24b$avMSE_10)
data_JUL24b$avMSE_11 <- as.numeric(data_JUL24b$avMSE_11)
data_JUL24b$avMSE_12 <- as.numeric(data_JUL24b$avMSE_12)
data_JUL24b$avMSE_15 <- as.numeric(data_JUL24b$avMSE_15)
data_JUL24b$Glucose_BBS <- scale(data_JUL24b$Glucose_BBS)
data_JUL24b$Glucose_BBS <- as.numeric(data_JUL24b$Glucose_BBS)
data_JUL24b$Insulin_BBS <- scale(data_JUL24b$Insulin_BBS)
data_JUL24b$Insulin_BBS <- as.numeric(data_JUL24b$Insulin_BBS)
data_JUL24b$HOMA_IR_BBS <- scale(data_JUL24b$HOMA_IR_BBS)
data_JUL24b$HOMA_IR_BBS <- as.numeric(data_JUL24b$HOMA_IR_BBS)
data_JUL24b$Adiponectin_f9 <- scale(data_JUL24b$Adiponectin_f9)
data_JUL24b$Adiponectin_f9 <- as.numeric(data_JUL24b$Adiponectin_f9)
data_JUL24b$f9ms026a <- scale(data_JUL24b$f9ms026a)
data_JUL24b$f9ms026a <- as.numeric(data_JUL24b$f9ms026a)
data_JUL24b$fddd303 <- scale(data_JUL24b$fddd303)
data_JUL24b$fddd303 <- as.numeric(data_JUL24b$fddd303)
data_JUL24b$HBA1C_f9 <- scale(data_JUL24b$HBA1C_f9)
data_JUL24b$HBA1C_f9 <- as.numeric(data_JUL24b$HBA1C_f9)
data_JUL24b <- data_JUL24b[which(data_JUL24b$visits >= 3),]
data_JUL24b$ccc900 <- as.numeric(data_JUL24b$ccc900)

impdata <- mice(data_JUL24b, m = 90, maxit = 5, method = "pmm", ridge = 0.01, seed = 500)
impdata1 <- complete(impdata, 1)

library(mice)

data2 <- impdata1[,c("ID", "ccc900", "Glucose_BBS", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensions <- dim(data2)[1]

data3 <- as.data.frame.array(melt(data2, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4 <- as.data.frame.array(melt(data2, id.vars = c("ID", "ccc900", "Glucose_BBS", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3)[1] <- "ID"
names(data3)[2] <- "avMSE"
names(data3)[3] <- "avMSE_D"
names(data4)[1] <- "ID"
names(data4)[2] <- "Sex"
names(data4)[3] <- "Fasting_Glucose"
names(data4)[5] <- "Age"
names(data4)[6] <- "Age_Years"

data5 <- data3[order(data3$ID),]
data6 <- data4[order(data4$ID),]

data5$visit     <- rep(1:5,dimensions)
data6$visit     <- rep(1:5,dimensions)

data7 <- merge(data5, data6, by=c("ID", "visit"))


head(data7)


lme1 <- lme(avMSE_D ~ Sex + Fasting_Glucose + poly(I(Age_Years - 7.53),4) + Fasting_Glucose:I(Age_Years - 7.53), 
            random=~I(Age_Years - 7.53) | ID,
            na.action = na.omit, 
            method="ML",
            correlation = corCAR1(form = ~ visit | ID),
            data=data7)


summary(lme1)

B <- summary(lme1)
atable <- B$tTable
atable

meannn           <- mean(data7$Fasting_Glucose)
cii              <- sd(data7$Fasting_Glucose)


minscores        <- min(data7$Fasting_Glucose)
fifthh           <- meannn-2*cii
maxscores        <- max(data7$Fasting_Glucose)
medscores        <- median(data7$Fasting_Glucose)
ninetyfifthh     <- meannn+2*cii

pdatas <- expand.grid(Age_Years=seq(7,15,by=1), Fasting_Glucose=c(fifthh, meannn,ninetyfifthh), Sex=c(1.5), avMSE_D="1")
pdatas

pdatas[,4]       <-predict(lme1, pdatas, level=0)

pdatas$Fasting_Glucose     <- as.factor(pdatas$Fasting_Glucose)
levels(pdatas$Fasting_Glucose) <- c("Low fasting glucose","Average fasting glucose","High fasting glucose")


figg1 <- ggplot(pdatas, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "Fasting glucose", colour = "Fasting glucose") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = Fasting_Glucose, colour = Fasting_Glucose), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = Fasting_Glucose)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))

figg1


data2in <- impdata1[,c("ID", "ccc900", "Insulin_BBS", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionsin <- dim(data2in)[1]

data3in <- as.data.frame.array(melt(data2in, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4in <- as.data.frame.array(melt(data2in, id.vars = c("ID", "ccc900", "Insulin_BBS", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3in)[1] <- "ID"
names(data3in)[2] <- "avMSE"
names(data3in)[3] <- "avMSE_D"
names(data4in)[1] <- "ID"
names(data4in)[2] <- "Sex"
names(data4in)[3] <- "Fasting_Insulin"
names(data4in)[5] <- "Age"
names(data4in)[6] <- "Age_Years"

data5in <- data3in[order(data3in$ID),]
data6in <- data4in[order(data4in$ID),]

data5in$visit     <- rep(1:5,dimensionsin)
data6in$visit     <- rep(1:5,dimensionsin)

data7in <- merge(data5in, data6in, by=c("ID", "visit"))

head(data7in)


lme1in <- lme(avMSE_D ~ Sex + Fasting_Insulin + poly(I(Age_Years - 7.53),4) + Fasting_Insulin:I(Age_Years - 7.53), 
              random=~I(Age_Years - 7.53) | ID,
              na.action = na.omit, 
              method="ML",
              correlation = corCAR1(form = ~ visit | ID),
              data=data7in)


summary(lme1in)

Bin <- summary(lme1in)
atablein <- Bin$tTable
atablein

meanin           <- mean(data7in$Fasting_Insulin)
ciin              <- sd(data7in$Fasting_Insulin)


minscorein        <- min(data7in$Fasting_Insulin)
fifthin           <- meanin-2*ciin
maxscorein        <- max(data7in$Fasting_Insulin)
medscorein       <- median(data7in$Fasting_Insulin)
ninetyfifthin     <- meanin+2*ciin

indata <- expand.grid(Age_Years=seq(7,15,by=1), Fasting_Insulin=c(fifthin, meanin,ninetyfifthin), Sex=c(1.5), avMSE_D="1")
indata

indata[,4]       <-predict(lme1in, indata, level=0)

indata$Fasting_Insulin     <- as.factor(indata$Fasting_Insulin)
levels(indata$Fasting_Insulin) <- c("Low fasting insulin","Average fasting insulin","High fasting insulin")

figin <- ggplot(indata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "Fasting insulin", colour = "Fasting insulin") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = Fasting_Insulin, colour = Fasting_Insulin), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = Fasting_Insulin)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))


figin


data2ho <- impdata1[,c("ID", "ccc900", "HOMA_IR_BBS", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionsho <- dim(data2ho)[1]

data3ho <- as.data.frame.array(melt(data2ho, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4ho <- as.data.frame.array(melt(data2ho, id.vars = c("ID", "ccc900", "HOMA_IR_BBS", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3ho)[1] <- "ID"
names(data3ho)[2] <- "avMSE"
names(data3ho)[3] <- "avMSE_D"
names(data4ho)[1] <- "ID"
names(data4ho)[2] <- "Sex"
names(data4ho)[3] <- "HOMA_IR"
names(data4ho)[5] <- "Age"
names(data4ho)[6] <- "Age_Years"

data5ho <- data3ho[order(data3ho$ID),]
data6ho <- data4ho[order(data4ho$ID),]

data5ho$visit     <- rep(1:5,dimensionsho)
data6ho$visit     <- rep(1:5,dimensionsho)

data7ho <- merge(data5ho, data6ho, by=c("ID", "visit"))

head(data7ho)

lme1ho <- lme(avMSE_D ~ Sex + HOMA_IR + poly(I(Age_Years - 7.53),4) + HOMA_IR:I(Age_Years - 7.53), 
              random=~I(Age_Years - 7.53) | ID,
              na.action = na.omit, 
              method="ML",
              correlation = corCAR1(form = ~ visit | ID),
              data=data7ho)


summary(lme1ho)

Bho <- summary(lme1ho)
atableho <- Bho$tTable
atableho

meanho           <- mean(data7ho$HOMA_IR)
ciho              <- sd(data7ho$HOMA_IR)


minscoreho        <- min(data7ho$HOMA_IR)
fifthho           <- meanho-2*ciho
maxscoreho        <- max(data7ho$HOMA_IR)
medscoreho       <- median(data7ho$HOMA_IR)
ninetyfifthho     <- meanho+2*ciho

hodata <- expand.grid(Age_Years=seq(7,15,by=1), HOMA_IR=c(fifthho, meanho,ninetyfifthho), Sex=c(1.5), avMSE_D="1")
hodata

hodata[,4]       <-predict(lme1ho, hodata, level=0)

hodata$HOMA_IR     <- as.factor(hodata$HOMA_IR)
levels(hodata$HOMA_IR) <- c("Low HOMA-IR","Average HOMA-IR","High HOMA-IR")

figho <- ggplot(hodata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "HOMA-IR", colour = "HOMA-IR") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = HOMA_IR, colour = HOMA_IR), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = HOMA_IR)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))
figho


data2ad <- impdata1[,c("ID", "ccc900", "Adiponectin_f9", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionsad <- dim(data2ad)[1]

data3ad <- as.data.frame.array(melt(data2ad, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4ad <- as.data.frame.array(melt(data2ad, id.vars = c("ID", "ccc900", "Adiponectin_f9", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3ad)[1] <- "ID"
names(data3ad)[2] <- "avMSE"
names(data3ad)[3] <- "avMSE_D"
names(data4ad)[1] <- "ID"
names(data4ad)[2] <- "Sex"
names(data4ad)[3] <- "Adiponectin"
names(data4ad)[5] <- "Age"
names(data4ad)[6] <- "Age_Years"

data5ad <- data3ad[order(data3ad$ID),]
data6ad <- data4ad[order(data4ad$ID),]

data5ad$visit     <- rep(1:5,dimensionsad)
data6ad$visit     <- rep(1:5,dimensionsad)

data7ad <- merge(data5ad, data6ad, by=c("ID", "visit"))

head(data7ad)


lme1ad <- lme(avMSE_D ~ Sex + Adiponectin + poly(I(Age_Years - 7.53),4) + Adiponectin:I(Age_Years - 7.53), 
              random=~I(Age_Years - 7.53) | ID,
              na.action = na.omit, 
              method="ML",
              correlation = corCAR1(form = ~ visit | ID),
              data=data7ad)


summary(lme1ad)

Bad <- summary(lme1ad)
atablead <- Bad$tTable
atablead

meanad           <- mean(data7ad$Adiponectin)
ciad             <- sd(data7ad$Adiponectin)


minscoread        <- min(data7ad$Adiponectin)
fifthad           <- meanad-2*ciad
maxscoread       <- max(data7ad$Adiponectin)
medscoread       <- median(data7ad$Adiponectin)
ninetyfifthad     <- meanad+2*ciad

addata <- expand.grid(Age_Years=seq(7,15,by=1), Adiponectin=c(fifthad, meanad,ninetyfifthad), Sex=c(1.5), avMSE_D="1")
addata

addata[,4]       <-predict(lme1ad, addata, level=0)

addata$Adiponectin     <- as.factor(addata$Adiponectin)
levels(addata$Adiponectin) <- c("Low adiponectin","Average adiponectin","High adiponectin")

figad <- ggplot(addata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "Adiponectin", colour = "Adiponectin") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = Adiponectin, colour = Adiponectin), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = Adiponectin)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))
figad

data2bmi <- impdata1[,c("ID", "ccc900", "f9ms026a", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionsbmi <- dim(data2bmi)[1]

data3bmi <- as.data.frame.array(melt(data2bmi, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4bmi <- as.data.frame.array(melt(data2bmi, id.vars = c("ID", "ccc900", "f9ms026a", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3bmi)[1] <- "ID"
names(data3bmi)[2] <- "avMSE"
names(data3bmi)[3] <- "avMSE_D"
names(data4bmi)[1] <- "ID"
names(data4bmi)[2] <- "Sex"
names(data4bmi)[3] <- "BMI"
names(data4bmi)[5] <- "Age"
names(data4bmi)[6] <- "Age_Years"

data5bmi <- data3bmi[order(data3bmi$ID),]
data6bmi <- data4bmi[order(data4bmi$ID),]

data5bmi$visit     <- rep(1:5,dimensionsbmi)
data6bmi$visit     <- rep(1:5,dimensionsbmi)

data7bmi <- merge(data5bmi, data6bmi, by=c("ID", "visit"))

head(data7bmi)


lme1bmi <- lme(avMSE_D ~ Sex + BMI + poly(I(Age_Years - 7.53),4) + BMI:I(Age_Years - 7.53), 
               random=~I(Age_Years - 7.53) | ID,
               na.action = na.omit, 
               method="ML",
               correlation = corCAR1(form = ~ visit | ID),
               data=data7bmi)


summary(lme1bmi)

BBMI <- summary(lme1bmi)
atablebmi <- BBMI$tTable
atablebmi

meanbmi           <- mean(data7bmi$BMI)
cibmi             <- sd(data7bmi$BMI)


minscorebmi       <- min(data7bmi$BMI)
fifthbmi           <- meanbmi-2*cibmi
maxscorebmi       <- max(data7bmi$BMI)
medscorebmi       <- median(data7bmi$BMI)
ninetyfifthbmi     <- meanbmi+2*cibmi

bmidata <- expand.grid(Age_Years=seq(7,15,by=1), BMI=c(fifthbmi, meanbmi,ninetyfifthbmi), Sex=c(1.5), avMSE_D="1")
bmidata

bmidata[,4]       <-predict(lme1bmi, bmidata, level=0)

bmidata$BMI     <- as.factor(bmidata$BMI)
levels(bmidata$BMI) <- c("Low BMI","Average BMI","High BMI")

figbmi <- ggplot(bmidata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "BMI", colour = "BMI") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = BMI, colour = BMI), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = BMI)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))
figbmi

data2ch2o <- impdata1[,c("ID", "ccc900", "fddd303", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionsch2o <- dim(data2ch2o)[1]

data3ch2o <- as.data.frame.array(melt(data2ch2o, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4ch2o <- as.data.frame.array(melt(data2ch2o, id.vars = c("ID", "ccc900", "fddd303", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3ch2o)[1] <- "ID"
names(data3ch2o)[2] <- "avMSE"
names(data3ch2o)[3] <- "avMSE_D"
names(data4ch2o)[1] <- "ID"
names(data4ch2o)[2] <- "Sex"
names(data4ch2o)[3] <- "Carbohydrate_Intake"
names(data4ch2o)[5] <- "Age"
names(data4ch2o)[6] <- "Age_Years"

data5ch2o <- data3ch2o[order(data3ch2o$ID),]
data6ch2o <- data4ch2o[order(data4ch2o$ID),]

data5ch2o$visit     <- rep(1:5,dimensionsch2o)
data6ch2o$visit     <- rep(1:5,dimensionsch2o)

data7ch2o <- merge(data5ch2o, data6ch2o, by=c("ID", "visit"))

head(data7ch2o)


lme1ch2o <- lme(avMSE_D ~ Sex + Carbohydrate_Intake + poly(I(Age_Years - 7.53),4) + Carbohydrate_Intake:I(Age_Years - 7.53), 
                random=~I(Age_Years - 7.53) | ID,
                na.action = na.omit, 
                method="ML",
                correlation = corCAR1(form = ~ visit | ID),
                data=data7ch2o)


summary(lme1ch2o)

Bch2o <- summary(lme1ch2o)
atablech2o <- Bch2o$tTable
atablech2o

meanch2o           <- mean(data7ch2o$Carbohydrate_Intake)
cich2o            <- sd(data7ch2o$Carbohydrate_Intake)


minscorech2o       <- min(data7ch2o$Carbohydrate_Intake)
fifthch2o         <- meanch2o-2*cich2o
maxscorech2o       <- max(data7ch2o$Carbohydrate_Intake)
medscorech2o      <- median(data7ch2o$Carbohydrate_Intake)
ninetyfifthch2o     <- meanch2o+2*cich2o

ch2odata <- expand.grid(Age_Years=seq(7,15,by=1), Carbohydrate_Intake=c(fifthch2o,meanch2o,ninetyfifthch2o), Sex=c(1.5), avMSE_D="1")
ch2odata

ch2odata[,4]       <-predict(lme1ch2o, ch2odata, level=0)

ch2odata$Carbohydrate_Intake    <- as.factor(ch2odata$Carbohydrate_Intake)
levels(ch2odata$Carbohydrate_Intake) <- c("Low carbohydrate intake","Average carbohydrate intake","High carbohydrate intake")

figch2o <- ggplot(ch2odata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "Carbohydrate intake", colour = "Carbohydrate intake") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = Carbohydrate_Intake, colour = Carbohydrate_Intake), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = Carbohydrate_Intake)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))

figch2o

data2hba <- impdata1[,c("ID", "ccc900", "HBA1C_f9", "visits", "age7","age10","age11","age12","age15", "avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")]

dimensionshba <- dim(data2hba)[1]

data3hba <- as.data.frame.array(melt(data2hba, id.vars = c("ID"), measure.vars = c("avMSE_7", "avMSE_10", "avMSE_11", "avMSE_12", "avMSE_15")))

data4hba <- as.data.frame.array(melt(data2hba, id.vars = c("ID", "ccc900", "HBA1C_f9", "visits"), measure.vars = c("age7", "age10", "age11", "age12", "age15")))

names(data3hba)[1] <- "ID"
names(data3hba)[2] <- "avMSE"
names(data3hba)[3] <- "avMSE_D"
names(data4hba)[1] <- "ID"
names(data4hba)[2] <- "Sex"
names(data4hba)[3] <- "HbA1c"
names(data4hba)[5] <- "Age"
names(data4hba)[6] <- "Age_Years"

data5hba <- data3hba[order(data3hba$ID),]
data6hba <- data4hba[order(data4hba$ID),]

data5hba$visit     <- rep(1:5,dimensionshba)
data6hba$visit     <- rep(1:5,dimensionshba)

data7hba <- merge(data5hba, data6hba, by=c("ID", "visit"))

head(data7hba)


lme1hba <- lme(avMSE_D ~ Sex + HbA1c + poly(I(Age_Years - 7.53),4) + HbA1c:I(Age_Years - 7.53), 
               random=~I(Age_Years - 7.53) | ID,
               na.action = na.omit, 
               method="ML",
               correlation = corCAR1(form = ~ visit | ID),
               data=data7hba)


summary(lme1hba)

Bhba <- summary(lme1hba)
atablehba <- Bhba$tTable
atablehba

meanhba           <- mean(data7hba$HbA1c)
cihba            <- sd(data7hba$HbA1c)


minscorehba       <- min(data7hba$HbA1c)
fifthhba         <- meanhba-2*cihba
maxscorehba       <- max(data7hba$HbA1c)
medscorehba      <- median(data7hba$HbA1c)
ninetyfifthhba     <- meanhba+2*cihba

hbadata <- expand.grid(Age_Years=seq(7,15,by=1), HbA1c=c(fifthhba,meanhba,ninetyfifthhba), Sex=c(1.5), avMSE_D="1")
hbadata

hbadata[,4]       <-predict(lme1hba, hbadata, level=0)

hbadata$HbA1c    <- as.factor(hbadata$HbA1c)
levels(hbadata$HbA1c) <- c("Low HbA1c","Average HbA1c","High HbA1c")

fighba <- ggplot(hbadata, aes(Age_Years, avMSE_D)) +
  labs(x = "Age (Years)", y = "Average MSE (D)", linetype = "HbA1c", colour = "HbA1c") +
  scale_linetype_manual(values = c(1, 1, 1)) +
  geom_line(linewidth = 2, aes(linetype = HbA1c, colour = HbA1c), show.legend = FALSE) +
  geom_point(size = 3, aes(colour = HbA1c)) +
  scale_shape_manual(values = c(16, 16, 16)) +
  scale_x_continuous(limits = c(7, 15.5), breaks = seq(7, 15, 1)) +
  scale_y_continuous(limits = c(-0.75, 0.50), breaks = seq(-0.75, 0.5, 0.25)) +
  scale_colour_manual(values = c("blue", "forestgreen", "orange")) +
  scale_fill_manual(values = c("#FFFFFF", "#00D800", "#000099", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#D55E00", "#CC79A7")) +
  theme(
    legend.text = element_text(size = 20, color = "black"),
    legend.title = element_text(size = 20, color = "black"),
    legend.position = c(0.75, 0.85),
    legend.key.size = unit(0.5, "lines"), # If you want to position the legend inside the plot area
    axis.text = element_text(size = 20, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key.width = unit(1, "cm") # Move legend key width adjustment here
  ) +
  guides(colour = guide_legend(reverse = TRUE))

fighba