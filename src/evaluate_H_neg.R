## This program is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or (at your option) any later
## version.

## This program is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## this program.  If not, see https://www.gnu.org/licenses/.

## Loading necessary packages
library("tidyverse")
library("xcms")
library("ggpubr")

ppm <- function(x, # Comparable quantity
                y) # Measure
{

    ppm <- 10^6 * (x - y) / y

    return(ppm)

}

## Annotate compounds in XCMS peaklists (e.g. internal standards)
annotate.compound  <-  function(data, # XCMS peaklist
                                compound.name, # name(s) of compound(s)
                                compound.mz, # exact mass(es) of compound(s)
                                compound.rt, # retention time(s) of compound(s)
                                ppmlim = 5, # tolerated mass deviation in parts per million
                                rtlim = 10) # tolerated retention time in seconds
{

    data$Compound <- NA
    
    for(i in 1:length(data[,1])) {

        delta.ppm <- NULL
        delta.rt <- NULL

        for(j in 1:length(compound.name)){

            delta.ppm[j] <- abs(ppm(data$mz[i], compound.mz[j]))
            delta.rt[j] <- abs(data$rt[i] - compound.rt[j])
        }

        index <- which(delta.ppm < ppmlim & delta.rt < rtlim)


        if(length(index) != 0) {

            data$Compound[i] <- as.character(compound.name[index])

        }

    }

    return(data)
    
}

##
## Preprocessing
##

## Load file paths and parameters
files <- list.files(path = getwd(), recursive = TRUE, full.names = TRUE, pattern = ".mzXML")
parameter <- c(8.5, 77, 1.5, 8, 0.004, 6, 10000, 1.5)

## Preprocess files
set <- xcmsSet(files, method = "centWave", peakwidth = c(parameter[1], parameter[2]),
               ppm = parameter[3], snthresh = parameter[4], mzdiff = parameter[5],
               prefilter = c(parameter[6], parameter[7]))
set <- group.density(set, bw = parameter[8])
set <- retcor(set, method = "obiwarp", plottype = "none")
set <- group.density(set, bw = parameter[8])
set <- retcor(set, method = "obiwarp", plottype = "none")
set <- group.density(set, bw = parameter[8])
set.filled <- fillPeaks(set)


##
## Start of evaluation
##

 load("evaluate.RData")

Classes <- set.filled$class
class.levels <- levels(Classes)
peaklist <- peakTable(set.filled)
peaklist.annotated <- annotate.compound(data = peaklist,
                                        compound.name = c("Trp-d5", "Palmitic", "Glucose"),
                                        compound.mz = c(208.1134, 255.2329, 179.0561),
                                        compound.rt = c(458, 86, 146),
                                        rtlim = 2)
write.csv(peaklist.annotated, file = "peaklist_annotated.csv")

rownames(peaklist) <- groupnames(set.filled)
class.levels <- levels(set.filled$class)

## Plot Trp-d5
is <- as.data.frame(matrix(ncol = 2, nrow = 23, dimnames = list(NULL, c("area", "experiment"))))
is$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Trp-d5"),
                                (8 + length(class.levels)):length(peaklist)])
is$experiment <- set.filled$class

p.is <- ggplot(is, aes(x = as.factor(experiment), y = is[,1] )) +
    geom_boxplot(aes(group = experiment)) +
    stat_compare_means(method = "anova", label.y = 3.5*10^6)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 3.3*10^6) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw() 

## Plot Palmitic Acid
palmitic <- as.data.frame(matrix(ncol = 3, nrow = 23, dimnames = list(NULL, c("area", "methanol", "experiment"))))
palmitic$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Palmitic"),(8 + length(class.levels)):length(peaklist)])
colnames(palmitic[,1]) <- "area"
palmitic$methanol <- set.filled$class

p.palmitic <- ggplot(palmitic, aes(x = as.factor(methanol), y = palmitic[,1] )) +
    geom_boxplot(aes(group = methanol)) +
    stat_compare_means(method = "anova", label.y = 7*10^6)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 6.5*10^6) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw()

## Plot Glucose
gluc <- as.data.frame(matrix(ncol = 3, nrow = 23, dimnames = list(NULL, c("area", "methanol", "experiment"))))
gluc$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Glucose"),(8 + length(class.levels)):length(peaklist)])
colnames(gluc[,1]) <- "area"
gluc$methanol <- set.filled$class

p.gluc <- ggplot(gluc, aes(x = as.factor(methanol), y = gluc[,1] )) +
    geom_boxplot(aes(group = methanol)) +
    stat_compare_means(method = "anova", label.y = 8*10^6)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 7.5*10^6) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw()

## Plot Peaknumber
peaknumber <- as.data.frame(t(peakTable(set)[,15:(14+length(set$class))]))
peaknumber[is.na(peaknumber) == FALSE] <- 1
peaknumber[is.na(peaknumber) == TRUE] <- 0
peaknumber$sum <- apply(peaknumber, 1, sum)
peaknumber$experiment <- set$class

p.peaknumber <- ggplot(peaknumber, aes(x = experiment, y = sum)) +
    geom_boxplot(aes(group = experiment)) +
    stat_compare_means(method = "anova", label.y = 1900)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 1800) +
    xlab("Methanol Fraction, %") +
    ylab("Number of Features") +
    theme_bw() 

ggsave("IS_Area.pdf", p.is, width = 8, height = 6, device = "pdf")
ggsave("Palmitic_Acid_Area.pdf", p.palmitic, width = 8, height = 6, device = "pdf")
ggsave("Glucose_Area.pdf", p.gluc, width = 8, height = 6, device = "pdf")
ggsave("Peaknumbers.pdf", p.peaknumber, width = 8, height = 6, device = "pdf")

## Plot internal standard Trp-d5
groupid <- which(peaklist.annotated$Compound == "Trp-d5")
eic <- getEIC(set, groupidx = groupid, rt = "corrected")
pdf("Trp-d5_EICs.pdf", width = 12, height = 8)
plot(eic, set, groupidx = groupnames(eic))
dev.off()

save.image("evaluate.RData")
