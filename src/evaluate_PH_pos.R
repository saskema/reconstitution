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

evalANOVA <- function (features,
                       names,
                       classes,
                       p.value = 0.05,
                       plot = TRUE,
                       labels = FALSE) {

    summary <- as.data.frame(matrix(ncol = 4,
                                    nrow = length(names),
                                    dimnames = list(NULL, c("feature",
                                                            "pvalue",
                                                            "Significant",
                                                            "label"))
                                    )
                             )
    
    summary$feature <- names

    for(i in 1:length(features[1,])) {

        sum <- unlist(summary(aov(features[,i]~classes)))
        summary$pvalue[i] <- sum["Pr(>F)1"]
    }
    summary$Significant <- ifelse(summary$pvalue < p.value, "TRUE", "FALSE")
    summary$label <- names
    summary$label[summary$Significant == "FALSE"] <- ""

    if(plot == TRUE) {

        plot

        coloursBoolean <- c("TRUE" = "red", "FALSE" = "grey")
        
        if(labels == TRUE) {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                geom_text_repel(point.padding = 0.2, colour = "black") +
                ggtitle("ANOVA Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        }else {

            print(
                ggplot(data = summary, aes(x = feature, y = -log10(pvalue),
                                           label = label, colour = Significant)) +
                geom_point(size = 2) +
                scale_color_manual(values = coloursBoolean) +
                geom_hline(yintercept = -log10(p.value), colour = "red") +
                ggtitle("ANOVA Results") +
                xlab("Features") +
                ylab("-log10(p-Values)") +
                theme_classic() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank()
                      )
            )

        }


        }

    return(summary)

}


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
parameter <- c(8.9, 42, 1.5, 35, -0.008, 7, 10000, 0.2)

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
                                        compound.name = c("Trp-d5", "Creatinine", "GPCho160160"),
                                        compound.mz = c(210.1291, 114.0662, 734.5694),
                                        compound.rt = c(146, 35, 616),
                                        rtlim = 2)
write.csv(peaklist.annotated, file = "peaklist_annotated.csv")

rownames(peaklist) <- groupnames(set.filled)
class.levels <- levels(set.filled$class)

## Plot Trp-d5
is <- as.data.frame(matrix(ncol = 3, nrow = 23, dimnames = list(NULL, c("area", "methanol", "experiment"))))
is$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Trp-d5"),(8 + length(class.levels)):length(peaklist)])
colnames(is[,1]) <- "area"
is$methanol <- set.filled$class

p.is <- ggplot(is, aes(x = as.factor(methanol), y = is[,1] )) +
    geom_boxplot(aes(group = methanol)) +
    stat_compare_means(method = "anova", label.y = 1.5*10^8)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 1.4*10^8) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw()

## Plot Creatinine
crea <- as.data.frame(matrix(ncol = 3, nrow = 23, dimnames = list(NULL, c("area", "methanol", "experiment"))))
crea$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Creatinine"),(8 + length(class.levels)):length(peaklist)])
colnames(crea[,1]) <- "area"
crea$methanol <- set.filled$class

p.crea <- ggplot(crea, aes(x = as.factor(methanol), y = crea[,1] )) +
    geom_boxplot(aes(group = methanol)) +
    stat_compare_means(method = "anova", label.y = 8.5*10^8)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 8.2*10^8) +
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
    stat_compare_means(method = "anova", label.y = 3.2*10^7)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 3*10^7) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw()


## Plot Peaknumber
peaknumber <- as.data.frame(t(peakTable(set)[,15:(14+length(set$class))]))
peaknumber[is.na(peaknumber) == FALSE] <- 1
peaknumber[is.na(peaknumber) == TRUE] <- 0
peaknumber$sum <- apply(peaknumber, 1, sum)
peaknumber$methanol <- set$class

p.peaknumber <- ggplot(peaknumber, aes(x = methanol, y = sum)) +
    geom_boxplot(aes(group = methanol)) +
    stat_compare_means(method = "anova", label.y = 1400)+
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "QC", label.y = 1300) +
    xlab("Methanol Fraction, %") +
    ylab("Number of Features") +
    theme_bw() 

ggsave("IS_Area.pdf", p.is, width = 8, height = 6, device = "pdf")
ggsave("Creatinine_Area.pdf", p.crea, width = 8, height = 6, device = "pdf")
ggsave("Glucose_Area.pdf", p.gluc, width = 8, height = 6, device = "pdf")
ggsave("Peaknumbers.pdf", p.peaknumber, width = 8, height = 6, device = "pdf")

## Plot internal standard Trp-d5
groupid <- which(peaklist.annotated$Compound == "Trp-d5")
eic <- getEIC(set, groupidx = groupid, rt = "corrected")
pdf("Trp-d5_EICs.pdf", width = 12, height = 8)
plot(eic, set, groupidx = groupnames(eic))
dev.off()

## Plot Creatinine
groupid <- which(peaklist.annotated$Compound == "Creatinine")
eic <- getEIC(set, groupidx = groupid, rt = "corrected")
pdf("Creatinine.pdf", width = 12, height = 8)
plot(eic, set, groupidx = groupnames(eic))
dev.off()


save.image("evaluate.RData")
