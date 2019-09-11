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
set <- fillPeaks(set)


##
## Start of evaluation
##

## load("evaluate.RData")

Classes <- set$class
class.levels <- levels(Classes)
peaklist <- peakTable(set)
peaklist.annotated <- annotate.compound(peaklist, "Trp-d5", 208.1134, 458, rtlim = 5)
write.csv(peaklist.annotated, file = "peaklist_annotated.csv")

rownames(peaklist) <- groupnames(set)
class.levels <- levels(set$class)
is <- as.data.frame(matrix(ncol = 2, nrow = 23, dimnames = list(NULL, c("area", "experiment"))))
is$area <- t(peaklist.annotated[which(peaklist.annotated$Compound == "Trp-d5"),
                                (8 + length(class.levels)):length(peaklist)])
is$experiment <- set$class

p.is <- ggplot(is, aes(x = as.factor(experiment), y = is[,1] )) +
    geom_boxplot(aes(group = experiment)) +
    xlab("Methanol Fraction, %") +
    ylab("Area") +
    theme_bw() 


peaknumber <- as.data.frame(t(peaklist[,8:14]))
peaknumber$sum <- apply(peaknumber, 1, sum)
peaknumber$methanol <- c("0", "10", "20", "30", "40", "50", "QC")

p.peaknumber <- ggplot(peaknumber, aes(x = methanol, y = sum)) +
    geom_col() +
    xlab("Methanol Fraction, %") +
    ylab("Number of Features") +
    theme_bw() 

ggsave("IS_Area.pdf", p.is, width = 8, height = 6, device = "pdf")
ggsave("Peaknumbers.pdf", p.peaknumber, width = 8, height = 6, device = "pdf")

## Plot internal standard Trp-d5
groupid <- which(peaklist.annotated$Compound == "Trp-d5")
eic <- getEIC(set, groupidx = groupid, rt = "corrected")
pdf("Trp-d5_EICs.pdf", width = 12, height = 8)
plot(eic, set, groupidx = groupnames(eic))
dev.off()

save.image("evaluate.RData")
