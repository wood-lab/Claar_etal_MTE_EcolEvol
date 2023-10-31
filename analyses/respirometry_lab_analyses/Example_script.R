library(FishResp)

setwd("/Users/treysasser/Desktop")

#convert files. The pyroscience file contains o2 data, but not when the flush pump is on. The aquaresp file contains the flush pump info, but not o2 data. We need to combine both for the package to work

pyroscience.aquaresp(pyroscience.file = "110120Bacteria.txt",
                     aquaresp.file = "110120BacteriaSummary.txt",
                     fishresp.file = "110120BacteriaFishresp.txt",
                     date.format = "MDY",
                     n.chamber = 4,
                     wait.phase = 180,
                     measure.phase = 1200)


pyroscience.aquaresp(pyroscience.file = "110320Bacteria.txt",
                     aquaresp.file = "110320BacteriaSummary.txt",
                     fishresp.file = "110320BacteriaFishresp.txt",
                     date.format = "MDY",
                     n.chamber = 4,
                     wait.phase = 180,
                     measure.phase = 1200)


#convert SMR files
pyroscience.aquaresp(pyroscience.file = "110220Fish.txt",
                     aquaresp.file = "110220FishSummary.txt",
                     fishresp.file = "110220FishFishresp.txt",
                     date.format = "MDY",
                     n.chamber = 4,
                     wait.phase = 180,
                     measure.phase = 1200)

# input necessary information about animals, respirometry chambers and DO units
info <-input.info(ID=c("GS1_3","GS3_3","SG1_2","SG1_3"),
                   Mass=c(2.5836,2.1534,1.9309,1.7846),
                   Volume=c(330,330,330,330),
                   DO.unit="mg/L")

# import background respiration rate (pre, post, or both)
pre <- import.test("110220BacteriaFishresp.txt",
                     info,
                     n.chamber=4,
                     logger= c("FishResp"),
                     plot.temperature = F,
                     plot.oxygen = TRUE)

post <- import.test("110320BacteriaFishresp.txt",
                     info,
                     n.chamber=4,
                     logger= c("FishResp"),
                     plot.temperature = F,
                     plot.oxygen = TRUE)

## import "raw" oxygen data
SMR.raw <- import.meas ("110220FishFishresp.txt",
                       info.data = info,
                       logger= "FishResp",
                       n.chamber=4,
                       date.format = "MDY",
                       start.measure = "11:10:23",
                       stop.measure = "09:43:23",
                       plot.temperature = F,
                       plot.oxygen = T)
                       
##correct for background resp. We decided to use a linear method. From the package description: subtracts a vector of progressively changing microbial consumptions from oxygen consumptions of meas.data. The values of oxygen consumption are linearly predicted from two reference points: oxygen consumption of pre.data and oxygen consumption of post.data.

SMR.clean <- correct.meas (info.data = info,
                          pre.data = pre,
                          post.data = post,
                          meas.data = SMR.raw,
                          method = "linear")

#Mass specific met rate before and after correction for background. QC step to visualize changes in met rate over the course of the experiment. Good for identifying if fish settled
QC.activity(SMR.clean, compare = T)

#Extract target slopes to calculate mass specific metrate
#In this case we will extract the "shallowest" 3 slopes
#e.g., the 3 cycles where the met rate was lowest for each chamber
SMR.slope <- extract.slope(SMR.clean,
						   method="min",
						   n.slope=3,)

#linear regression of the extracted slopes. 
QC.slope(SMR.slope, SMR.clean, chamber = "CH1", current = 1200, alter = 900)
QC.slope(SMR.slope, SMR.clean, chamber = "CH2", current = 1200, alter = 900)
QC.slope(SMR.slope, SMR.clean, chamber = "CH3", current = 1200, alter = 900)
QC.slope(SMR.slope, SMR.clean, chamber = "CH4", current = 1200, alter = 900)

#Calculation of background, absolute, and mass spec metrates
SMR <- calculate.MR(SMR.slope,
					density=1000,
					plot.BR=F,
					plot.MR.abs=F,
					plot.MR.mass=T)

#Export the final dataset as a .txt file. 

results <- export.MR(SMR,
					 file= "G3UW_wet.txt",
					 simplify=T,
					 plot.MS.abs=F,
					 plot.MS.mass=T,
					 plot.MS.fact=F)








