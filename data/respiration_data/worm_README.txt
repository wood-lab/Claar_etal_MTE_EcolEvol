Worm ID: I named the worms using a nomenclature that indicates the fish they came out of. FishID_wormnumber, so GS4_1_1 is the 1st worm out of fish GS4_1.
Host pop: The population the fish came from.
Wet mass: Wet mass of worm in grams.
Tube mass: Mass of drying tube in grams.
Tube+worm mass: Mass of drying tube and dried worm in grams.
Dry mass: Mass of dried worm in grams.
ChVolMinusWormVol: The volume of RPMI in each measurement vial accounting for displacement by the worm. Used to calculate the metabolic rate.
Default settings: The wet mass specific metabolic rate (mg of O2/kg of worm/hour) generated using the default settings of RespR. I included this to see how much of an effect the linear subsetting method had. I wouldn’t use this. Basically, the default setting of the package just does a linear regression over the whole trial and that doesn’t work well when the worms use all of the oxygen in some of the trials (imagine a regression drawn over a slope that looks something like \_).
16C linear method: The wet mass specific metabolic rate (mg of O2/kg of worm/hour) at 16C generated using a function in RespR that performs a rolling regression to find the period in which the change in oxygen content was most linear. I feel confident that this is a good method to assess the actual metabolic rate of the worms. The default for RespR is to output oxygen consumption with a negative sign in front unlike FishResp which doesn’t. So, a met rate of -75 on the worm page and 75 on the fish page are equivalent. I should’ve changed the worm measurements so that they looked the same as the fish measurements, but it didn’t happen. Sorry!
20C linear method: Same as 16C linear, but the measurements were taken when the temp was 20C.
 
I think the most important columns for your analyses will be WormID, WetMass (grams), and 16C Linear method (metabolic rate in mg of O2/kg worm/hour). Hopefully that all makes at least some sense :)
