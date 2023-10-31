# Install necessary libraries
library(tidyverse)
library(wesanderson)
library(lme4)
library(patchwork)
library(emmeans)
library(infer)

# Load necessary data
resp_15C_SGGS_all <- read.csv("data/respiration_data/MTECollab_final_corrected.csv")

resp_15C_SGGS_all %>% 
  group_by(Exposed,Pop) %>%  
  summarize(exposed = n())

# Remove flagged individuals
resp_15C_SGGS_all <- resp_15C_SGGS_all %>% 
  filter(Notes=="")

resp_15C_SGGS_all %>% group_by(Exposed,Pop) %>%  summarize(exposed = n())

# Standardize Infection column ## NOTE: This is now only Y/N, and unexposed are lumped with exposed but uninfected.
resp_15C_SGGS_all$Infected[resp_15C_SGGS_all$Infected=="N "]<-"N"
resp_15C_SGGS_all$Infected[is.na(resp_15C_SGGS_all$Infected)==TRUE]<-"N"
resp_15C_SGGS_all<-droplevels(resp_15C_SGGS_all)


# Remove infected individuals
resp_15C_SGGS_all_Inf <- resp_15C_SGGS_all %>% 
  filter(Infected=="Y")
resp_15C_SGGS_all_Uninf <- resp_15C_SGGS_all %>% 
  filter(Infected!="Y")

E=0.63
k=8.62*10^(-5)
Temp=273.15+15
conv=0.00000390625 # 450kJ/molO2 converted to kJ/s

resp_15C_SGGS_all_Uninf <- resp_15C_SGGS_all_Uninf %>% 
  mutate(T1_avg=as.numeric(as.character(W_Mavg))) %>% 
  rename(Population=Pop) %>% 
  mutate(Mass=WetMass) %>% 
  mutate(B=T1_avg*conv) %>% 
  mutate(I=B*(WetMass)) %>% # I = B*M 
  mutate(y1=log(I*exp(E/(k*Temp))),
         T1_mass_conv=log(WetMass))

resp_15C_SGGS_all_Inf <- resp_15C_SGGS_all_Inf %>% 
  mutate(T1_avg=as.numeric(as.character(W_Mavg))) %>% 
  mutate_at(vars(WormWetMass), ~replace_na(., 0)) %>% 
  rename(Population=Pop) %>% 
  mutate(Mass=WetMass-WormWetMass) %>% 
  mutate(B=T1_avg*conv) %>% 
  mutate(I=B*(WetMass-WormWetMass)) %>% # I = B*M 
  mutate(y1=log(I*exp(E/(k*Temp))),
         T1_mass_conv=log(WetMass-WormWetMass))

resp_15C_SGGS_all <- resp_15C_SGGS_all %>% 
  mutate(T1_avg=as.numeric(as.character(W_Mavg))) %>%
  mutate_at(vars(WormWetMass), ~replace_na(., 0)) %>% 
  rename(Population=Pop) %>% 
  mutate(Mass=WetMass-WormWetMass) %>% 
  mutate(B=T1_avg*conv) %>% 
  mutate(I=B*(WetMass-WormWetMass)) %>% # I = B*M 
  mutate(y1=log(I*exp(E/(k*Temp))),
         T1_mass_conv=log(WetMass-WormWetMass))

resp_15C_SGGS_all_Uninf %>% # mass-corrected respiration for manuscript
  group_by(Population) %>%#each population separately
  summarise(n=n(),
            mean.resp=mean(T1_avg),sd.resp=sd(T1_avg),
            mean.size=mean(T1_mass_conv),sd.mass=sd(T1_mass_conv),
            mean.y=mean(y1),sd.y=sd(y1))

resp_15C_SGGS_all_Inf %>% # mass-corrected respiration for manuscript
  group_by(Population) %>% #each population separately
  summarise(n=n(),
            mean.resp=mean(T1_avg,na.rm=TRUE),
            sd.resp=sd(T1_avg,na.rm=TRUE),
            max.resp=max(T1_avg,na.rm=TRUE),
            min.resp=min(T1_avg,na.rm=TRUE))

resp_15C_SGGS_all_Uninf %>% # mass-corrected respiration for manuscript
  summarise(n=n(),
            min.resp=min(T1_avg),max.resp=max(T1_avg),
            mean.resp=mean(T1_avg),sd.resp=sd(T1_avg),
            mean.size=mean(T1_mass_conv),sd.mass=sd(T1_mass_conv),
            mean.y=mean(y1),sd.y=sd(y1))

resp_15C_SGGS_all_Inf %>% # mass-corrected respiration for manuscript
  summarise(n=n(),
            mean.resp=mean(T1_avg,na.rm=TRUE),
            sd.resp=sd(T1_avg,na.rm=TRUE),
            max.resp=max(T1_avg,na.rm=TRUE),
            min.resp=min(T1_avg,na.rm=TRUE))

resp_15C_SG_Uninf <-
  resp_15C_SGGS_all_Uninf %>% filter(Population=="SG")
resp_15C_GS_Uninf <-
  resp_15C_SGGS_all_Uninf %>% filter(Population=="GS")

resp_15C_SG_Inf <-
  resp_15C_SGGS_all_Inf %>% filter(Population=="SG")
resp_15C_GS_Inf <-
  resp_15C_SGGS_all_Inf %>% filter(Population=="GS")

resp_15C_SGGS_all_Uninf %>% # whole organism respiration for manuscript
  group_by(Population) %>%
  summarise(n=n(),
            mean.resp.c=mean(T1_avg*Mass,na.rm=TRUE),
            sd.resp.c=sd(T1_avg*Mass,na.rm=TRUE),
            max.resp.c=max(T1_avg*Mass,na.rm=TRUE),
            min.resp.c=min(T1_avg*Mass,na.rm=TRUE))

resp_15C_SGGS_all_Inf %>% # whole organism respiration for manuscript
  group_by(Population) %>%
  summarise(n=n(),
            mean.resp.c=mean(T1_avg*Mass,na.rm=TRUE),
            sd.resp.c=sd(T1_avg*Mass,na.rm=TRUE),
            max.resp.c=max(T1_avg*Mass,na.rm=TRUE),
            min.resp.c=min(T1_avg*Mass,na.rm=TRUE))

resp_15C_SGGS_all_Uninf %>% # whole organism respiration for manuscript
  summarise(n=n(), # all populations pooled
            min.resp=min(T1_avg*Mass),max.resp=max(T1_avg*Mass),
            mean.resp.c=mean(T1_avg*Mass,na.rm=TRUE),
            sd.resp.c=sd(T1_avg*Mass,na.rm=TRUE),
            max.resp.c=max(T1_avg*Mass,na.rm=TRUE),
            min.resp.c=min(T1_avg*Mass,na.rm=TRUE))

resp_15C_SGGS_all_Inf %>% # whole organism respiration for manuscript
  summarise(n=n(), # all populations pooled
            mean.resp.c=mean(T1_avg*Mass,na.rm=TRUE),
            sd.resp.c=sd(T1_avg*Mass,na.rm=TRUE),
            max.resp.c=max(T1_avg*Mass,na.rm=TRUE),
            min.resp.c=min(T1_avg*Mass,na.rm=TRUE))


lm1 <- lm(y1 ~ T1_mass_conv, data=resp_15C_SGGS_all_Uninf)
summary(lm1) # 
lm2 <- lm(y1 ~ T1_mass_conv, data=resp_15C_SGGS_all_Inf)
summary(lm2) # 
lm3 <- lm(y1 ~ T1_mass_conv+Population+Infected, data=resp_15C_SGGS_all)
summary(lm3) # 


ggplot(resp_15C_SGGS_all_Uninf, aes(x = T1_mass_conv, y = y1, color = Population) ) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
  ggtitle("15C SGGS")

ggplot(resp_15C_SGGS_all, aes(x = T1_mass_conv, y = y1, 
                               color=Infected, shape=Population)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
  ggtitle("15C GG & SGGS All")

ggplot(resp_15C_SGGS_all, aes(x = T1_mass_conv, y = y1, 
                         color=Infected)) +
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = rev(wes_palette("Darjeeling1", n = 2)))+
  ggtitle("15C SGGS All")

# jpeg("analyses/respirometry/plot_15C_SGGS_Uninf.jpg")
# ggplot(resp_15C_SGGS_all_Uninf, aes(x = T1_mass_conv, y = y1, 
#                                color=Population)) +
#   geom_point()+theme_bw()+
#   theme(legend.position=c(0.9,0.1))+
#   geom_smooth(method = "lm", se = FALSE)+
#   scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
#   NULL
# dev.off()

# jpeg("analyses/respirometry/plot_15C_SGGS_Uninf_all.jpg")
# ggplot(resp_15C_SGGS_all_Uninf, 
#        aes(x = T1_mass_conv, y = y1)) +
#   geom_point()+theme_bw()+
#   geom_smooth(method = "lm",color="darkgrey")+
#   scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
#   ggtitle("15C All Populations (Uninfected)")
# dev.off()

resp_15C_SGGS_all_Uninf %>% group_by(Population) %>% summarize(max=max(Mass))

bypop <- ggplot(resp_15C_SGGS_all_Uninf, 
                aes(x = T1_mass_conv, y = y1, 
                               color=Population)) +
  geom_point()+theme_bw()+
  theme(legend.position=c(0.85,0.15))+
  geom_smooth(method = "lm", se = FALSE)+
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
  xlab("ln(mass)")+ylab("")+
  ylim(17.75,19.3)+
  NULL

allstickles <-ggplot(resp_15C_SGGS_all_Uninf, 
                     aes(x = T1_mass_conv, y = y1)) +
  geom_point()+theme_bw()+
  geom_smooth(method = "lm",color="darkgrey")+
  scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
  # ggtitle("15C All Populations (Uninfected)")
  xlab("ln(mass)")+ylab("")+
  NULL

bypop+allstickles

# jpeg("analyses/respirometry/empirical_resp.jpg",width=9.475,height=4.5,units="in",res=300)
# bypop+allstickles
# dev.off()

infuninf <- ggplot(resp_15C_SGGS_all, aes(x = T1_mass_conv, y = y1, 
                         color=Infected)) +
  geom_point()+theme_bw()+
  theme(legend.position=c(0.85,0.15))+
  geom_smooth(method = "lm")+
  scale_color_manual(values = rev(wes_palette("Darjeeling1", 
                                              n = 2)))+
  xlab("ln(mass)")+
  ylab(expression(paste("Metabolic rate: ln(Ie"^"E/kT",")")))+
  ylim(17.75,19.3)+
  NULL


# jpeg("analyses/respirometry/empirical_resp_all.jpg",width=14.2,height=4.5,units="in",res=300)
# bypop+allstickles+infuninf
# dev.off()

# jpeg("analyses/respirometry/empirical_resp_fig.jpg",width=9.5,height=4.5,units="in",res=300)
# bypop+infuninf
# dev.off()

resp_15C_all$WormWetMass[is.na(resp_15C_all$WormWetMass)]<-0
resp_15C_all$wormtofish <- (resp_15C_all$WetMass-resp_15C_all$WormWetMass)/resp_15C_all$WetMass

resps <- ggplot(resp_15C_SGGS_all,aes(x=Infected,y=y1,
                        fill=Infected,color=Infected))+
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = rev(wes_palette("Darjeeling1", 
                                              n = 2)))+
  scale_fill_manual(values = rev(wes_palette("Darjeeling1", 
                                              n = 2)))+
  geom_violin(alpha=0.3)+
  geom_boxplot(alpha=0.5)+
  geom_point(shape=21, position = position_jitterdodge(),color="gray40")+
  xlab("Infection Status")+ 
  scale_x_discrete(labels=c("Uninfected","Infected"))+
  ylab(expression(paste("Metabolic rate: ln(Ie"^"E/kT",")")))+
  ylim(17.75,19.3)+
  NULL
resps

# respsbypop <- ggplot(resp_15C_SGGS_all_Uninf,aes(x=Population,y=y1,
#                                  fill=Population,color=Population))+
#   theme_classic()+
#   theme(legend.position = "none")+
#   scale_color_manual(values = wes_palette("FantasticFox1", n = 3))+
#   scale_fill_manual(values = wes_palette("FantasticFox1", n = 3))+
#   geom_violin(alpha=0.3)+
#   geom_boxplot(alpha=0.5)+
#   geom_point(shape=21, position = position_jitterdodge(),color="gray40")+
#   xlab("Population")+ 
#   # scale_x_discrete(labels=c("Uninfected","Infected"))+
#   # ylab("")+
#   ylim(17.75,19.3)+
#   ylab(expression(paste("Metabolic rate: ln(Ie"^"E/kT",")")))+
#   NULL
# respsbypop


jpeg("analyses/respirometry/empirical_resp_Figure3.jpg",
     width=7.3,height=4.5,units="in",res=300)
resps+infuninf+
  plot_layout(widths = c(0.6,1))+
  plot_annotation(tag_levels = "A")
dev.off()

pdf("analyses/respirometry/empirical_resp_Figure3.pdf",
     width=7.3,height=4.5)
resps+infuninf+
  plot_layout(widths = c(0.6,1))+
  plot_annotation(tag_levels = "A")
dev.off()

# jpeg("analyses/respirometry/empirical_resp_suppfig.jpg",
#      width=7.3,height=4.5,units="in",res=300)
# respsbypop+bypop+
#   plot_layout(widths = c(0.8,1))+
#   plot_annotation(tag_levels = "A")
# dev.off()


resp_15C_SGGS_all_Uninf %>% 
  group_by(Population) %>% #each population separately
  summarise(n=n(),
            mean.wt=mean(Mass,na.rm=TRUE),
            sd.wt=sd(Mass,na.rm=TRUE),
            max.wt=max(Mass,na.rm=TRUE),
            min.wt=min(Mass,na.rm=TRUE))

resp_15C_SGGS_all_Uninf %>% 
  summarise(n=n(),
            mean.wt=mean(Mass,na.rm=TRUE),
            sd.wt=sd(Mass,na.rm=TRUE),
            max.wt=max(Mass,na.rm=TRUE),
            min.wt=min(Mass,na.rm=TRUE))

resp_15C_SGGS_all %>% group_by(Infected) %>% tally()

set.seed(2020)
resp_15C_SGGS_all_sub <- resp_15C_SGGS_all %>% group_by(Infected) %>% sample_n(12) 

lm3sub <- lm(y1 ~ T1_mass_conv+Population+Infected, 
             data=resp_15C_SGGS_all_sub)
summary(lm3sub) # 
summary(lm3)
