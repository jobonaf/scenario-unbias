
#library(tiff)
library(raster)
library(ncdf4) 

# BOXPLOT OF THE SCENARIO METHODS
scenario_boxplot <- function(var,state,dir_fairmode_data,dir_scen_nc,dir_out) {
  
  ######################################################################
  #- USAGE:
  ## var = variable PM25, O3, NO2
  ## state =c("Italia", "Germania", "Romania"," Belgio", "Francia", "Polonia" ) and so on - states of which you want the single boxplot
  ## dir_fairmode_data = directory in which the faimode excercise data are located
  ## dir_scen_nc =directory in which the the elaborated scenario data are located
  ## dir out= directory where you find the elaborated boxplot
  #
  # YOU NEED A FILE WITH THE ORDER THAT YOU WANT FOR THE SCENARIOS IN THE BLOXPLOT --> ROW: 65
  # YOU NEED A FILE WITH THE LOCATION OF THE STATION --> ROW: 71
  ######################################################################
  
  if (var=="PM25") {
    id<-"SURF_ug_PM25_rh50"
    id2<-"PM25_rh50"
  } else if (var=="NO2") {
    id<-"SURF_ug_NO2"
    id2<-"NO2"
  } else if (var=="O3") {
    id<-"SURF_ppb_O3"
    id2<-"O3"
  } else {
    stop(paste0(" *** no available variable: ",var," *** admitted vars: PM25, O3, NO2"))
  }

  
# LIST OF ORDERED CONSIDERED SCENARIO  
unbiased_name<- read.delim (paste0('../data/elenco_',var,'.txt'),header=F)
num_scen<-length(t(unbiased_name))

# FILE WITH STATION POINTS
matr<-read.csv("../data/staz_eu.csv")
d<-dim(matr)
  
# UPLOAD FAIRMODE DATA:
  #BaseCase
BCfile<-nc_open(paste0(dir_fairmode_data,'/BaseCase_Perturbed_Gridded/BaseCase_PERT_',id2,'_YEARLY.nc'))
BCdata<-ncvar_get(BCfile,id)
LON<-ncvar_get(BCfile,"lon")
LAT<-ncvar_get(BCfile,"lat")
DD<-dim(BCdata)
nc_close(BCfile)
  #Reference Points
SURFdata<-read.csv(paste0(dir_fairmode_data,'/BaseCase_Reference_Points/yearly_',id,'.csv'), sep = ',', header = TRUE)
  #Scenario Pert
SCENfile<-nc_open(paste0(dir_fairmode_data,'/Scenario_Perturbed_Gridded/SCEN_PERT_',id2,'_YEARLY.nc'))
SCENdata<-ncvar_get(SCENfile,id)
nc_close(SCENfile)
dd<-dim(SCENdata)
ref_case<- matrix(nrow=length(SURFdata$lat),ncol=3)

#Extract data on the point stations
for (k in seq(1:length(SURFdata$lat))) {
  l1<-which.min(abs(SURFdata$lon[k]-LON))
  l2<-which.min(abs(SURFdata$lat[k]-LAT))
  ref_case[k,1]<-BCdata[l1,l2]
  ref_case[k,2]<-SURFdata[k,3]
  ref_case[k,3]<-SCENdata[l1,l2]
  
}
colnames(ref_case)<-c("BASE-CASE","REFERENCE","SCENARIO")


# REDING THE FILE CONTENING THE UNBIAS SCENARIO
v<-matrix(nrow=DD[1]*DD[2], ncol=num_scen)

for (i in seq(1:num_scen)) {
  str_name<-unbiased_name[i,1]
  file_name<-paste0(dir_scen_nc,"/",var,"/Scen_ITAWG_",var,"_",str_name,"_CORR_YEARLY.nc")
  if (file.exists(file_name)) {
    ncfile<-nc_open(file_name)
    data_scen<-ncvar_get(ncfile,id)
    lat<-ncvar_get(ncfile,"lat")
    lon<-ncvar_get(ncfile,"lon")
    
    matr[1:d[1],d[2]+i]<-NA  
    colnames(matr)[d[2]+i]<-str_name
  
    # EXTRACT POINT CORRESPONDING TO POINT STATIONS
    for (j in seq(1:d[1])) {
      l1<-which.min(abs(matr$lon[j]-lon))
      l2<-which.min(abs(matr$lat[j]-lat))
      matr[j,d[2]+i]<-data_scen[l1,l2]
    }
    # FOR BLOXPLOT ON ALL THE GRIDDED DATA
    v[,i]<-as.vector(data_scen)
    
    nc_close(ncfile)
  } else {
    print(paste0(dir_scen_nc,"/",var,"/Scen_ITAWG_",var,"_",str_name,"_CORR_YEARLY.nc"))
    print (paste(str_name,": this method is not in the directory"))
    matr[1:d[1],d[2]+i]<-NA  
    colnames(matr)[d[2]+i]<-str_name
    v[,i]<-NA
  }

}

# BOXPLOT EUROPE - on station points
png(paste0(dir_out,"/boxplot_points_",var,"_EU.png"),width = 1200, height = 800)
par(mar = c(9, 5, 4, 2) + 0.1)
M<-max(apply(matr[,(d[2]+1):(d[2]+i)], 2, function(x) quantile(x, probs = 0.75,na.rm=T)))
boxplot(cbind(ref_case,matr[,(d[2]+1):(d[2]+i)]),col = c("yellow","magenta","green",rep("lightblue",num_scen)), border = c("yellow3","darkmagenta","darkgreen",rep("darkblue",num_scen)),
        outline = FALSE,
        ylim = c(0, 1.5*M), ylab="ug/m3", las=2, main= paste0(var," - Station Points - Europe"))
dev.off()


# BOXPLOT ON ALL EUROPE DOMAIN
colnames(v)<-colnames(matr)[(d[2]+1):(d[2]+i)]
png(paste0(dir_out,"/boxplot_ALL_",var,"_EU.png"),width = 1200, height = 800)
par(mar = c(9, 5, 4, 2) + 0.1)
M<-max(apply(v, 2, function(x) quantile(x, probs = 0.75,na.rm=T)))
box<-boxplot(v[,],col = "deepskyblue", border = rep("royalblue4",num_scen),
        outline = FALSE,
        ylim = c(0,1.2*M), ylab="ug/m3", las=2, main= paste0(var," - ALL Europe"))
dev.off()


#BOXPLOT ON SINGLE STATES - on station points
for (stato in state) {
  par(mar = c(9, 5, 4, 2) + 0.1)
  dati_stato<-subset(matr, matr[,5] == stato)
  print(paste("state:", stato))
  if (nrow(dati_stato) == 0) {
    print(paste("Not correct name for state:", stato))
    print("*** LIST OF ADMITTED STATE: ***")                                                                                                            
    print(unique(matr[, 5]))
    return()
  }
  png(paste0(dir_out,"/boxplot_points_",var,"_",stato,".png"),width = 1200, height = 800)
  ref_case_stato<-subset(ref_case, matr[,5] == stato)
  M<-max(apply(dati_stato[,(d[2]+1):(d[2]+i)], 2, function(x) quantile(x, probs = 0.75,na.rm=T)))
  box<-boxplot(cbind(ref_case_stato,dati_stato[,(d[2]+1):(d[2]+i)]),col = c("yellow","magenta","green",rep("grey",num_scen)), border = c("yellow3","darkmagenta","darkgreen",rep("black",num_scen)),
               outline = FALSE,
               ylim = c(0, 1.5*M), ylab="ug/m3", las=2, main= paste0(var," - Station Points ",stato))
  dev.off()
  
}

}
