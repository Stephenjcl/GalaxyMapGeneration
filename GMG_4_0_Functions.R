####Naming Functions####
nam.rand <- function(coordx = 6000, coordy = 4000, Sysnum = 64, n=1,clock = F,gen = F) {
  
  rngs <- vector(length = n)
  
  m <- 2 ** 32
  a <- 1103515245 - coordx
  c <- 12345 + coordy
  
  d <- (Sysnum*((coordx/2)+coordy/2))
  if(clock == T) {d <- as.numeric(Sys.time()) * 1000}
  
  
  for (i in 1:n) {
    d <- (a * d + c) %% m
    rngs[i] <- d / m
    if (gen == T){rngs[i] <- round((rngs[i]*63)+1,0)}
  }
  
  return(rngs)
}

truffleshuffle <- function(coordx,coordy,Sysnum,val=1){
  coordx = coordx + Sysnum;
  coordy = coordy + coordx
  coordx = (coordx * (2 ^ 3))
  while(is.infinite(coordx)==T){coordx <- coordy/Sysnum}; while(is.infinite(coordy)==T){coordy <- coordx/Sysnum}; while(is.infinite(Sysnum)==T){Sysnum <- coordx/coordy}
  coordx = coordx + coordy
  coordy = (coordy * (2 ^ 5))
  coordy = coordy + coordx
  coordy = (coordy * (2 ^ 4))
  coordx = (coordx * (2 ^ Sysnum))
  coordx = coordx + coordy
  while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
  while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
  if(val==1){return(coordx)}
  if(val==2){return(coordy)}
  if(val==3){return(Sysnum)}
}


####Syllable Construction####
#Done procedurally, since we have already set phoneme frequencies. This will inherently increase how natural the language sounds as the odd pairings will be reduced. 

GetSystemName <- function(coordx, coordy, Sysnum) {
  #Here Sysnum will itself be altered by the galactic coordinates. Ideally this would prevent the trend of similarily positioned stars having similar names.
  Sysnum <- nam.rand(coordx, coordy, Sysnum)
  
  coordx = coordx + Sysnum
  coordy = coordy + coordx
  coordx = (coordx * (2 ^ 3))
  coordx = coordx + coordy
  coordy = (coordy * (2 ^ 5))
  coordy = coordy + coordx
  coordy = (coordy * (2 ^ 4))
  coordx = (coordx * (2 ^ Sysnum))
  coordx = coordx + coordy
  while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
  while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
  
  sylcount <- syllablechance[bitAnd(coordx, 31) + 1] #Decide number of syllables
  if (sylcount == 0 | is.na(sylcount) == TRUE) {sylcount <- 2} #Set syllables to 2 if we error out and get 0 as a result
  
  if ((nam.rand(coordx,coordy,Sysnum)) > 0.7) {
    #70% Chance of utilizing Eyan A names
    tmp1 <- as.character(EyanB[nam.rand(coordx,coordy,Sysnum)*2502,1]) #Word selection
  } else {tmp1 <- as.character(EyanA[nam.rand(coordx,coordy,Sysnum)*2502,1])} #Word selection
  coordx <- truffleshuffle(coordx,coordy,Sysnum,1)
  
  #Second word
  coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
  coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
  coordy <- truffleshuffle(coordx, coordy, Sysnum,2)
  coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
  
  if ((nam.rand(coordx,coordy,Sysnum)) > 0.7) {
    #70% Chance of utilizing Eyan A names
    tmp2 <- as.character(EyanB[nam.rand(coordx,coordy,Sysnum)*2502,1]) #Word selection
  } else {tmp2 <- as.character(EyanA[nam.rand(coordx,coordy,Sysnum)*2502,1])} #Word selection}
  
  while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
  while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
  
  if(sylcount == 1){tmp2 <- ""}
  
  sysname <- paste0(tmp1, tmp2)
  sysname <- paste(toupper(substr(sysname, 1, 1)), substr(sysname, 2, nchar(sysname)), sep="")
  tmp1 <- NULL; tmp2<- NULL;
  #print(sysname)
  return(sysname)
}

####Renaming for 9 Nation Style Map####
RenameSystem <- function(df) {
  
  pb2 <- txtProgressBar(min = 0, max = sum(df$Name != "", na.rm=TRUE), style = 3); c <- 0;
  
  res <- character(nrow(df)) #Optimization by preallocating space for results
  
  for (i in 1:nrow(df)){
    if((df$Name[[i]]) != ""){
      coordx <- df$x[[i]]; coordy <- df$y[[i]]; Sysnum <- df$z[[i]]; col <- df$clusterNum[[i]] #Get baseline data
      subdict <- dict[,col]
      subdict <- as.data.table(as.character(subdict))
      subdict <- subdict[!apply(subdict == "", 1, all),]
      maxvalue <- as.numeric(nrow(subdict))-1
      #Here is where the loop will return to modify Sysnum again in the hopes of producing a new unique name
      repeat {
        sysname <- NULL #This will empty sysname if we are looping. Missed this at first, ended out with 4000 character or more names.
        Sysnum <- nam.rand(coordx, coordy, Sysnum)
        
        #Extra randomizing
        coordx = coordx + Sysnum
        coordy = coordy + coordx
        coordx = (coordx * (2 ^ 3))
        coordx = coordx + coordy
        coordy = (coordy * (2 ^ 5))
        coordy = coordy + coordx
        coordy = (coordy * (2 ^ 4))
        coordx = (coordx * (2 ^ Sysnum))
        coordx = coordx + coordy
        while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
        while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
        
        sylcount <- syllablechance[bitAnd(coordx, 31) + 1] #Decide number of syllables
        if (sylcount == 0 | is.na(sylcount) == TRUE) {sylcount <- 2} #Set syllables to 2 if we error out and get 0 as a result
        
        tmp1 <- as.character(subdict[(1+nam.rand(coordx,coordy,Sysnum)*maxvalue),1]) #Word selection
        
        #Second word
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        coordy <- truffleshuffle(coordx, coordy, Sysnum,2)
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        
        tmp2 <- as.character(subdict[(1+nam.rand(coordx,coordy,Sysnum)*maxvalue),1]) #Word selection
        
        while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
        while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
        
        if(sylcount == 1){tmp2 <- ""}
        
        sysname <- paste0(tmp1, tmp2)
        #sysname <- paste(toupper(substr(sysname, 1, 1)), substr(sysname, 2, nchar(sysname)), sep="")
        sysname <- stringr::str_to_title(sysname)
        tmp1 <- NULL; tmp2<- NULL;
        
        #Test for uniqueness of name
        if(sysname %in% res == F){break}
      }
      
      res[i] <- sysname
      
      
      sysname <- NULL
      c <- c +1; setTxtProgressBar(pb2,c)
    } #End While loop
  }
  close(pb2)
  df$Name <- res
  return(df)
}

FillSystemName <- function(df) {
  
  pb2 <- txtProgressBar(min = 0, max = nrow(df), style = 3); c <- 0;
  
  res <- character(nrow(df)) #Optimization by preallocating space for results
  check <- df$Name
  
  for (i in 1:nrow(df)){
    if((df$Name[[i]]) != ""){res[i] <- df$Name[[i]]}
    if((df$Name[[i]]) == ""){
      coordx <- df$x[[i]]; coordy <- df$y[[i]]; Sysnum <- df$z[[i]]; col <- df$clusterNum[[i]] #Get baseline data
      subdict <- dict[,col]
      subdict <- as.data.table(as.character(subdict))
      subdict <- subdict[!apply(subdict == "", 1, all),]
      maxvalue <- as.numeric(nrow(subdict))-1
      #Here is where the loop will return to modify Sysnum again in the hopes of producing a new unique name
      repeat {
        sysname <- NULL #This will empty sysname if we are looping. Missed this at first, ended out with 4000 character or more names.
        Sysnum <- nam.rand(coordx, coordy, Sysnum)
        
        #Extra randomizing
        coordx = coordx + Sysnum
        coordy = coordy + coordx
        coordx = (coordx * (2 ^ 3))
        coordx = coordx + coordy
        coordy = (coordy * (2 ^ 5))
        coordy = coordy + coordx
        coordy = (coordy * (2 ^ 4))
        coordx = (coordx * (2 ^ Sysnum))
        coordx = coordx + coordy
        while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
        while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
        
        sylcount <- syllablechance[bitAnd(coordx, 31) + 1] #Decide number of syllables
        if (sylcount == 0 | is.na(sylcount) == TRUE) {sylcount <- 2} #Set syllables to 2 if we error out and get 0 as a result
        
        tmp1 <- as.character(subdict[(1+nam.rand(coordx,coordy,Sysnum)*maxvalue),1]) #Word selection
        
        #Second word
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        coordy <- truffleshuffle(coordx, coordy, Sysnum,2)
        coordx <- truffleshuffle(coordx, coordy, Sysnum,1)
        
        tmp2 <- as.character(subdict[(1+nam.rand(coordx,coordy,Sysnum)*maxvalue),1]) #Word selection
        
        while (coordx >= (2^31)-1) {coordx <- coordx/((2^33)+coordy)}
        while (coordy >= (2^31)-1) {coordy <- coordy/((2^33)+coordx)}
        
        if(sylcount == 1){tmp2 <- ""}
        
        sysname <- paste0(tmp1, tmp2)
        #sysname <- paste(toupper(substr(sysname, 1, 1)), substr(sysname, 2, nchar(sysname)), sep="")
        sysname <- stringr::str_to_title(sysname)
        tmp1 <- NULL; tmp2<- NULL;
        
        #Test for uniqueness of name
        if(sysname %in% check == F){break}
      }
      
      res[i] <- sysname
      
      
      sysname <- NULL
      c <- c +1; setTxtProgressBar(pb2,c)
    } #End While loop
  }
  close(pb2)
  df$FullName <- res
  return(df)
}

####Generator Functions####
#Acquire number of stars for a given sector
generate_sectorhr <- function(coordx=4000, coordy=4000, galaxyScale = 0) {
  if (coordx > 0x1fff || coordy > 0x1fff) {
    return(0)
    break
  }
  
  pixelval <- (coordx)+2*bitAnd(coordy,0x1fc0);
  
  p1 <- TheMilkyWayHR[coordx, coordy];    #Current center
  p3 <- TheMilkyWayHR[coordx, (coordy+1)];    #Next row
  p4 <- TheMilkyWayHR[(coordx+1), (coordy+1)];    #Next row, next column
  p2 <- TheMilkyWayHR[(coordx+1), (coordy)];    #Next column
  
  coordx <- bitAnd((coordx * 0x200), 0x7e00);
  coordy <- bitAnd((coordy * 0x200), 0x7e00);
  
  ebx <- (p2-p1)*coordx + (p3-p1)*coordy;    #-0x0FB0400...0x0FB0400
  esi <- (coordx*coordy)/(2^15);        #0..0x7C08
  edi <- p4 - p3 - p2 + p1;
  esi <- esi*edi;
  ebx <- ebx + esi;
  ebx <- ebx + ebx;
  p1 <- p1 * (2^16);                #p1 is now VVVV.VVVV.0000.0000
  ebx <- ebx + p1;
  ecx <- ebx/(2^8);
  if (galaxyScale < 16) {
    ebx = coordx + ecx;
    eax =coordx * coordy;
    eax = eax/(2^15);
    ebx = bitOr(ebx, eax);
    ebx = ebx/(2^5);
    ebx = bitAnd(ebx,0x7f);
    if (ebx == 0) {
      ecx <- 0
    } else {
      eax = SystemDensity[ebx];
      if (exists("galaxyScale")==TRUE) {
        edx <- 16-galaxyScale;
        edx <- edx*eax;
        edx <- edx/(2^5);
        eax <- 65535-edx;
        ecx <- ecx*eax;
        ecx <- ecx/(2^16);
      } else{
        ecx <- ecx*eax;
        ecx <- ecx/(2^16);
      }
    }
  }
  p1 <- ecx;
  p1 <- p1/(2^10);
  p1 <- round((p1*290.6587),0) #Convert output to a number of starts between 0-64
  return(p1)
}

#Randomizing function
rotate_some <- function(){
  Tmp1 = bitOr((SystemParam_0*(2)), (SystemParam_0/(2^10)))
  Tmp2 = SystemParam_0 + SystemParam_1;
  Tmp1 = Tmp1 + Tmp2;
  SystemParam_0 <<- Tmp1;
  SystemParam_1 <<- bitOr((Tmp2*(2)), (Tmp2/(2^10)));
  if (SystemParam_0 >= 100000000) {SystemParam_0 <<- SystemParam_0/10000000}
  if (SystemParam_1 >= 100000000) {SystemParam_1 <<- SystemParam_1/10000001}
}

#Returns basic star information (number, position, color, type) for given coordinate
Shuffle_Coordinates <- function(coordx = 4000, coordy=4000, number_of_systems = 50,trim=T) {
  dataframe <- NULL
  x <- c("x", "y", "z", "Stardesc", "Multiple", "Colour"); dataframe <- as.data.frame(matrix(nrow=64,ncol=length(x)));   names(dataframe) <- x; x <- NULL
  SystemParam_0 <<- bitShiftR(coordx,16) + coordy;
  SystemParam_1 <<- bitShiftR(coordy,16) + coordx;
  
  rotate_some();
  rotate_some();
  rotate_some();
  
  for (i in 1:number_of_systems) {
    rotate_some()
    z <- (((bitAnd(SystemParam_0,0xFF0000)/(2^16))));
    if (nam.rand(coordx,coordy,number_of_systems)>0.5){z <- z*(-1)}
    dataframe[i,3] = z; #z-axis
    dataframe[i,2] = bitAnd(SystemParam_1,0x0000FF)*.92; #y-axis
    piv <- bitAnd(SystemParam_0,0x0001FE)/16
    for (j in 1:piv){rotate_some()}
    dataframe[i,1] = (bitAnd(SystemParam_0,0x0001FE)/2)-3; #x-axis
    dataframe[i,5] = StarChance_Multiples[bitAnd(SystemParam_1,0x1f)+1]; #Multi-star
    dataframe[i,4] = StarChance_Type[bitAnd((SystemParam_1/(2^16)), 0x1f)+1]; #Star Type
    plc <- as.numeric(dataframe[i,4]) + 1
    dataframe[i,6] = ColorForStar[plc]; #Star Color
    dataframe[i,4] = StarDesc[plc]
  }
  if (trim == T) {dataframe <- na.omit(dataframe)}
  return(dataframe)
}

#Returns the same as Shuffle_Coordinates, but with star names attached as well.
sectorfullgen <- function(coordx, coordy){
  #Initialize variables
  coords_t <- NULL
  x <- c("x", "y", "z", "Stardesc", "Multiple", "Colour"); coords_t <- as.data.frame(matrix(nrow=64,ncol=length(x)));   names(coords_t) <- x; x <- NULL
  x <- coordx
  y <- coordy
  n <- as.numeric(generate_sectorhr(x,y)) #Generate number of systems in this sector
  coords_t <- Shuffle_Coordinates(x,y,n)
  
  #The Naming stuff (doesn't have it's own function yet)
  namestorage <- as.data.frame(matrix(nrow=as.numeric(nrow(coords_t)), ncol=1))
  namestorage$V1 <- as.character(namestorage$V1)
  row <- 0
  for (i in 1:nrow(namestorage)){
    row <- row + 1
    row <- as.numeric(row)
    namestorage[i,1] <- GetSystemName(x,y,row)
  }
  
  #Check for mono character names and rerun
  namestorage <- within(namestorage, V1[nchar(namestorage$V1) == 1] <- GetSystemName(x,y,65)) #Replaces all cells in namestorage column V1 which have a character length of 1 and replaces them with a name generated on system position "65" (which is typically impossible). This only works if the is only 1 star in the system with a monocharacter name.
  #Apply the names to the cells
  coords_t$Name <- ""
  for (i in 1:nrow(namestorage)){
    coords_t[i,7] <- namestorage[i,1]
  }
  return(coords_t)
}

#Returns Star information for a given range of the galaxy. Planets = T will also append data for the planets of interest around said stars. Garden decides if you would like to filter said data.
snapgaldata <- function(reach = 4, planets=F, garden = F){
  #Loop to grab many stars data
  expanse <- reach #Value about the centre to travel
  stretch <- expanse/2
  #Min-max calculation
  llx <- 318-stretch
  lly <- 404-stretch
  ulx <- 318+stretch
  uly <- 404+stretch
  
  #Data frame to hold all the things
  galaxysnapshot <- NULL
  temp1 <- NULL
  print("Connecting to database...")
  brag <- ((reach+1)^2)*32
  brag <- paste0("Recruiting sector information for approximately ",brag," star systems. Please stand by.")
  if (Sys.info()["sysname"]=="Darwin"){system("say -v Karen -r 200 System data incoming. Stand by.")}
  print(brag)
  pb <- txtProgressBar(min = 0, max = ((reach+1)^2), style = 3)
  c <- 0
  for (x in llx:ulx){
    for (y in lly:uly){
      if (planets == T){
        temp1 <- sectorfullscan(x,y,garden)
      } else {
        temp1 <- sectorfullgen(x,y)
      }
      temp1$GalacticX <- (x*250)+temp1$x
      temp1$GalacticY <- (y*250)+temp1$y
      temp1 <- as.data.frame(temp1)
      galaxysnapshot <- suppressWarnings(bind_rows(galaxysnapshot, temp1))
      c <- c + 1
      setTxtProgressBar(pb,c)
    }
  }
  close(pb)
  brag <- (nrow(galaxysnapshot))
  brag <- paste0("Data acquired for ",brag," star systems with suitability for human presence.")
  brag2 <- paste0("say -v Karen -r 200 Download complete. ",(nrow(galaxysnapshot)), " suitable planets found.")
  if (Sys.info()["sysname"]=="Darwin"){system(brag2)}
  print(brag)
  
  galaxysnapshot <- galaxysnapshot[!grepl("No Garden Worlds Found", galaxysnapshot[,1]),]
  return(galaxysnapshot)
}

#For a given sector outputted by sectorfullgen, returns planetary info for that star.
planetscan <- function(starmap,coordx,coordy,garden=F){
  stdata <- NULL
  stdata <- matrix(nrow=nrow(starmap),ncol=8)
  stdata[,1] <- as.character(stdata[,1])
  
  for (s in 1:nrow(stdata)){ #Find number of planets around a given star
    stdata[s,1] <- as.numeric((round((nam.rand(coordx,coordy,s)*4),0))+1)
  }
  stdata[,2:8] <- ""
  
  for (i in (1:nrow(stdata))){
    if (stdata[i,1] > 0){
      rnd1 <- i*coordx
      for (n in 1:(stdata[i,1])){
        rnd2 <- n*coordy
        stdata[i,2] <- paste0(stdata[i,2], Atmosphere[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        rnd1 <- truffleshuffle(rnd1,rnd2,1,1) #Randomization shuffle
        
        stdata[i,3] <- paste0(stdata[i,3],Temperature[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        rnd1 <- rnd1*(nam.rand(rnd1,rnd2,1))
        
        stdata[i,4] <- paste0(stdata[i,4],Biosphere[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        rnd1 <- rnd1/(nam.rand(rnd1,rnd2,1))
        
        stdata[i,5] <- paste0(stdata[i,5],Population[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        rnd1 <- rnd1*(nam.rand(rnd1,rnd2,1))
        
        stdata[i,6] <- paste0(stdata[i,6],TechLevel[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        rnd1 <- rnd1/(nam.rand(rnd1,rnd2,1))
        stdata[i,7] <- paste0(stdata[i,7],WorldTags[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
        nd1 <- rnd1*(nam.rand(rnd1,rnd2,1))
        stdata[i,8] <- paste0(stdata[i,8],WorldTags[[((round((nam.rand(rnd1,rnd2,1)*10),0))+1)]], sep = "|", collapse = "|")
      }
    }
  }
  colnames(stdata) <- c("Planet Count", "Atmosphere","Temperature","Biosphere","Population","TechLevel","World Tag 1", "World Tag 2")
  #Added these re-classing areas to ensure clean merging later
  stdata[,1] <- as.numeric(stdata[,1])
  stdata[,2:8] <- as.character(stdata[,2:8])
  starmap <- cbind(starmap,stdata)
  
  #Below here will truncate the data to only show systems with potential "garden worlds" i.e., atmosphere, temp, and biosphere condusive to life. Not perfect however!
  if (garden == T){
    #Pull only rows which contain the key words for the various planet states that support life.
    garden <- starmap[grep("Breathable Mix", starmap$Atmosphere), ]
    garden <- garden[grep("Temperate", garden$Temperature), ]
    garden <- garden[grep('Human-miscible biosphere', garden$Biosphere), ]
    
    #Truly parse the data now to only produce the rows which contain not only these states, but these states in the correct combination
    #EJECT if there is nothing in garden
    if(nrow(garden)==0){return()}
    truelife <- garden
    #Prepare the dataframe by first splitting the columns apart on a per-planet basis
    suppressWarnings(truelife <- separate(truelife, Atmosphere, c("Atm1", "Atm2", "Atm3", "Atm4", "Atm5"), sep = "\\|", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn"))
    suppressWarnings(truelife <- separate(truelife, Biosphere, c("Bio1", "Bio2", "Bio3", "Bio4", "Bio5"), sep = "\\|", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn"))
    #Setup the key rows
    truelife$GardenRating <- 0
    
    for (i in 1:nrow(truelife)){
      n <- matrix(0, 1, 5)
      n[1,1]<- ifelse((grepl("Breathable Mix", truelife$Atm1)) && (grepl("Human-miscible biosphere", truelife$Bio1)), +1, 0)
      n[1,2] <- ifelse((grepl("Breathable Mix", truelife$Atm2)) && (grepl("Human-miscible biosphere", truelife$Bio2)), +1, 0)
      n[1,3] <- ifelse((grepl("Breathable Mix", truelife$Atm3)) && (grepl("Human-miscible biosphere", truelife$Bio3)), +1, 0)
      n[1,4] <- ifelse((grepl("Breathable Mix", truelife$Atm4)) && (grepl("Human-miscible biosphere", truelife$Bio4)), +1, 0)
      n[1,5] <- ifelse((grepl("Breathable Mix", truelife$Atm5)) && (grepl("Human-miscible biosphere", truelife$Bio5)), +1, 0)
      truelife[i,24] <- sum(n)
      n <- NULL
    }
    
    garden$Rating <- truelife$GardenRating
    garden<-garden[!(garden$Rating==0),]
    garden$Rating <- NULL
    starmap <- garden
    
  }
  return(starmap)
}

sectorfullscan <- function(coordx,coordy,garden=F){
  stardata <- sectorfullgen(coordx,coordy)
  #Grab planet info. Garden True will filter the creation by only viable planets.
  stardata <- planetscan(stardata,coordx,coordy,garden)
  return(stardata)
}

starimages <- function(df){
  #Time to get some star images.
  img <- NULL; mask <- NULL; sun <- NULL
  sun <- image_read("8k_sun-2.jpg")
  
  #Produce mask which will be used to find circular star texture
  png(tf <- tempfile(fileext = ".png"), 600, 600)
  par(mar = rep(0,4), yaxs="i", xaxs="i")
  plot(0, type = "n", ylim = c(0,1), xlim=c(0,1), axes=F, xlab=NA, ylab=NA)
  plotrix::draw.circle(.5,0.5,.5, col="black", border= "black")
  dev.off()
  mask <- image_read(tf)
  
  #Produce image for each star
  pb <- txtProgressBar(min = 0, max = (nrow(df)), style = 3)
  c <- 0
  for (i in 1:nrow(df)){
    #Find random crop of the sun to produce a star from.
    rndcrp1 <- runif(1, 0,(4096-850)) #Find random x-coordinate
    rndcrp2 <- runif(1, 0,(2048-850)) #Find random y-coordinate
    crploc <- paste0("800x800+",rndcrp1,"+",rndcrp2) #Crop to those points
    img <- image_crop(sun, crploc)
    img <- image_convert(img, format = "png")
    
    #Use mask to make image circular and retrieve
    img <- image_composite(mask, img, "minus") #Mask to circle with all black border
    img <- image_composite(mask, img, "plus") #Remask to remove the black border
    img <- image_rotate(img, runif(1,0,360)) #Rotate Randomly
    
    #Produce image for each star.
    img <- image_convert(img, colorspace = "gray") #Make B&W
    img <- image_convert(img, colorspace = "sRGB") #Recolor
    img <- image_trim(img) #Remove excess whitespace
    img <- image_transparent(img, "#ffffff") #Make pure white transparent. Since we've made our source image 1 RGB off of pure white, the star won't be affected.
    overlay <- image_colorize(img,100,df[i,6]) #Produce a colour overlay
    img <- image_composite(overlay,img, "blend") #Colorize with the stars color value
    
    #Produce a filename
    fname <- paste0(getwd(),"/Star_Images")
    if(dir.exists(fname)==F){dir.create(fname)}
    fname <- paste0(fname,"/","X",df[i,8],"-","Y",df[i,9],":",df[i,7],".png")
    
    image_write(img, path = fname, format = "png")
    c <- c + 1
    setTxtProgressBar(pb,c)
  }
  close(pb)
  print("Example image from stellar data on right:")
  return(img)
}

jumpvector <- function(snapgaldata = NULL, ps = NULL, test =F){
  print("Calculating 3-Dimensional vectors for lightspeed travel...")
  #calctime <-0.0002*(nrow(snapgaldata))^2.0059 #New version calculates based on number of rows being calculated in a power function. Only 2.3% error on average.
  #calctime <- round(calctime/60,1)
  #calctime <- paste0("This will take approximately ",calctime," minutes. Please stand by.")
  #print(calctime)
  if(test==T){
  start_time <- Sys.time()
  print(start_time)
  
  
  #The combn part was the slowest step. Using BioConductor's version will reportedly cut the time to a tenth!
  ps <- data.frame(t(apply(combnPrim(seq_len(nrow(snapgaldata)), 2), 2,
                           function(x) c(snapgaldata[x, ]$GalacticX, snapgaldata[x, ]$GalacticY, snapgaldata[x,]$z))))
  } else {
  
  ps <- data.frame(t(apply(combn(seq_len(nrow(snapgaldata)), 2), 2,
                           function(x) c(snapgaldata[x, ]$GalacticX, snapgaldata[x, ]$GalacticY, snapgaldata[x,]$z))))
  
  }
  print("Combinatorial Matrix Completed.")
  end_time <- Sys.time()
  timeelapsed <- end_time - start_time
  print(timeelapsed)
  
  #3D Vector
  ps$Vector <- sqrt((sqrt((ps$X2-ps$X1)^2+(ps$X4-ps$X3)^2)^2)+((ps$X6-ps$X5)^2))
  print("Completed.")
  
  if(test==T){
  end_time <- Sys.time()
  timeelapsed <- end_time - start_time
  print(timeelapsed)
  }
  
  return(ps)
}

nationmapfullmap <- function(snapgaldata, civilianrange=50, militaryrange=75,jumprange=150.5,ljumprange=300.5,print=T,size=24, save=F, filename=NULL, gardenworlddata=NULL, ps = NULL, savecsv = T){
  #Timing the code
  start_time <- Sys.time()
  print(start_time)
  #Adding a progress bar.
  pb <- txtProgressBar(min = 0, max = 5, style = 3); c <- 0;
  
  sf <- (max(snapgaldata$GalacticX)-min(snapgaldata$GalacticX))*(0.0003654971) #Scaling factor
  
  #ggplot solution
  if (is.null(ps)==T){ps <- jumpvector(snapgaldata)}
  
  print("Finding routes by engine class...")
  civilianpaths <- filter(ps, Vector <=civilianrange)
  militarypaths <- filter(ps, Vector <=militaryrange); militarypaths <- filter(militarypaths, Vector > civilianrange) #This prevents double lining paths that the military could clearly transport through since they have better engines.
  jumppaths <- filter(ps, Vector <= jumprange); jumppaths <- filter(jumppaths, Vector > jumprange-1)
  ljumppaths <- filter(ps, Vector <= ljumprange); ljumppaths <- filter(ljumppaths, Vector > ljumprange-1)
  
  #Remove non-garden worlds
  if (is.null(gardenworlddata)==F){
    print("Pruning dataset to contain only Garden worlds...")
    c <- c +1; setTxtProgressBar(pb,c)
    gardenworlddata$Key <- 1
    suppressMessages(x <- left_join(snapgaldata, gardenworlddata))
    x <- x %>% mutate(Name = ifelse(is.na(Key)==T, "", Name))
    x$Key <- NULL
    snapgaldata <- x
  }
  
  #Cluster Analysis
  print("Detecting current national borders...")
  c <- c +1; setTxtProgressBar(pb,c)
  sgdc <- nationclusterhier(snapgaldata)
  mdfsdc <- sgdc
  mdfsdc$clusterNum <- as.factor(mdfsdc$clusterNum)
  levelprint <- paste0("Plotting travel map for all nations of sector selected...")
  c <- c +1; setTxtProgressBar(pb,c)
  print(levelprint)
  
  #Add colour vector 
  mapcolour <- c('#ffe119', '#4363d8', '#f58231', '#fabebe', '#e6beff', '#800000', '#000075', '#a9a9a9', "#3cb44b")
  key <- data.frame(mapcolour = mapcolour, clusterNum = seq(1, nlevels(as.factor(mdfsdc$clusterNum)), by = 1))
  nm <- "mapcolour"
  mdfsdc['mapcolour'] <- lapply(nm, function(x) key[[x]][match(mdfsdc$clusterNum, key$clusterNum)])
  
  #Rename everything 
  mdfsdc <- RenameSystem(mdfsdc)
  
  #Scaling through the magic of math, to find the right size setting.
  borderboundary <- ((militaryrange*504.092)/10)/1.00375
  
  #Unify to two geom_point setups
  jumpmap <- ggplot(mdfsdc, aes(x = GalacticX, y = GalacticY, color = mapcolour)) + 
    geom_point(size = I(sf*0.2), shape = I(1)) + 
    #geom_point(size = I(borderboundary*sf*0.003), shape = I(16), alpha = 0.05) + #Working border attempt
    #stat_density2d(aes(fill = mapcolour, colour=mapcolour), alpha = 0.05, geom = "polygon", contour = F) + 
    # stat_density2d(aes(fill = mapcolour, colour=mapcolour), alpha = 0.05, geom = "polygon", contour = F) + 
    
    geom_segment(data = civilianpaths, mapping = aes(x = X1, xend = X2, y = X3,yend = X4), alpha = 0.3, linetype = 1, colour = "darkgreen", size = (sf*.75)) +
    geom_segment(data = militarypaths, mapping = aes(x = X1, xend = X2, y = X3,yend = X4), alpha = 0.2, linetype = 1, colour = "red", size = (sf*.75)) +
    geom_segment(data = jumppaths, mapping = aes(x = X1, xend = X2, y = X3,yend = X4), alpha = 0.3, linetype = 1, colour = "grey", size = (sf*.75)) +
    geom_segment(data = ljumppaths, mapping = aes(x = X1, xend = X2, y = X3,yend = X4), alpha = 0.3, linetype = 1, colour = "grey", size = (sf*.75)) +
    theme_bw() + 
    geom_text_repel(label=mdfsdc$Name, size=(sf*1.8),nudge_x = (sf*10),nudge_y = -(sf*3),segment.alpha = 0.6, segment.size = (sf*.25), box.padding = .125, color = "black") + 
    theme(axis.ticks=element_blank(),
          axis.text=element_blank(),
          plot.margin=grid::unit(c(0,0,0,0),"cm"),
          axis.title=element_blank(),
          legend.position = "none",
          panel.grid = element_blank()) + 
    labs(x = NULL, y = NULL)
  
  print("Saving jump-map to data drive... This will take several minutes, and is the slowest step. Standby.")
  c <- c +1; setTxtProgressBar(pb,c)
  
  if(is.null(filename == T)) {
    filename <- paste0(min(mdfsdc$GalacticX),",",min(mdfsdc$GalacticY),"-",max(mdfsdc$GalacticX),",",max(mdfsdc$GalacticY)," Sector Jump Map", ".png")}
  ggsave(filename, plot = jumpmap, width = size, height = size, units = "in", dpi = 300,device="png")
  print("Completed.")
  c <- c +1; setTxtProgressBar(pb,c)
  close(pb)
  end_time <- Sys.time()
  timeelapsed <- end_time - start_time
  print(timeelapsed)
  beep(1)
  if (savecsv == T) {
    mdfsdc <- FillSystemName(mdfsdc)
    filename <- substr(filename, 1, nchar(filename)-4)
    filename <- paste0(filename, ".csv")
    write.csv(mdfsdc,file=filename, row.names = F, fileEncoding="UTF-8")
  }
  return(jumpmap)
}

####Cluster Analysis for Different Cultures
nationclusterhier <- function(df, number = 9){
  set.seed(19930612)
  cdf <- df
  cdf <- na.omit(cdf)
  cdf$GalacticX <- scale(cdf$GalacticX); cdf$GalacticY <- scale(cdf$GalacticY); cdf$z <- scale(cdf$z)
  d <- dist(cdf[,c(16,17)], method = "euclidean")
  hc1 <- hclust(d, method = "complete")
  sub_grp <- cutree(hc1, k = number)
  clusterNum <- sub_grp
  df <- cbind(df, clusterNum = clusterNum)
  return(df)
}

generateMap <- function(rangevalue = 14, outputname = "WallMap_9_Nations.png"){
  #Initial data
  wall <- snapgaldata(rangevalue, planets=T)
  wallc<- (snapgaldata(rangevalue, planets=T, garden=T))
  jumpdata <- jumpvector(wall)

  jmpp <- nationmapfullmap(wall, gardenworlddata = wallc, ps = jumpdata, size = 24,save = T, filename = outputname, savecsv=T)
}

####Data Reading Functions####
readGalaxy <- function(filename = "MilkyWayScaled.png"){
  x <- readPNG(filename) #Read in the milky way map
  y <- rgb(x[,,1], x[,,2], x[,,3]) #Organize the channels
  yg <- desaturate(y) #Desaturate the image to colour saturation values as a hex character
  yn <- col2rgb(yg)[1, ]/255 #Convert Hex values to /numeric/ saturation values
  dim(y) <- dim(yg) <- dim(yn) <- dim(x)[1:2] #Rearrange back into single matrix
  TheMilkyWayHR <- yn; y <- NULL; yg <- NULL; yn <- NULL
  return(TheMilkyWayHR)
}

gRinstall <- function(){
  source("http://bioconductor.org/biocLite.R")
  biocLite(c("graph","RBGL","Rgraphviz"))
  
  install.packages("gRbase", dependencies=TRUE)
  # install.packages("gRain", dependencies=TRUE)
  # install.packages("gRim", dependencies=TRUE)
}