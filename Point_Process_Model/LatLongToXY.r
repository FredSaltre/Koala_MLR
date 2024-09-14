     transfo.fct <- function(Lon,Lat)
      {
tab <- data.frame(Latitude = Lat, Longitude = Lon) 
	#La premi�re �tape est de crer un SpatialPointsDataFrame 
        #(�quivalent d'un shape sous R) � partir du tableau
	 library(rgdal)
	coords <- as.data.frame(cbind(Lon, Lat)) 
        #petite subtilit�, la premi�re colonne contient les Longitudes!
	SP.temp <- SpatialPointsDataFrame(coords = coords ,data=as.data.frame(rep(1,length(Lon))) ,proj4string = CRS("+init=epsg:4326"))
		#comme dans un shape on renseigne dans la fonction 
                #les coordonn�es des points (argument coords), 
                #la table attributaire (argument data, non indispensable) 
                #et surtout le syt�me de coordonn�es (arguement proj4string). 
                #Ici le syst�me g�ographique est le WGS84, 
                #soit on �crit la formule compl�te: CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), 
                #soit on �crit le code init: CRS("+init=epsg:4326").
		
	attributes(SP.temp) #l'objet cr�e est de la classe S4, il s'agit d'un syst�me de listes imbriqu�es.
	coordinates(SP.temp) 
          #les coordonn�es sont toujours en Longitude Latitude, 
          #pour les transformer il faut utiliser la fonction spTransform
	
	SP.proj <-  spTransform(SP.temp, CRS("+init=epsg:3034")) 
            #avec pour arguments le SpatialPointsDataFrame cr�e pr�c�dement et 
            #le syst�me de projection choisi. Dans l'exemple ce doit �tre de l'ETRS89.
	coordinates(SP.proj) #maintenant les coordonn�es sont projet�es et tu peux les r�cup�rer.
	
	tab$X <- coordinates(SP.proj)[, 1]
	tab$Y <- coordinates(SP.proj)[, 2]
	return(tab)
}

invtransfo.fct <- function(x,y)
{
 library(rgdal)
# Source data
	coords <- as.data.frame(cbind(x, y)) 
	SP.temp <- SpatialPointsDataFrame(coords =coords ,
     data=as.data.frame(rep(1,length(y))),
     proj4string = CRS("+init=epsg:3034"))
	SP.proj <-  spTransform(SP.temp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) 
ww <- coordinates(SP.proj)
return(ww)
}

