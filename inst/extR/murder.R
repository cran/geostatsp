# https://data.torontopolice.on.ca/datasets/TorontoPS::homicides-open-data-asr-rc-tbl-002/about
myFile = tempfile(fileext='.csv')
download.file(
'https://opendata.arcgis.com/api/v3/datasets/d96bf5b67c1c49879f354dad51cf81f9_0/downloads/data?format=csv&spatialRefId=3857&where=1%3D1',
myFile)

x1 = read.csv(myFile)

# no age and sex info.  

x2 = x1[,grep("^X|^Y|^LAT|^LONG|HOOD|ID$|DOY|DIVISION", names(x1), invert=TRUE)]
library(terra)
x = vect(as.matrix(x1[,c("LONG_WGS84", "LAT_WGS84")]), 
	crs = mapmisc::crsLL,
	atts = x2)

load("~/research/diseasemapping/pkg/geostatsp/inst/extdata/oldmurder.Rdata")
names(murder@attributes)

# https://homicidecanada.com/toronto-ontario-2023-homicide-victim-list/
# https://globalnews.ca/tag/toronto-homicide/
#https://www.instagram.com/p/C2NN_dUOYmZ/?hl=fr
x1$url = x1$age = x1$sex = NULL

x1[x1$OCC_YEAR==2023 & x1$OCC_MONTH=='October', ]


x1[x1$EVENT_UNIQUE_ID == 'GO-20232320141', 'url'] = 'https://globalnews.ca/news/10007531/toronto-fatal-stabbing-sheppard-avenue-wilmington/'
x1[x1$EVENT_UNIQUE_ID == 'GO-20232320141', 'age'] = 53
x1[x1$EVENT_UNIQUE_ID == 'GO-20232320141', 'sex'] = 'M'

x1[x1$OCC_YEAR==2023 & x1$OCC_MONTH=='September', ]

theId = 'GO-20232227404' 
x1[x1$EVENT_UNIQUE_ID == theId, 'url'] = 'https://www.cp24.com/news/2-killed-2-injured-in-north-etobicoke-shooting-1.6575349'
x1[x1$EVENT_UNIQUE_ID == theId, 'age'] = c(25, NA)
x1[x1$EVENT_UNIQUE_ID == theId, 'sex'] = 'M'

theId = 'GO-20232573487'
x1[x1$EVENT_UNIQUE_ID == theId, 'url'] = 'https://www.tps.ca/organizational-chart/specialized-operations-command/detective-operations/investigative-services/homicide/case/58/2023/'
x1[x1$EVENT_UNIQUE_ID == theId, 'age'] = 24
x1[x1$EVENT_UNIQUE_ID == theId, 'sex'] = 'M'


toAdd = list(
	'GO-20232935410' = list(url = 'https://www.cp24.com/news/1-man-dead-1-man-with-serious-injuries-after-overnight-shooting-in-north-york-1.6699701',
		age = 26, sex = 'M'),
	'GO-20232980418' = list(url = 'https://homicidecanada.com/toronto-ontario-2023-homicide-victim-list/',
		age = 41, sex = 'M'),
	'GO-20232954061' = list(url = 'https://www.cbc.ca/news/canada/toronto/homicide-toronto-investigation-1.7071672',
		age = 51, sex = "M"),
	'GO-20231309309' = list(url = 'https://www.instagram.com/p/CzZRRLqOOyH/?hl=fr',
		age = 53, sex = 'M'),
	'GO-20232400648' = list(url = 'https://www.instagram.com/p/CzZPCR2ush2/?hl=fr',
		age = 57, sex = 'F'),
	'GO-20232234677' = list(url='https://www.instagram.com/p/CyTcVfauIOE/?hl=fr',
		age = 21, sex = 'M'),
	'GO-20232320141' = list(url = 'https://homicidecanada.com/toronto-homicide-54-victorio-adriatico-charged-in-relation-to-the-fatal-stabbing-of-enrique-vinluan/',
		age =53, sex = 'M'),
	'GO-20232234677' = list(url = 'https://homicidecanada.com/toronto-homicide-53-police-investigate-the-fatal-stabbing-of-joyous-magdirila-in-north-york/',
		age = 23, sex= 'M')
)

