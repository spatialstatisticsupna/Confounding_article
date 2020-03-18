# Alleviating confounding in spatio-temporal areal models
This repository contains the Supplementary Material and R code to fit with INLA the spatio-temporal models considered in the data analysis section of the work entitle _"Alleviating confounding in spatio-temporal areal models with an application on crimes against women in India"_ (Adin et al., 2020). It also contains the necessary functions to reproduce all the figures and tables of the present article.


## Table of contents

- [Data](#Data)
- [R code](#R-code)


# Data
Dowry deaths and socio-demographic covariates in 70 districts of Uttar Pradesh, India, during the period 2001-2014. The data is publically available online without any form of restriction or copyright.

- File name: [**DowryDeaths_UttarPradesh.txt**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/DowryDeaths_UttarPradesh.txt)
  
  This .txt file contains a data set with the following variables:
	- _dist_: Districts.
	- _year_: Year (from 2001 to 2014).
	- _state_: Uttar Pradesh.
	- _obs_: Number of dowry deaths.
	- _pop_linear_: Female population between 15 and 49 years (linear interpolation).
	- **x1**: Sex ratio. Number of women per 1,000 men. Source: Office of the Registrar General and Census Commissioner, India (http://censusindia.gov.in).
	- **x2**: Population density (people/km2). Source: Office of the Registrar General and Census Commissioner, India (http://censusindia.gov.in).
	- **x3**: Female literacy rate. Office of the Registrar General and Census Commissioner, India (http://censusindia.gov.in).
	- **x4**: Per capita income referenced to year 2004. Source: Directorate of Economics and Statistics Government of Uttar Pradesh  (http://updes.up.nic.in).
	- **x5**: Murder rate. Number of murders per 100,000 inhabitants. Source: Open Government Data Platform India https://data.gov.in).
	- **x6**: Burglary rate. Number of burglaries per 100,000 inhabitants. Source: Open Government Data Platform India https://data.gov.in).


- File name: [**Uttar_Pradesh_nb.graph**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/Uttar_Pradesh_nb.graph)
  
  An inla.graph object with the spatial neighborhood structure of the 70 districts of Uttar Pradesh.


- File name: [**carto_up.zip**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/carto_up.zip)

  Zipfile containing the cartography of the 70 districts of Uttar Pradesh (in shapefile format).
