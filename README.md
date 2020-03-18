# Alleviating confounding in spatio-temporal areal models
This repository contains the Supplementary Material and R code to fit with INLA the spatio-temporal models considered in the data analysis section of the work entitle _"Alleviating confounding in spatio-temporal areal models with an application on crimes against women in India"_ (Adin et al., 2020). It also contains the necessary functions to reproduce all the figures and tables of the present article.


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
Dowry deaths and socio-demographic covariates in 70 districts of Uttar Pradesh, India, during the period 2001-2014. The data is publically available online without any form of restriction or copyright.

- [**DowryDeaths_UttarPradesh.txt**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/DowryDeaths_UttarPradesh.txt)
  
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
	- **x5**: Murder rate. Number of murders per 100,000 inhabitants. Source: Open Government Data Platform India (https://data.gov.in).
	- **x6**: Burglary rate. Number of burglaries per 100,000 inhabitants. Source: Open Government Data Platform India (https://data.gov.in).


- [**Uttar_Pradesh_nb.graph**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/Uttar_Pradesh_nb.graph)
  
  An inla.graph object with the spatial neighborhood structure of the 70 districts of Uttar Pradesh.


- [**carto_up.shp**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/data/carto_up/)

  Shapefile containing the cartography of the 70 districts of Uttar Pradesh.


# R code
R code to fit with INLA (http://www.r-inla.org/) the spatio-temporal models considered in the data analysis section of the present work and code to reproduce all the figures and tables. All the R files are written by the authors of the paper using R version 3.6.2 (2019-12-12).

- [**DataAnalysis_INLA.R**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/R/DataAnalysis_INLA.R)

  This R script contains the necessary functions to replicate the model fitting with INLA of the spatio-temporal models considered in the data analysis section of the present work, using the dowry deaths data in Uttar Pradesh considered here, or replacing it with any other data with similar structure.
  
- [**Figures_and_Tables.R**](https://github.com/spatialstatisticsupna/Confounding_article/blob/master/R/Figures_and_Tables.R)
 
 
  This R script contains the necessary functions to reproduce all the figures and tables of the data analysis section of the present work. The fitted models with INLA and PQL can be download from


# References
Adin, A., Goicoa, T., Hodges, J.S., Schnell, P., and Ugarte, M.D. (2020). Alleviating confounding in spatio-temporal areal models with an application on crimes against women in India. https://arxiv.org/abs/2003.01946
