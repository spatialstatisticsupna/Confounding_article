# Alleviating confounding in spatio-temporal areal models
This repository contains the Supplementary Material and R code to fit with INLA the spatio-temporal models considered in the data analysis section of the work entitle _"Alleviating confounding in spatio-temporal areal models with an application on crimes against women in India"_ (Adin et al., 2020). It also contains the necessary functions to reproduce all the figures and tables of the present article.


## Table of contents

- [Data](#Data)
- [R code](#R-code)


# Data
Dowry deaths and socio-demographic covariates in 70 districts of Uttar Pradesh, India, during the period 2001-2014.

- **DowryDeaths_UttarPradesh.txt**:
	- **x0**: Political party of the Chief Minister ruling Uttar Pradesh during the study period: Bharatiya Janata Party (BJP) during 2001; Bahujan Samaj Party (BSP) during 2002-2003 and 2007-2011; Samajwadi Party (SP) during 2004-2006 and 2012-2014 (Source: https://www.mapsofindia.com/uttar-pradesh/chief-ministers.html or https://en.wikipedia.org/wiki/List_of_chief_ministers_of_Uttar_Pradesh)
	- **x1**: sex ratio. Number of females per 1000 males (Source: Office of the Registrar General and  Census Commissioner, India. (http://censusindia.gov.in)
	- **x2**: population density (People/Km2) (Source: Office of the Registrar General and  Census Commissioner, India. http://censusindia.gov.in)
	- **x3**: female literacy rate (Source: Office of the Registrar General and  Census Commissioner, India. (http://censusindia.gov.in)
	- **x4**: per capita income referenced to year 2004 (Source: Directorate of Economics And Statistics Government Of Uttar Pradesh. (http://updes.up.nic.in)
	- **x5**: number of murders per 100000 inhabitants (Source: Open Government Data Platform India. https://data.gov.in)
	- **x6**: number of burglaries per 100000 inhabitants (Source: Open Government Data Platform India. https://data.gov.in)
