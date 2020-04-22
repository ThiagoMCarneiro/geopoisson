---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# GeoPoisson

<!-- badges: start -->
[![Travis build status](https://travis-ci.com/ThiagoMCarneiro/GeoPoisson.svg?branch=master)](https://travis-ci.com/ThiagoMCarneiro/GeoPoisson)
<!-- badges: end -->

The goal of GeoPoisson is to find the estimates from the Non-homogeneous Poisson model in Morales et al. (2016) and make a interpolation map for a neighboring region.

## Installation



You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ThiagoMCarneiro/GeoPoisson")
```
## Example


A short example of the application of the model:

Load the package:


>library(GeoPoisson)


Download the txt archives "X1.txt", "dados1.txt", and "loca1.txt" in data-raw folder.This data is made up by observations of 29 stations from 1980 to 2010. This rainfall data (daily observations in milimeters) from ANA (Agência Nacional de Águas - National Water Agency).

To find the parameters of the model:


>install.packages("mapproj")
>install.packages("maptools")
>install.packages("graphics")
>install.packages("ggmap")
>require(mapproj)
>require(maptools)
>require(graphics)
>require(MASS)
>require(Matrix)
>require(ggmap)


>X<-read.table("X1.txt",head=T)

>data<-read.table("dados1.txt",head=T)

>loca<-read.table("loca1.txt",head=T)


>resultado <- GeoPoisson(data,48.8,0.001,0.001,200000,150000,X,loca)


>Interpolationgrid(resultado,loca,data,size=19,intensity=0.75)

## basic example code


At the end of the code, one should get a table with estimates for the parameters and an Interpolationmap.





