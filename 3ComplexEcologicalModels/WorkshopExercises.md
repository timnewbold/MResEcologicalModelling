# Session 3: Complex Ecological Models

This session will be based entirely around the Madinlgey general ecosystem model (Harfoot et al., 2014). Using the MadingleyR R package, you will download the model, run and plot basic simulations of the model, try altering the parameters of the model, and run scenarios of human impacts in the model.

## Setting up the model

You need to install an R package that works with NetCDF files (the format of many of the output files of the Madingley Model). Because the installer for this package is available only from my Dropbox, the installation process is a bit complicated.

First, you need to create a directory for downloading the installer, for example:

```R
installationDir <- "C:/Users/Tim/Documents/"
```

Now download the installer to this directory:

```R
download.file(url = "https://dl.dropboxusercontent.com/u/23999497/sds_0.1-5.zip",
              destfile = paste(installationDir,"sds_0.1-5.zip",sep=""))
```

Now install the package:

```R
install.packages(paste(installationDir,"sds_0.1-5.zip",sep=""), 
                 repos = NULL, type = "win.binary")
```

Now you need to install the R package that interfaces with the Madingley Model, and the R package for plotting Madingley simulations.

```R
library(devtools)
install_github("timnewbold/MadingleyR",subdir = "MadingleyR")
install_github("timnewbold/MadingleyPlots",subdir = "MadingleyPlots")

library(MadingleyR)
library(MadingleyPlots)
```

Now, define a directory on your computer where all of the code, data and eventual simulations will be stored, for example:

```R
codeDir <- "C:/Temp"
```

Now you can download the Madingley code and necessary data to this location (this can take a while, depending on the speed of your internet connection):

```R
DownloadRelease(codeDir = codeDir)
GetMadingleyData(codeDir = codeDir)
```

## Exercise 1: Running some basic simulations

In this exercise, we will run and plot some basic simulations from the Madingley Model.

To begin with, let's run a simulation for a single location.

The first step is to initialize the model with the parameters we want:

```R
init <- MadingleyInitialisation(duration = 50,
                                specific_locations = list(latitudes=2,longitudes=33))
params <- MadingleyParameters()
```

The duration=50 parameter means that the model will run for 50 years, and the specific_locations parameter defines the latitude and longitude of the location we wish to simulate. The empty call to MadingleyParameters means that default values for all of the parameters will be used (see <a href="https://github.com/timnewbold/MadingleyR/blob/master/EcologicalParameters.md">here</a> for the definitions and default values of the parameters, as described in Harfoot et al., 2014).

Having initialized the model, we can now run the specified simulation:

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
```

Finally, we can plot the time series of simulated biomass densities across the whole duration of this simulation:

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

You can see that plants (shown by the dark green line) have a higher biomass density than herbivores (light green), omnivores (purple) or carnivores (red), which is expected. But the estimated biomasses from a single simulation of a model as complex as the Madingley Model can be quite variable. Therefore, it is better to run multiple simulations, and then to calculate average estimated biomasses and confidence limits.

To do this we need to re-initialize the Madingley Model, specifying that the simulations be run in parallel (to make the most of all of the cores in your computer).

```R
init <- MadingleyInitialisation(duration = 50,parallel_simulations = "yes",
                                specific_locations = list(latitudes=2,longitudes=33))
```

Now when we run the model, we specify numSims=5 to tell the model to run 5 replicate simulations.

```R
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params,numSims = 5)
```
When we plot the model this time, we will specify plotConfidence=TRUE to show the variation in estimated biomasses around the mean:

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity",
             plotConfidence = TRUE)
```

As well as plotting the simulated biomass densities across the whole length of the simulation, we can also plot properties of the ecosystem at the end of the simulation. Useful properties of the ecosystem include the relative biomasses of different trophic levels, which can be plotted as a biomass pyramid, or the relationship between the body mass of organisms and their abundance density (empirical evidence leads us to expect a negative power-law relationship between mass and density). There are functions for making both of these figures in MadingleyPlots (PlotPyramids and PlotMassDensity). We can plot a biomass pyramid from the simulation we just ran:

```R
PlotPyramids(resultsDir = madingleySim$outputDir,plotName = "BiomassPyramids")
``` 

Plotting the relationship between body mass and abundance density requires running the Madingley Model with a high level of output detail. This makes the simulations rather slow, so we won't do this right now, but if you want to try running this later, you just need to add output_detail="high" to the call to the MadingleyInitialisation function.

Now, let's try running a simulation for two different locations (the original cell in Uganda and a cell with much drier climatic conditions in northern Algeria), and compare the results. We will switch back to running just a single simulation in the interests of time, but you could try running an ensemble of simulations later if you are interested.

```R
init <- MadingleyInitialisation(duration = 50,parallel_cells = "yes",
                                specific_locations = list(latitudes=c(2,35),longitudes=c(33,1)))
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
```

You can plot a map of the cells simulated using the following command:

```R
PlotCells(resultsDir = madingleySim$outputDir,map = "Africa")
```

Finally, plot the time series of estimated biomasses as before.

```R
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

You will see that this ecosystem is estimated to have lower biomasses overall, especially of carnivores (the relatively low abundance of carnivores in ecosystems with low plant biomass matches theoretical expectations - see e.g. Post, 2004).

## Exercise 2: Running simulations with altered ecological parameters

So far, we have been running simulations using the default ecological parameters. But it can be interesting to assess the effect on simulated dynamics of changing certain parameters.

Remember from yesterday's lecture that functional responses for a single predator or herbivore organism are composed of two main parameters: attack rates and handling times (the Madingley Model also contains parameters describing variation in attack rates and handling times across organisms, as a function of their body mass). Let's test the effect of increasing the handling time for all predator organisms in the model (the <a href="https://github.com/timnewbold/MadingleyR/blob/master/EcologicalParameters.md">default</a> handling-time scalar is 0.5; we will increase this by a factor of 4 to 2.0):

```R
init <- MadingleyInitialisation(duration = 50,
                                specific_locations = list(latitudes=2,longitudes=33))
params <- MadingleyParameters(pred_terr_handlingscalar = 2.0)
madingleySim <- RunMadingley(codeDir = codeDir,init = init,params = params)
PlotTemporal(resultsDir = madingleySim$outputDir,plotName = "BiomassDensity")
```

Any parameters not referred to in your call to MadingleyParameters will assume their default values.

You will see that there is now a much smaller biomass of carnivores and omnivores (because they can only predate much more slowly), and consequently a higher biomass of herbivores and a lower biomass of plants.

Now try instead increasing the handling times for  herbivory. The default value for the terrestrial handling time scalar is 0.7, so try increasing it to 2.8.

## References

* Harfoot, M.B.J., Newbold, T., Tittensor, D.P., Emmott, S., Hutton, J., Lyutsarev, V., Smith, M.J., Scharlemann, J.P.W. & Purves, D.W. (2014). Emergent global patterns of ecosystem structure and function from a mechanistic general ecosystem model. <i>PLoS Biology</i> <b>12</b>: e1001841.
* Post, D.M. (2004). The long and short of food-chain length. <i>Trends in Ecology & Evolution</i> <b>17</b>: 269-277.