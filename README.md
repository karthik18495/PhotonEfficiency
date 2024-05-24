# What is in the data

The Fall 2018 data was analyzed and stored. Now, this datafile was produced by making histograms with various cuts from a really big tree that is put up in the farm. This was done using the script `Scripts/MakePlots.C`. This made the root file located in the folder `Data`

The data is split into various historgrams. Each histogram is Invariant mass of both reconstructed and expected (omega) particles. The difference is that reconstructed particles exactly have all the showers reconstructed photons versus the expected particles that has either 3 or 4 photons. Check out the presentation from [November 2023](OmegaPi0Efficiency.pdf) on how these photons are numbered and reconstructed. 

## The fitting script

From the trimmed data which has histograms of invariant mass a function of momentum and theta bins. We have to fit them with various difference functions. Once we fit them we then put it in the efficiency formula, which basically, compared how many omega's are reconstricted versus how many "ideally" expect. THis is run using the script `FitOmegaSpectrum.C`. The way to run this is simply do 

```
cd Scripts
root -q -b 'FitOmegaSpectrum.C("../Data/Fall2018_PhotonEff.root")'
```

This script goes over every histogram and fits them with 4 different functions. and makes a pdf. Additionally, it computes the yield after fitting for both the reconstructed and expected histograms and makes the efficiency plot, This may be too much to read and understand from the script but I can go over this.

## What would i do to start

PLease open up a root browser and load the file `Data/Fall2018_PhotonEff.root` 

```
root -l Data/Fall2018_PhotonEff.root
TBrowser a
```

This should open up a dialog box and browse each of the histogram. Browse through all the histograms that starts with the name `Omega14MassCut7_*`. This are the histogram of the Mass of reconstructed particles and the names `Omega14InvariantMassCut7_` should be ideal expectations.

Now, Also, check out the first 200 lines of the script `FitOmegaSpectrum.C` to see how the histogram 
# to do 
* add the location of the raw data files and instruction

