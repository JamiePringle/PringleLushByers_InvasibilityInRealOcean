# Warning, this Repo is not yet ready… Do not use until this warning has gone away, or at least do not complain if you have any problems doing so

# Code and figures Pringle, Lush & Byers

Interactive figures and code for the Pringle, Lush and Byers paper "Alongshore variation in currents make some regions of the coastal ocean more vulnerable to the introduction of non-native species with planktonic dispersal." More information will appear here once the manuscript makes it through the review process and we have adjusted this page to reflect any revisions. 

**Important:** More details of the python code, and the larval connectivity data it generates, can be found in the python section of the EZfate documentation. EZfate's GitHub page can be found [here](https://github.com/JamiePringle/EZfate) and the documentation for the python version of the code can be found [here.](https://jamiepringle.github.io/EZfate/EZfate_python_guide.html)

### Interactive figure 3

The interactive version of [Figure 3 in the paper is here](https://jamiepringle.github.io/PringleLushByers_InvasibilityInRealOcean/).

### Code to replicate the results of Section ****Discussion: Comparisons with observed introductions:_Carcinus maenas_ in the Canadian Maratimes after the mid 1980s****

The results in this paper can be replicated for other regions of the globe, other planktonic durations and larval vertical behaviors. The code to do so is described here. 

The specific files are configured for the section on _Carcinus maenas_ because it is less computationally demanding to reproduce than the results for all of coastal North and South America. To configure the code for a different problem, it is important to step through this guide in order. 

#### Preparation for running code

The programs below are all written in python, and a working python environment with the `numpy`, `matplotlib`, `cartopy`, `zarr`, `xarray`, `scikit-learn`, `numba`, `s3fs`, and `tdqm` packages must be installed. If you have the miniconda python environment installed (available [here](https://www.anaconda.com/docs/getting-started/miniconda/main)), you can create and enable an appropriate environment called LSRfind with the following command:
```
conda create -n LSRfind numpy ipython matplotlib cartopy zarr xarray scikit-learn numba s3fs tqdm
conda activate LSRfind
```

When the code release in this GitHub repository [LINK TO RELEASE] is downloaded and uncompressed, it should create the necessary empty directories that this code expects: `OSNdataDir`, `modelOutputNeutral` and `modelOutputRelativeFitness`. The programs should be run in a directory which has these directories within it. 

All of the data required to run the models has been created as part of the [EZfate project](https://github.com/JamiePringle/EZfate). It will automatically be downloaded by the `getEZfateFromOSN.py` module.  

#### Step One: configure the "habitat" with `00_makeConnectivityMatrices_trimByDistance.py` or `00_makeConnectivityMatrices_trimByDepth.py`

The "habitat" is defined as the collection of locations where a species can exist. It is defined on the same 1/12th of a degree grid as the ocean model used to parameterize the larval dispersal, the Mercator Ocean GLORYS model. More details on the computation of the Lagrangian particle tracks, including the available choices of vertical behavior, can be found at the [EZfate project web page.](https://github.com/JamiePringle/EZfate)

In all cases, the habitat is confined to a spatial extent with a polygon defined in the variable `regionPoly`. Furthermore, you should alter the loop which defines the variable `regionName` to only include the regions will be using – this avoids unnecessarily downloading large data files. The possible regions are `theAmericas`, `AsiaPacific`, `EuropeAfricaMiddleEast` and `Antarctica` as shown in the EZfate docmentation [here.](https://jamiepringle.github.io/EZfate/03_GetData_Subset_and_Combine.html)

You must also define the depth at which the larvae are released (from 1, 10, 20 and 40m depth) with the variable `depth` and the vertical behavior of the larvae (either fixed to a depth or drifting in all three dimensions) with the variable `vertBehavior`.  You can also specify what month(s) the larvae are released in the array `inMonths`.

If you want to define your habitat with water depth, use `00_makeConnectivityMatrices_trimByDepth.py` and define the variables `minHabDepth` and `maxHabDepth`. If you want to define the habitat with distance from the coast, use `00_makeConnectivityMatrices_trimByDistance.py` and define the variable `gridRadius`. 

When one or both of these codes are run, they will save connectivity data into a directory called `transposes`. This directory is created if it does not exist.  The connectivity data structures are defined [here.](https://jamiepringle.github.io/EZfate/EZfate_python_guide.html)

#### Step Two: Define the size of the introduction regions with `01_makeInitialIntroductionRanges.py`

In a perfect world with infinite computational resources, each habitat point would be analyzed individually. However, it is often more tractable to lump adjacent 1/12th of a degree habitat points together into a larger region, and then introduce novel species into the larger regions. There is a tradeoff involved – make the larger region too large, and spatial resolution is lost. Make it too small, and the models can take a very long time to run (especially those in which there can be a fitness difference between species). Some experimentation is required.  

To make the regions, and save the data defining the regions into the directory `initialConditions`, all that needs to be done is to define the variable `Nclose` to the size of the region in `01_makeInitalIntroductionRanges.py`, and also define the connectivity you wish to analyze in the variable `ConnectivityModelName`, and then run `01_makeInitalIntroductionRanges.py`. 

For the _Carcinus maenas_ model, `Nclose` was 5, for the full North and South American domain, it was 25. 

This code takes only a few moments to run, and then displays the resulting set of regions. The regions are saved in the directory `initialConditions` with the same name as the connectivity data defined in `ConnectivityModelName`.
