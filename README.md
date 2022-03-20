# SDC2-tianlai

SKA Science Data Chanlleges(SDCs) is organized by SKA Observatory(SKAO), of which the purpose is to prepare for the data processing of SKA in the future.

This repository contains a python wrapper script for processing [SKA science data challenge 2(SDC2)](https://sdc2.astronomers.skatelescope.org/sdc2-challenge) data from one group of team @NAOC-Tianlai. The main data processing software for source finding and characterization is :tada: [SoFiA-2](https://github.com/SoFiA-Admin/SoFiA-2):tada:.


---
### Introduction

There are one templete SoFiA-2 parameters file ` sofia_ska.par` and two python files `util_func.py` and `run.py`.

- `sofia_ska.par` is the template parameter file of SoFiA-2, the detail can be found in https://github.com/SoFiA-Admin/SoFiA-2/wiki/SoFiA-2-Control-Parameters. 

- `util_func.py` contains most functions which will be applied. There are four main functions for four tasks.

  - Function `params_search()` in `util_func.py` is to explore the parameters space. You can set a list of discrete values of each parameter you try to turn, and each parameter combination will be a independent input for SoFiA-2 to process the  same development datacube of which the source catalog is known. Following the description of scoring in https://sdc2.astronomers.skatelescope.org/sdc2-challenge/scoring-code, a simple scoring function can be applied to each output catalog of SoFiA-2. There will be N(the number of total parameter combination) subdirectories in the root output directory with the name like *sofia_dev_params_[0, 1, …, N-1]*, and each subdirectory contains not only the standard output files from SoFiA-2, but also a parameter file, a `matched_source.txt`file and a emply file named `[matched number]_[total number founded]_[score]_sources`. The `matched_soruce.txt` file contains the information of sources whose positions have been matched with the true catalog. With the `multiprocessing` package, this step can be parallelized.
  
  - Function `params_search_summary()` is to make a summary of the parameter exploration results from function `params_search`. It will collect the simple statistical results of each parameter combination(i.e. the `[matched number]_[total number founded]_[score]_sources` file), and calculate more statistics such as *match rate*([number matched] / [total number founded]), the *corrected score*([score] - [unmatched number]), the *average score per source*([corrected score] / [matched number]), and it is also encouraged to take more statistics into account. Then sort the results by these statistics and take the "optimal" parameter combination. 
  
  - After getting the "optimal" input parameters according to your sorting and selection, with function `sofia2_pipe()`, you can apply this parameter combination to the target datacube which is assumed to have a consistent statistics as the develepment data. And the data division is also operated in this function.
  
  - Function `submit()` is to combine and format the output from function `sofia2_pipe()` for submission. 
  
- The above four main functions in `util_func.py` will be imported and called in `run.py`, you can specify arguments of each in this file.


### Requirements

- SoFiA-2
  > (Retest and make some minor changes for the latest version v2.4.1. Farther revision might still be needed if the ACSII format of SoFiA-2 output catalogue changes, mainly for the function `read_sofia_cat()` in *util_func.py*)
- python >= 3.6 (f-string is great!)
- pathos (replace the built-in *multiprocessing* package)
- numpy
- pandas
- astropy
