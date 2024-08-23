# condensate-fusion-aspect-ratio

****** README.md *****

# Title: Aspect Ratio Image Analysis

### Description
This code uses the Python (v3.10.11) and packages Scikit-image, NumPy, Pandas, Scipy 
and MatplotLib to identify single fusion events of two condensates in pre-processed 
image files. Each fusion event is fit to the equation AR = (AR0 - 1)*e^(-t/tau) +1,  
where AR is the aspect ratio of the fusing droplets, AR0 is the aspect ratio at 
time 0, t is time and tau is the time constant of decay. Aspect ratio, approximate 
length scale of the fusion event, and tau are recorded and output. This process is
repeated across multiple events and the resulting data is compiled into a csv file. 

### Code Contributors
Written by: Jenna Tom \
Deniz Lab, www.scripps.edu/deniz

Corresponding Author: Prof. Ashok Deniz\
Institution: The Scripps Research Institute\
ORCID: 0000-0003-2819-4049 \
GitHub: https://github.com/Deniz-Lab \
Email: deniz@scripps.edu

Published and used in:\
Ravi Chawla*, Jenna K. A. Tom*, Tumara Boyd, Nicholas Tu, Tanxi Bai, Danielle A. Grotjahn, 
Donghyun Park, Ashok A. Deniz, Lisa R. Racki. "Reentrant DNA shells tune polyphosphate 
condensate size."  

For use, please cite the above publication.

## Requirements (Version used)
- Python (v 3.10.11)

- Packages: NumPy, Pandas, Matplotlib, Scikit-image, Scipy
- Custom code: aspectratio.py, parsekey_wrangling.py (developed by Emily Bentley)
- Data files: "PolyPFusion_tifFiles" (available on GitHub: https://github.com/Deniz-Lab/condensate-fusion-aspect-ratio)


## File Overview

- aspectratio.py: handles single droplet fusion events. 
                  Measures aspect ratio, returns fit parameters, displays segmentation
- parsekey_wrangling.py: translates titles into dictionary of experimental parameters
										 


## HOW TO USE: 
1. Open PolyP_AspectRatio_v1.ipynb in a program that can run a Jupyter Notebook file. 
2. Follow notebook instructions for running code. You should, at the end of running all 
blocks, have a compiled csv output. Graphs and figure exports and displays can be added 
and/or adjusted according to notebook instructions. 


