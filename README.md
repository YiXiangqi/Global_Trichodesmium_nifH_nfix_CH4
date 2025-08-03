# Global_Trichodesmium_nifH_nfix_CH4
This repository contains code for analyzing *Trichodesmium*-associated methane production in oceanic environments, supporting the findings of our research(under submission). The work demonstrates the significant role of this nitrogen-fixing cyanobacterium in marine methane emissions, particularly in tropical ocean regions.

## Reproducibility
Processed data needed for reproducing this work is not deposited here because of file size limitation. A complete frozen snapshot of this research, including all code, processed data and key intermidiate results, needed for exact reproduction of our work, is permanently archived on [Figshare](http://doi.org/10.6084/m9.figshare.29815921):

> Note: The Figshare archive will not receive code updates, while this GitHub repository may be actively developed.

## Origianl data

- Field *Trichodesmium* nifH gene copies and nitrogen fixation: [Shao et al., 2023](http://doi.org/10.5194/essd-15-3673-2023)
- SST, SSS, surface nitrate, surface phosphate, MLD and minimum dissolved oxygen in 0-500m: [World Ocean Altas 2023](https://www.ncei.noaa.gov/products/world-ocean-atlas)
- PAR and Chl-*a*: [NASA Ocean Color](https://oceandata.sci.gsfc.nasa.gov/l3/)

## Coputational Environment
### Hardware
- Dell PowerEdge R940xa
- 160 CPU cores
- 3TB RAM

### Software Stack
- OS: Ubuntu 20.04 LTS
- R: version 4.4.2 (see renv.lock file for used packages)
- Matlab: R2024b and [m_map](https://www-old.eoas.ubc.ca/~rich/map.html) package (mainly used for global mapping figure)

## Contact
For questions regarding the analysis or data access:  
Xiangqi Yi - yixq@jmu.edu.cn  
Polar and Marine Research Institute, College of Harbor and Coastal Engineering, Jimei University