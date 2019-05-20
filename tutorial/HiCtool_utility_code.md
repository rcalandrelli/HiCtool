# HiCtool utility code

We provide a series of utility functions included into [HiCtool_utilities.py](/scripts/HiCtool_utilities.py) that allow to work with HiCtool generated data in the Python environment, such as:

- Loading and saving HiCtool contact matrices (either in compressed and tab-separated format).
- Loading DI values or HMM states.
- Loading topological domain coordinates.

These functions will allow you to use HiCtool data for additional analyses, and eventually re-save to file your processed data in order to be plotted using HiCtool. You could even have contact matrices generated with other software, parse them into a numpy matrix format, and then save them using the saving functions in order to be plotted.

First, open your Python or iPython shell and execute the script:
```Python
execfile('HiCtool_utilities.py')
```






