This script uses a modified version of R package "copynumber"

If you copy all the scripts in this folder to a new directory, ensure that your R working directory is set to this same location, and then run the main script source("SCNA_script.R") it should work

Included in the this folder is some example toy data so you can see the file formats and check it runs. The raw data comes from ASCAT copy number tools, but you should be able to
use this script also with facets output, you just will just need to change the code around lines 12-20, and 30-34 to work with the column names of the different tool

The output should be 3 PDFs written in that same directory, and then 1 PDF is outputted to the R environment. Please note the output will not match the figure from our paper, as only example toy data is shared here, not the full dataset


