# dendritic_distance
Compute the distance between pairs of spines captured in a single field of view for different types of distance - both dednritic path following distance, and euclidian "as the crow flies" distance." It is currently intended for planar images, but could be adapted to Z stacks and distance in 3d space in the future. Currently requires annotation of the dendrite, but certainly an algorithm could be implemented to replace this in the future. 


## Installation instructions

1. DD requires a  a recent Python distribution.

2. Open a terminal and type:

```bash
cd ..\documents\code #or similar as relevant for your machine
git clone git@github.com:GreggHeller1/dendritic_distance.git
cd dendritic_distance
conda env create -f environment.yml
Conda activate dendritic_distance
pip install read_roi

```

3.  should now be installed and the dendritic distance environment should still be activated. 
4. to run the algorithm, change the path in config.py to a data directory containing a tiff and corresponding ROIs, or contianing nested directories with tiffs and ROIs
```bash
cd path/to/dendritic_distance_installation
python dendritic_distance.py
```



## Usage instructions
### Input Files:
The code expects to find a directory or tree of nested directories where some the following files in each directory:
1. 1 file whose name includes the string "dend" that contains annotations of the dendrites and branch points
1. 1 file whose name includes the string ".zip" and not the string ".dend" that contains the annotations of rois for each spine head and corresponding nearest dendritic segment
1. A file whose name includes the string '.tif' that contains the average image from which the ROIs were generated (necessary only to produce QC images, not to compute the distance itself, but this is HIGHLY recomended as a way to validate both the annotations and the algorithm). If there are multiple filenames with '.tif' then the code will preferentially select a tif that includes the string 'stabalized.tif'

For more information on the expected annotations see annotations.md and ROI_selection.pdf

### Usage steps
1. Rename "default_config.py" to "config.py" (First time only)
1. Change the paths in "config.py" as desired
1. In a terminal with the conda environment active, run "python dendritic_distance.py"

### Output Files:
Each session directory should contain a subdirectory (named from "config.py") containing the following files.
N refers to the unmber of annotated spines. All distance matricies should be symmetric with 0s along the diagonal, with the ixjth entry representing the distance from the ith spine to the jth spine. If there is no dendritic path between the two spines (dendritic segments are not connected within the frame) then the distance will be nan. 
1. "dendritic_distance.csv" - CSV with NxN matrix of the distance between spines measured following the dendrite. The spine neck and head are not included (measurement is taken from the neck junction with the dendtire, as inferred by the closest distance to the spine ROI. 
1. "euclidian distance between Dendrite-Neck Junctions.csv" - CSV with NxN matrix, same as above (measures from neck junction to neck junction) except distance is measured "as the crow flies," euclidian. (Assumes same Z plane)
1. "euclidian distance between Spine Centers.csv" - same as above, but measures euclidan distance from the mean of each spine ROI
1. "euclidian distance between Spine Centers.csv" - same as above, but measures euclidan distance from the point on each spine ROI that is furthest from the dendrite segment
1. "seperated_by_branch_point.csv" - CSV with NxN binary matrix, 0 if the pair of spines is on the same dendritic segment and 1 otherwise. 
1. "annotated_dendrite.png" - png overlaying the ROIs of interest with the avereage image. Two pairs of spines are connected indicating the two pairs with the largest differential between "dendritic distance" and "euclidian distance" with the corresponding distances indicated in the legend. The connection is drawn from the means of the ROIs only for visualization purposes. The actual dendritic measurement is made using the neck junction. The euclidean distance used here is for the spine centers (means) because otherwise the euclidean distance can never be longer than the dendritic distance. 
1. "stem_stats.csv" - CSV containing computed statistics about the spines, mainly for QC purposes (long necks, residuals or angles far from 90 can be indicative of incorrect annotations)

Please see "demo_data" for examples

Some summary files are also produced, each with a datestring corresponding to when the code finished running
1. processed_check_directories_datestring.json - this contains a list all the directories that were sucessfully processed 
1. missing_annotation_check_directories_datestring.json - This file contains a list of directories that did not have the 3 necessary files to run the algorithm.
1. failed_list_check_directories_datestring.json - This file contains a list of the directories that failed to run the algorithm (possibly for other reasons than missing files, not having 3 dendritic segments associated with a branch point is a common one)
1. asymetric_list_check_directories_datestring.json - the distance matricies should all be symmetrical. If they are not, directories will be listed here. This file will also contain sessions where not all dendritic segments are connected and therefore Nans are found in the distance matricies. 

The remaining summaries each consist of a pair of files. 1 png containing histograms of the relevant statsitic, and 1 JSON with directories and values of the highest 5 and lowest 5 sessions for each stastic. 
1. Neck length - The distance from the mean of the spine ROI to the neck junction with the dendtite segment. Generally we are concerend about neck lengths that are too long indicating the spine was not mapped to the correct dendritic segment.
1. Neck Angle - The angle between the dendritig segment and the spine neck. A Neck far from 1.5 (radians) can indicate that the spine was not mapped to the correct dendritic segment
1. Dendrite Residual - The distance from the mean of the "dendrite ROI" to the nearest point on the nearest dendritic segment. Large Dendrite residuals can indicate that the same image was not used for both annotations, a dendritic segment was not annotated, or that the spine ROIs are not ordered properly (one may be missing causing an offset and spine ROIs to be inferred as dendrite ROIs
1. Spine Circumference - sum of all the outer segments in the spine ROI. If this is excessively large it could indicate that the wrong pixel/um mapping was used (perhaps a different objective or pixel size was chosen) or that the spine ROIs are not ordered properly (one may be missing causing an offset and dendrite ROIs to be inferred as spine ROIs which have a smaller cicumferance)


Please see "demo_data/summary_plots" for examples

# Quality Control (QC)

It is important that you QC the results of this code and do not trust it blindly. Here is my process for QC.
1. Look through the pngs in annotated images. The dendrite should be annotated appropriately, the spine ROIs should look like spines (not background or dendritic ROIs) and spines should be connected to the appropriate dendritic segment. Pay attention near branch poits. 
1. Look through the "missing_annotation" and "failed" JSONs. If there are sesssion directories, check them for what may have gone wrong, look for all the files, check the branch points for 3 dendritic segments, etc. (this may involve re-running the code on that directory alone to do some debugging.)
1. look at the statistics plots and Jsons. The histograms should give you a feel for what are extreme outliers. I normally closely inspect the sessions for the 5 "max_all" directories for each statistic. If I find errors in all of them it is necessarry to rerun the code on all of the images to regenerate the top 5. (Maybe I should provide afull list ordered so the user can choose when to stop?. "Max all" is typically the most useful, although pixel size errors would show up in the median statistic, and it can be good to check the minimum neck angles (should just take diff from pi/2 here...). 
1. Correct the ROIs if necessary and then re-run (the easiest way to do this can be to delete the dendritic_distance subfolder, and then re-run while leaving the flag "re-run = false" in the config this will skip andirectories that already contain the dendritic_distance subfolder). 

# Credit

This code was created for the surlab at MIT by Gregg Heller. It was created using images provided by Kyle Jenks, who also helped with inspiration and validation. 

