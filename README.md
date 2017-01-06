# NeuroSeg

NeuroSeg is an open-source MATLAB program used for automatic detecting and segmenting cells in two-photon calcium imaging data.

###Installation
1. Extract the ZIP file (or clone the GitHub repository) somewhere you can easily reach it.
2. Add the NeuroSeg/ folder to your path in MATLAB by e.g.,

 * using the "Set Path" dialog in MATLAB;
 * running the “addpath” function from your command window or startup script;
 * running the script GUI: NeuroSeg.m.


###Usage
####1. Load imaging data: 
The “Load Image” function is to load an image file (such as *.tiff file), and the “Load Video” function is to load a video file (e.g., *.avi file).
####2. Set the parameter values of “Min sigma”, “Max sigma” and “Min area”: 
The “Min sigma” is the value of the minim radius of cells; the “Max sigma” is the value of the maximum radius of cells; the “Min area” is the minimum area for selecting cells above a defined area. 
####3. Detect the location of cells: 
The “Seed Detection” function is to detect the location of cells, and the results are marked with red plus symbols.
####4. Segment the cells: 
The “Segmentation” function is to segment the boundaries of cells, the results are shown as red contours.
####5. Label ROIs manually: 
The “Draw ROI” function is to add a freehand ROI, and the “Add Circle” function is to add a circle ROI with the defined radius (“Radius” parameter). The “Delete ROI” function is to remove the unwanted ROI, and it can be done by moving the cursor in the region of the ROI and then clicking.
####6. Plot calcium traces, the ΔF/F value (only optional for video data): 
The “Plot DF/F” function is to calculate the relative calcium signals (ΔF/F) and plot the calcium traces extracted from each ROI. The “Save DF/F” function is to save calcium traces (ΔF/F) as a mat file. You can set the values of the frame rate (Hz) and the display offset for the ROIs before you plot the calcium traces.
####7. Save and load ROI data: 
The “Save ROI & Image” function is to save the information of the ROIs and the image into a mat file, and the “Load ROI & Image” is to load the mat file and show the ROIs in the image.
####8. Select a specific ROI: 
Click the index number of ROI in the “ROI no.” listbox, and then the contour of the ROI will flicker in green three times. So the specific ROI can be deleted.


*The software was tested on MATLAB R2014a version.

###Contact
Jiangheng Guan (735341369@qq.com)

Xiang Liao (xiang.liao@aliyun.com)





