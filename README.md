# Colored Objects Recognition in OpenCV
Machine vision (image processing) project for detection and counting colored objects (Haribo jellies) on scene restricted to area delimited by 3 types of markers (OpenCV, C++)


  Given the sequence of images representing Haribo jellies on the piece of paper with the area delimited by 3 types of black markers, program firstly had to filter the images (by means of adaptive thresholding, dilation, erosion, shape approximation, contours reduction, etc.) so that there are only 4 contours left representing markers' contours. Then the cartesian coordinates of the contours had to be converted to polar coordinates so that they could be next sorted (quicksort algorithm) in the clockwise order. After that the suitable warp perspective transform was applied depending on the detected type of marker. Now when we extracted from the input image desired ROI (region of interest) we have to convert the initial image from RGB space to HSV space (hue, saturation, value) so we could count how many jellies in each color are visible on the scene (we had 6 different jelly colors: white, yellow, orange, light red, dark red, green). In order to accomplish it, InRange operator had to be used with the experimentally matched HSV ranges representing each of this colors. Next some morphological operations were needed to enhence the quality of detected contours. Number of jellies in each color was estimated by sorting the remaining  contours by they area size. After that, the median of these area sizes was calculated and the middle element was selected. Each area is then divided by this middle element and rounded to the nearest integer value. The same operation is carried out for all detected contour's areas and the results are then added up. The arisen value represent the number of jellies with the specific color. The same operation is repeated for all 6 colors and the result for each color is saved to .txt output file.


PURPOSE:
  The goal of the project is to count colorful haribo jellies on the various scenes inside particular area delimited by 3 types of markers. Every scene is captured from different perspectives, called frames. We assume that:
- Longer edge of the marker takes over 1/3 width of a picture
- The angle between camera optical axis and surface isn't greater than 30 degrees.

INPUT FILES:
  All folders have to be in the same catalog as the executable file. The names of the pictures are located in the file 'pictures_names.txt', inside 'pictures_names' folder. Pictures are placed in 'pictures' folder. 
Format of the text file containing pictures' names is as following:
- scene_XXX_frame_AAA.jpg
- scene_XXX_frame_BBB.jpg
- scene_YYY_frame_AAA.jpg
- scene_YYY_frame_BBB.jpg
XXX, YYY are the following scene numbers. AAA, BBB are the number of the following frames. The numbers always start with '001' and always has 3 digits.

OUTPUT FILES:
  The result text file is located in 'results' folder, which is in the same catalog as other folders.
The name of the executable file: Surname_Name.exe
The name of output text file: Surname_Name.txt
Numbers of gummy bears for every scene have to be in a separate line. They should be listed in a specified order:
dark_red, bright_red, green, orange, white, yellow<CR><LF>
and so on for every scene.

  To sum up, the folder structure is as follows:
  - main folder
  ----- images
  --------- scene_001_grame_001.jpg
  --------- scene_001_grame_002.jpg
  ...
  ---- pictures_names
  --------- pictures_names.txt
  ---- results
  ---------- Surname_Name.txt
  ---- Surname_Name.exe














