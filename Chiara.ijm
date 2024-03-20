// README
// This macro was designed to help Chiara Micchelli determine what threshold value to input
// into Imaris to get standardized segmentation of mitochondria, apicoplasts, and nuclei so
// that she can determine how the distance between these objects is affected by various
// treatments.

// This macro subtracts a somewhat ambiguous "background" level from each channel in order to
// improve Otsu thresholding.  The background level is determined by measuring half of the 
// average pixel value of the top X percent of pixels, with X being determined by user input.
// It has been empirically determined that the range of best values for X appears to be between
// 0.1 - 1%, but should be thoroughly tested by the end user.

// For the macro to function correctly, the following must be true:
// 1) Images should be three channel, in .czi format, and stored in the same folder.
// 2) Each channel should be 16 bit. If other bit depths are to be used, the code must be modified.
// 3) The order of channels must be Red (channel 1), Blue (channel 2), Green (channel 3).
// 3) Images should contain only one parasite. Images with multiple parasites should be cropped.

// This macro has two outputs - first, a series of tif files with the calculated background level
// subtracted from each channel; and second, a CSV file called "Results" with the subtracted
// background values and the calculated Otsu thresholds for each image and channel.  An IMPORTANT
// NOTE: the Otsu thresholds are calculated for the background subtracted image!  If you want to
// apply the thresholding to the original image, you must add the background subtracted value
// to the Otsu threshold and use that number.

// If you want to change the threshold method, change the word "Otsu" in line 150.

// Code starts here.
setBatchMode(true);

//Let user determine threshold levels for each channel
red = 0;
blue = 0;
green = 0;
Dialog.create("Set Interesting Pixel Level in Percent Form (ex, for top 1%, enter '1') AND directory of images.");
Dialog.addNumber("Red Threshold Level:", red);
Dialog.addNumber("Blue Threshold Level:", blue);
Dialog.addNumber("Green Threshold Level:", green);
Dialog.addDirectory("Directory:", File.getDefaultDir);
Dialog.show();
red = Dialog.getNumber();
blue = Dialog.getNumber();
green = Dialog.getNumber();
path = Dialog.getString();

//run processing function on every czi or tif in the directory, save results in one .csv
files = getFileList(path);
tabletitle = "Total Results";
Table.create(tabletitle);
Table.set("Image Name", 0, "Chosen Threshold Levels");
Table.set("Value Subtracted From Red", 0, "N/A", tabletitle);
Table.set("Red Otsu Threshold", 0, "Top "+red+"% of pixels");
Table.set("Value Subtracted From Blue", 0, "N/A", tabletitle);
Table.set("Blue Otsu Threshold", 0, "Top "+blue+"% of pixels");
Table.set("Value Subtracted From Green", 0, "N/A", tabletitle);
Table.set("Green Otsu Threshold", 0, "Top "+green+"% of pixels");
row = 1;
for (i=0; i<files.length; i++) {
	if (endsWith(files[i], ".czi") || endsWith(files[i], ".tif")) {
		imagenoext = File.getNameWithoutExtension(path+files[i]);
		result = Processor(path+files[i], path, files[i], imagenoext, red, blue, green);
		Table.set("Image Name", row, imagenoext, tabletitle);
		Table.set("Value Subtracted From Red", row, result[0], tabletitle);
		Table.set("Red Otsu Threshold", row, result[1], tabletitle);
		Table.set("Value Subtracted From Blue", row, result[2], tabletitle);
		Table.set("Blue Otsu Threshold", row, result[3], tabletitle);
		Table.set("Value Subtracted From Green", row, result[4], tabletitle);
		Table.set("Green Otsu Threshold", row, result[5], tabletitle);
		row = row + 1;
	}
}

Table.save(path+"Results.csv", tabletitle);

winlist = getList("window.titles");
for (i=0; i<winlist.length; i++){
     winame = winlist[i];
     selectWindow(winame);
     run("Close");
}
close("*");
print("Finished!");


function Processor(fullpath, path, image, imagenoext, red, blue, green) {
// Subtracts background from images and determines Otsu threshold

channels = newArray("red", "blue", "green");
channelsets = newArray("C1-", "C2-", "C3-");
channelthresh = newArray(red, blue, green);
result = newArray(0);
for (k = 0; k < channels.length; k++) {

//Open image
run("Bio-Formats Windowless Importer", "open=[fullpath]");
imagename = "Image.tif";
rename(imagename);
run("Split Channels");
selectWindow(channelsets[k]+imagename);
close("\\Others");

//Set desired threshold value
bottom = (100 - channelthresh[k])/100;

//Get Counts for Each Pixel Value in each slice
totalcounts = newArray(65536);
for (i = 0; i < nSlices; i++) {
	setSlice(i+1);
	getHistogram(values, slicecounts, 65536);
	for (j = 0; j < totalcounts.length; j++) {
		newcount = totalcounts[j] + slicecounts[j];
		totalcounts[j] = newcount;
	}
}

//Find out what pixel value (threshold) constitutes bottom whatever %
height = getHeight();
width = getWidth();
numpixels = width*height*nSlices;

sum = 0;
threshold = 0;
while (sum < numpixels*bottom) {
	sum = sum + totalcounts[threshold];
	threshold = threshold + 1;
  }

//set all counts below threshold to 0
thresholdedcounts = Array.copy(totalcounts);
for (i = 0; i < threshold; i++) {
	thresholdedcounts[i] = 0;
}

//Get average of all pixel values in top whatever percent

value = 0;
totalsum = 0;
for (i = 0; i < thresholdedcounts.length; i++) {
	sum = thresholdedcounts[i]*value;
	totalsum = totalsum + sum;
	value = value + 1;
}

mean = totalsum / (numpixels*(1-bottom));
background = mean*0.5;

//Subtract background (half the average of the top whatever % of pixels) from image and save
run("Subtract...", "value="+background+" stack");
saveAs("Tiff", path+imagenoext+"_"+channels[k]+"_background_subtracted.tif");

//Determine Otsu Threshold Value
setAutoThreshold("Otsu dark stack");
getThreshold(lower, upper);

//Save values in results array
myarray = newArray(background, lower);
result = Array.concat(result,myarray);

close("*");

}

//make combined channel images
open(path+imagenoext+"_blue_background_subtracted.tif");
open(path+imagenoext+"_green_background_subtracted.tif");
open(path+imagenoext+"_red_background_subtracted.tif");
run("Merge Channels...", "c1=["+imagenoext+"_red_background_subtracted.tif] c2=["+imagenoext+"_green_background_subtracted.tif] c3=["+imagenoext+"_blue_background_subtracted.tif] create ignore");
saveAs("Tiff", path+imagenoext+"_merged.tif");

return result;

}
