function Thresholds(imagefilepath) { 
//gets maximal Otsu threshold for each channel of an image
//open, split, and rename channels
	run("Bio-Formats Windowless Importer", "open=[imagefilepath]");
	run("Split Channels");
	titlelist = getList("image.titles");
	
//obtain threshold information from each channel
	resultsarray = newArray(titlelist.length);
	for (i = 0; i < titlelist.length; i++) {
		selectWindow(titlelist[i]);
		thresholdarray = newArray(nSlices);
		for (j = 1; j <= nSlices; j++) {
			setSlice(j);
			setAutoThreshold("MaxEntropy dark");
			getThreshold(lower, upper);
			thresholdarray[j-1] = lower;
		}
		Array.getStatistics(thresholdarray, min, max, mean, stdDev);
		resultsarray[i] = max;
	}

//end and return results
	close("*");
	Table.setColumn(image, resultsarray);
	
}



//Run Thresholds function on all czis in the directory
path = getDirectory("Choose a Directory ");
files = getFileList(path);
Table.create("Thresholds");

for (i=0; i<files.length; i++) {
	image = files[i];
	if (endsWith(image, ".czi")) {
	Thresholds(path+image);
	}
}

//Save table, close all windows, and end
Table.save(path + "Results.csv");
list = getList("window.titles");
    	for (i=0; i<list.length; i++){
    		winame = list[i];
			selectWindow(winame);
			run("Close");
     	}

print("Finished!");