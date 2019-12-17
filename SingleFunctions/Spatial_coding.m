%Author: Tori-Lynn Temple, started Feburary 13th, 2018
%Spatial Coding for JB and Guillaume 

%%BINNING EACH OF THE CELLS IN A VIDEO  
%Here we pass in ms, behav to bin each cell

%here I read a video and only take the first frame for reference : 
videoObj = VideoReader ( from behav to get field);
first_frame = read(videoObj,1);

BinningTheCells(ms, behav, first_frame);

%%MEAN FIRING RATE INSIDE THE BIN


%%MEAN FIRING RATE OUTSIDE THE BIN 



%%PREDICTING THE MOUSES PLACEMENT BASED ON STATISTICS 


%%BAYESIAN DECODING 
