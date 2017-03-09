clc; clear; close all;

imdata_tmp = imread('C:\Users\xnote\Desktop\result_1.jpg');
figure(); imshow(imdata_tmp);

threshold = graythresh(imdata_tmp);
imdata_tmp = im2bw(imdata_tmp,threshold);
figure(); imshow(imdata_tmp);


h = fspecial('average',5);
imfilt_tmp = imfilter(imdata_tmp,h);
figure(); imshow(imfilt_tmp);
