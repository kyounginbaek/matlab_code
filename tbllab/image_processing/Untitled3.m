close all; clear all; clc;

img = imread('coins.png');
figure(); imshow(img);
H = imhist(img);

level = graythresh(img);
imgBin = im2bw(img, level);
figure('Name', 'otsu');
imshow(imgBin);