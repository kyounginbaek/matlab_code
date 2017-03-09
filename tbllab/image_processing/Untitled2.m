close all; clear all; clc;

img = imread('3.jpg');
figure(); imshow(rgb2gray(img));
H = imhist(rgb2gray(img));
h1 = cumsum(H);
figure(); area(H);

imgAdj = histeq(rgb2gray(img)); % more contrast
figure(); imshow(imgAdj);
H = imhist(imgAdj);
h2 = cumsum(H);
figure(); area(H);

figure();
    hold on;
    plot(h1, 'r');
    plot(h2, 'g'); % more linear