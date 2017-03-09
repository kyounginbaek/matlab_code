close all; clear all; clc;
img = imread('rgb2gray2.jpg');
figure(); imshow(img);
figure(); imshow(rgb2gray(img));
imgHSV = rgb2hsv(img);
S = imgHSV(:,:,2);
figure(); imshow(S);

bw = im2bw(S,0.3);
figure(); imshow(bw);

bwCl = bwmorph(bw,'close');
figure(); imshow(bwCl);

se = strel('rectangle', [4 4]);
bwCl2 = imclose(bw,se);
figure(); imshow(bwCl2);
bwOpened = bwareaopen(bwCl2, 40);
figure(); imshow(bwOpened);

[L,num] = bwlabel(bwOpened);
figure(); imshow(L);
num % number of objects

figure(); imshow(L,[]);
colormap('summer');
colormap('jet');

imshow(img);
hold on
[r,c] = find(L == 1);
plot(c,r,'.g');
[r,c] = find(L == 2);
plot(c,r,'.b');

prop = regionprops(L);
prop(1)
prop(2)

box = round(prop(1).BoundingBox);
box

rectangle('Position', box, 'EdgeColor', [1 0 0]);

box2 = round(prop(2).BoundingBox);

rectangle('Position', box2, 'EdgeColor', [0 0 1]);
