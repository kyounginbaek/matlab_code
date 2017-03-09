close all; clear all; clc;
% keywords: bwmorpth, imclose, imopen, imerode, indilate, imfill, strel,
% bwareaopen

% loading image
img = imread('rgb2gray2.jpg');

% converting to HSV
imHSV = rgb2hsv(img);

% use saturation channel for thresholding
S = imHSV(:,:,2);
% figure('Name','Saturation'); imshow(S);
bw = im2bw(S,0.3);
figure('Name','Binarized image'); imshow(bw);

% apply 'closing' operation
bwClosed = bwmorph(bw,'close');
figure('Name','Morphology "Closing"'); imshow(bwClosed);

% another way to apply 'closing' with different structure element
SE = strel('rectangle',[5 5]);
bwCl2 = imclose(bw,SE);
figure('Name','Morphology "Bridge" + "Closing"'); imshow(bwCl2);

% 'bridge' followed by 'closing'
bwCl3 = bwmorph(bwmorph(bw,'bridge'),'close');
figure('Name','Morphology "Bridge" + "Closing"'); imshow(bwCl3);

bwFil = bwareaopen(bwCl2,60);
figure('Name','Removed small regions'); imshow( bwFil);

[L,num] = bwlabel(bwFil);
figure('Name','Labeled image'); imshow(L,[]); colormap('summer');

% displaying 
figure(); imshow(img);
% computing blob properties
prop = regionprops(L,'BoundingBox');
hold on;
cmap = hsv(num); %preparing color map;
for i=1:num
    [r,c] = find(L == i);
    plot(c,r,'.','Color',cmap(i,:));
    rectangle('Position',prop(i).BoundingBox,'EdgeColor','b','LineWidth',2);
end


% filling example
coins   = im2bw(imread('coins.png'));
figure('Name','Coins image'); imshow(coins);

coinsBW = imfill(coins,'holes');
figure('Name','Filled image'); imshow(coinsBW);

coinsBound = bwmorph(coinsBW,'remove');
figure('Name','Filled image'); imshow(coinsBound);


[L,num] = bwlabel(coinsBW);

% displaying 
figure(); imshow(coins);
% computing blob properties
prop = regionprops(L,'BoundingBox');
hold on;
cmap = hsv(num); %preparing color map;
for i=1:num
    [r,c] = find(L == i);
    plot(c,r,'.','Color',cmap(i,:));
    rectangle('Position',prop(i).BoundingBox,'EdgeColor','b','LineWidth',2);
end

imskel = bwmorph(coinsBW,'skel',Inf);
[r,c] = find(imskel > 0);
plot(c,r,'.w');