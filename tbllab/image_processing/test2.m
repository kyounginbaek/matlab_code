close all; clear all; clc;
img = imread('coins.png');
figure(); imshow(img);
bw = im2bw(img,0.5);
figure(); imshow(bw);

bw2 = imfill(bw, 'holes');
figure(); imshow(bw2);

bw3 = bwmorph(bw2, 'remove');
figure(); imshow(bw3);

[L,num] = bwlabel(bw2);
props = regionprops(L);

figure();
imshow(img);

hold on;
cmap = hsv(num);
cmap2 = summer(num);
for i=1:num
    [r,c] = find(L == i);
    plot(c,r,'.', 'Color', cmap(i,:));
    box = round(props(i).BoundingBox);
    rectangle('Position',box, 'EdgeColor', cmap2(i,:), 'LineWidth',2);
    
    center = round(props(i).Centroid);
    plot(center(1), center(2),'xb','MarkerSize', 15);
end