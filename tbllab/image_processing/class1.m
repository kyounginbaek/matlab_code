clear all; close all; clc;

fileName = 'rgb2gray4.jpg';
% fileName = 'img1.jpg';

%% basic image manipulations
% getting image info
imgInfo = imfinfo(fileName)

% reading image
img = imread(fileName);

% getting image size
[H,W] = size (img);

% information about array
whos img

% dispaying image
% imshow(img);

% saving to disk
imwrite(img,'newfile.jpg','quality',100);


%% convert to grayscale
img = imread(fileName);
figure('Name','Original image');imshow(img);

imgGray = rgb2gray(img);
figure('Name','Grayscale image (rgb2gray)');imshow(imgGray);

imgGray2 = img(:,:,1)/3 + img(:,:,2)/3 + img(:,:,3)/3;
figure('Name','Grayscale image (R/3 + G/3 + B/3');imshow(imgGray2);


%% intensity histogram
[imgHist,x] = imhist(rgb2gray(img));
figure('Name','Image histogram'); bar(imgHist);


%% some image manipulation techniques;
% brightness change
figure('Name','Brightness manipulation');
hold on;
subplot(1,3,1); imshow(img,[1 255]); title('Original');
subplot(1,3,2); imshow(img+70,[1 255]); title('Brighter');
subplot(1,3,3); imshow(img-70,[1 255]); title('Darker');

% contrast change
figure('Name','Contrast manipulation');
hold on;
subplot(1,3,1); imshow(img,[1 255]); title('Original');
subplot(1,3,2); imshow(img*1.5,[1 255]); title('Increased contrast');
subplot(1,3,3); imshow(img*0.5,[1 255]); title('Decreased contrast');


%% single channel manipulation;
% getting RGB channels
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

% ... and histograms
hR = imhist(R);
hG = imhist(G);
hB = imhist(B);

figure('Name','Image histograms');
    hold on;
    area(imgHist,'FaceColor',[0.6 0.6 0.6]);
    plot(hR,'r');
    plot(hG,'g');
    plot(hB,'b');

% changing brigtness of single channel
figure('Name','Green channel brightness manipulation');
    hold on;
    subplot(1,3,1); imshow(img,[1 255]); title('Original');
    subplot(1,3,2); imshow(cat(3,R,G-100,B)); title('Decreased brightness');
    subplot(1,3,3); imshow(cat(3,R,G+100,B)); title('Decreased brightness');


%% color space changing
HSVimg = rgb2hsv(img);
% getting Hue, Saturation and Intensity channels
H = HSVimg(:,:,1);
S = HSVimg(:,:,2);
V = HSVimg(:,:,3);

% HSV channels
figure('Name','HSV channels');
    hold on;
    subplot(2,2,1); imshow(img,[1 255]); title('Original image');
    subplot(2,2,2); imshow(H,[]); colormap(hsv); title('Hue');
    subplot(2,2,3); imshow(cat(3,S,S,S),[]); title('Saturation');
    subplot(2,2,4); imshow(cat(3,V,V,V),[]); title('Value');


%% image binarization
% binImg = im2bw(0.5);
figure('Name','Image thresholding');
    hold on;
    subplot(2,2,1); imshow(img,[1 255]); title('Original image');
    subplot(2,2,2); imshow(H,[]); colormap(hsv); title('Hue');
    subplot(2,2,3); imshow(cat(3,S,S,S),[]); title('Saturation');
    subplot(2,2,4); imshow(cat(3,V,V,V),[]); title('Value');


%% working with figures
% displayin multiple images in one figure
out = figure('Name',fileName);
    hold on;
    sp = subplot(2,6,[1 2, 7 8]); subimage(img); title('Original image');
    subplot(2,6,3); imshow(img(:,:,1)); title('Red');
    subplot(2,6,4); imshow(img(:,:,2)); title('Green');
    subplot(2,6,9); imshow(img(:,:,3)); title('Blue');
    subplot(2,6,10); bar(imgHist); title('Histogram');
    subplot(2,6,[5 6, 11 12]); subimage(rgb2gray(img)); title('Grayscale');

%getting screen parameters
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
figPos = [ (scnsize(3)-1200) / 2, (scnsize(4)-350) / 2, 1200, 330];
set(out,'Position',figPos);