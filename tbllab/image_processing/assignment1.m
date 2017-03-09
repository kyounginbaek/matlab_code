close all; clear all; clc;

original = imread('img1.jpg');

%% assignment1

figure('Name', 'assignment1');

R = original(:,:,1);
G = original(:,:,2);
B = original(:,:,3);

HR = imhist(R);
HG = imhist(G);
HB = imhist(B);

subplot(3,2,1); imshow(R); title('red image');
subplot(3,2,2); plot(HR,'r'); title('red histogram');
subplot(3,2,3); imshow(G); title('green image');
subplot(3,2,4); plot(HG,'g'); title('green histogram');
subplot(3,2,5); imshow(B); title('blue image');
subplot(3,2,6); plot(HB,'b'); title('blue histogram');

%% assignment2

figure('Name', 'assignment2');

gray = rgb2gray(original);

subplot(3,5,[1 2, 6 7, 11 12]); imshow(original); title('original image');
subplot(3,5,3); imshow(R); title('red image');
subplot(3,5,8); imshow(G); title('green image');
subplot(3,5,13); imshow(B); title('blue image');
subplot(3,5,[4 5, 9 10, 14 15]); imshow(gray); title('gray image');

%% assignment3

figure('Name', 'assignment3');

subplot(3,6,[1 2 3, 7 8 9, 13 14 15]); imshow(original);
subplot(3,6,[4 5 6]); plot(HR, 'r');
subplot(3,6,[10 11 12]); plot(HG, 'g');
subplot(3,6,[16 17 18]); plot(HB, 'b');

% assignment4

figure('Name', 'assignment4');

original2 = imread('rgb2gray4.jpg');

HSV = rgb2hsv(original2);
H = HSV(:,:,1);

hold on 
subplot(1,2,1); imshow(original2); title('original image');
subplot(1,2,2); imshow(H); title('hue image');










