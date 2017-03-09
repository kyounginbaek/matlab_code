close all; clear all; clc;

% code optimization by vectorizing loops
% 
% M = 2012;
% N = 2012;
% u0 = 1/(4*pi);
% v0 = 1/(4*pi);
% A = 1;
% g1 = zeros(M,N);
% tic
% g2 = zeros(M,N);
% toc
% 
% tic
% r = 0:M-1;
% c = 0:N-1;
% [C,R] = meshgrid(c,r);
% g1 = A*sin(u0*R+v0*C);
% 
% t1 = toc
% 
% tic
% for i = 1:M;
%     for j = 1:N;
%         g2(i,j) = A*sin(u0*(i-1) + v0*(j-1));
%     end
% end
% t2 = toc

% img = imread('3.jpg');
img = imread('coins.png');
% % figure('Name','Original image'); imshow(img);
% showHist(img);
% 
% % basic image adjustments
% imga = imadjust(img,[0 0 0; 1 1 1],[],0.5);
% % figure('Name','Adusted image'); imshow(ima);
% showHist(imga);
% 
% % image negative
% imgNeg = imcomplement(img);
% showHist(imgNeg);

% img = rgb2gray(img);
% imh  = imhist(img);
% imch = cumsum(imh);
% figure(); 
%     hold on;
%     plot(imh);
%     plot(imch,'r');
%     
%     
% imgeq = histeq(img);
% figure(); imshow(imgeq);
% imh  = imhist(imgeq);
% imch = cumsum(imh);
% figure(); 
%     hold on;
%     plot(imh);
%     plot(imch,'r');

    
% thresholding
imgBin = im2bw(img,0.5);
figure();imshow(imgBin);
figure('Name','Binarization examples');
    subplot(2,2,1); imshow(img); title('Original');
    subplot(2,2,2); imshow(im2bw(img,0.3)); title('T=0.3');
    subplot(2,2,3); imshow(im2bw(img,0.6)); title('T=0.6');
    subplot(2,2,4); imshow(im2bw(img,0.7)); title('T=0.8');
    
% Otsu
level = graythresh(img);
imgBin = im2bw(img,level);
figure('Name','Binarization with Otsu threshold');
imshow(imgBin);

close all;
T = 0;
TNew = 255/2;
while abs(TNew - T) > 1
    indLow = find(img < TNew);
    indHi  = find(img >= TNew);
    
    m1 = mean(img(indLow));
    m2 = mean(img(indHi));    
    T = TNew;
    TNew = (m1+m2)/2;
end;

figure('Name','Histogram-based thresholding');
imshow(im2bw(img,T/255));


% balanced histogram method
H = imhist(img);
Il = 1;
Ir = 256;
Im = round((Ir+Il)/2);
Wl = sum(H(Il:Im));
Wr = sum(H(Im+1:Ir));
while Il ~= Ir
    if Wr > Wl
        Ir = Ir - 1;        
    else 
        Il = Il + 1;  
    end
    Im = round((Ir+Il)/2);
    Wl = sum(H(Il:Im));
    Wr = sum(H(Im+1:Ir));
end

figure('Name','Balanced histogram thresholding');
imshow(im2bw(img,Im/255));