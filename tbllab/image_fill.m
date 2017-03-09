close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 36; % 받아들인 사진의 개수

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : alpha%beta
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\hole3\1 (',FILEnum,').jpg');
    imdata_raw = imread(FILEname);
    
    im_check = imdata_raw(61:793, 541:659,1);
    im_bw = im2bw(im_check, graythresh(im_check));
    figure(1); imshow(im_bw);
    [hole_bound, im_hole] = bwboundaries(im_bw, 'holes');
    hole_boundary = hole_bound{1};
%     hole_boundary(:, 1)
%     hole_boundary(:, 2)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     name = strcat('result',num2str(index),'.bmp');
%     imwrite(im_fill,name);

    close all;

end