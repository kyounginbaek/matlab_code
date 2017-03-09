close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 8; % 받아들인 사진의 개수 (계속 바꿔줘야함!!!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\TBL LAB\Sample - label max\1 (',FILEnum,').bmp');
    imdata_raw = imread(FILEname); % read image
  
    boundary = edge(imdata_raw, 'canny');
    
    name = strcat('result',num2str(index),'.bmp');
    imwrite(boundary,name);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
