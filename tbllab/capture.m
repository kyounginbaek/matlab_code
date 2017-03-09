close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 8; % 받아들인 사진의 개수

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta;
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\3D structure picture\andy\3rd hole - 복사본\1 (',FILEnum,').jpg');
    imdata_raw = imread(FILEname); % read image
    
    im_capture = imdata_raw(1:606, 481:950, :); % 좌표 추출 (height, width, channel)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = strcat('result',num2str(index),'.jpg');
    imwrite(im_capture,name);

    close all;

end
