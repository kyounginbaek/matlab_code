close all; clear all; clc;

alpha = 1; % �̹��� ���� �ѹ�
beta = 8; % �޾Ƶ��� ������ ����

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta;
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\3D structure picture\andy\3rd hole - ���纻\1 (',FILEnum,').jpg');
    imdata_raw = imread(FILEname); % read image
    
    im_capture = imdata_raw(1:606, 481:950, :); % ��ǥ ���� (height, width, channel)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = strcat('result',num2str(index),'.jpg');
    imwrite(im_capture,name);

    close all;

end
