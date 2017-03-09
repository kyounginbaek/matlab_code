close all; clear all; clc;

alpha = 1; % �̹��� ���� �ѹ�
beta = 8; % �޾Ƶ��� ������ ���� (��� �ٲ������!!!)

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
