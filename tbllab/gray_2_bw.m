close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 8; % 받아들인 사진의 개수 (계속 바꿔줘야함!!!)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index1 = alpha : beta;
    
    number_pixel = 0; % 픽셀 개수를 세기 위한 초기값 지정
    
    FILEnum = num2str(index1);
    FILEname = strcat('D:\TLB Lab (경인 - MATLAB)\경인 Sample (Volume Calculation)\1 (',FILEnum,').bmp');
    imdata_raw = imread(FILEname); % read image
    [m,n] = size(imdata_raw); % 가로,세로 픽셀 수를 받아들인다
    
    imdata_double = im2double (imdata_raw);
    
    for ii = 1:m
        for jj = 1:n
            if imdata_double(ii,jj) < 1;
                imdata_double(ii,jj) = 0;
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [l,n]=bwlabel(imdata_double);

    max_points=find(l==1);

    for index2 = 2:n
        
        points = find(l==index2);
        
        if(length(points) > length(max_points))
            max_points = points;
            
        end
        
    end
    
    
    [r,c]=size(imdata_double);
    label_result=zeros(r,c);
    
    for index3 = 1:length(max_points)
        
        label_result(max_points(index3)) = 1;
        
    end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = strcat('result',num2str(index1),'.bmp');
imwrite(label_result,name);

close all;

end
