clc; clear; close all;

imdata_tmp = imread('C:\Users\xnote\Desktop\result_1.jpg');
figure(1); imshow(imdata_tmp);
figure(2); imhist(imdata_tmp); % 이미지를 gray scale로 바꿔준다

[m,n] = size(imdata_tmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_double = im2double(imdata_tmp);
figure(); imhist(imdata_double);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filter = 1/25; % 5X5 average filter의 mask를 생성
mask = [filter filter filter filter filter; 
    filter filter filter filter filter; 
    filter filter filter filter filter;
    filter filter filter filter filter;
    filter filter filter filter filter];

d = zeros(m+4,n+4); % 5X5 average filter의 원리
d(1:m,1:n) = imdata_double;

imdata_filter = zeros(m,n);

for row = 1:m
    for clm = 1:n
        tmp = d(row:row+4,clm:clm+4).*mask; % 각 행열 간의 곱
        tmp2 = sum(tmp); % 각 행 간의 합
        tmp3 = sum(tmp2); % 각 열 간의 합
        imdata_filter(row,clm) = tmp3;
    end
end

figure(); imhist(imdata_filter);
figure(); imshow(imdata_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_threshold = zeros(m,n);
threshold_line = 0.3; % thereshold_line 기준을 설정한다

for row = 1:m
    for clm = 1:n
        if imdata_filter(row, clm) < threshold_line
            imdata_threshold(row, clm) = 0;
        else 
            imdata_threshold(row, clm) = 1;            
        end
    end
end

figure(); imhist(imdata_threshold); 
figure(); imshow(imdata_threshold);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%