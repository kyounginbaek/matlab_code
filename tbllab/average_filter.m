clc; clear; close all;

imdata_tmp = imread('C:\Users\xnote\Desktop\result_1.jpg');
figure(1); imshow(imdata_tmp);
figure(2); imhist(imdata_tmp); % �̹����� gray scale�� �ٲ��ش�

[m,n] = size(imdata_tmp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_double = im2double(imdata_tmp);
figure(); imhist(imdata_double);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filter = 1/25; % 5X5 average filter�� mask�� ����
mask = [filter filter filter filter filter; 
    filter filter filter filter filter; 
    filter filter filter filter filter;
    filter filter filter filter filter;
    filter filter filter filter filter];

d = zeros(m+4,n+4); % 5X5 average filter�� ����
d(1:m,1:n) = imdata_double;

imdata_filter = zeros(m,n);

for row = 1:m
    for clm = 1:n
        tmp = d(row:row+4,clm:clm+4).*mask; % �� �࿭ ���� ��
        tmp2 = sum(tmp); % �� �� ���� ��
        tmp3 = sum(tmp2); % �� �� ���� ��
        imdata_filter(row,clm) = tmp3;
    end
end

figure(); imhist(imdata_filter);
figure(); imshow(imdata_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_threshold = zeros(m,n);
threshold_line = 0.3; % thereshold_line ������ �����Ѵ�

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