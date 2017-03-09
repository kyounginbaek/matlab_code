clc; clear all; close all;

image_tmp = imread('C:\Users\xnote\Desktop\result_1.jpg');
figure(1); imshow(image_tmp);

mask = fspecial('average', 3); % mask를 생성
image_filter = imfilter(image_tmp,mask);
figure(2); imshow(image_filter);

image_gray = graythresh(image_filter); % gray 스케일을 bw로 바꾼다
image_binary = im2bw(image_filter,image_gray);
figure(3); imshow(image_binary);

f = strel('disk', 4);
image_dilate = imdilate(image_binary,f);
figure(4); imshow(image_dilate);
% g = imdilate(e,f);

image_erode = imerode(image_dilate,f);
figure(5); imshow(image_erode);

image_fill = imfill(image_erode, 'holes');
figure(6); imshow(image_fill);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[l,n]=bwlabel(image_fill);

max_points=find(l==1);

for index2 = 2:n
    
    points = find(l==index2);
    
    if(length(points) > length(max_points))
        max_points = points;
        
    end
    
end


[r,c]=size(image_fill);
label_result=zeros(r,c);

for index3 = 1:length(max_points)
    
    label_result(max_points(index3)) = 1;
    
end
    
figure(7); imshow(label_result);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
[m,n] = size(label_result);

image_result = zeros(m,n);

image_double = im2double(image_filter);
figure(8); imhist(image_filter);
figure(9); imhist(image_double);

for row = 1:m
    for clm = 1:n
        image_result(row,clm) = image_double(row,clm) - ~label_result(row,clm);
        if image_result(row,clm) < 0
            image_result(row,clm) = 0;
        end
        
    end
end

figure(10); imshow(image_result);

% I = imopen(k,f);
% I = imclose(k,f);
% dilate와 erode 과정을 따로 안 쓰고 생략이 가능