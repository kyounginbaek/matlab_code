close all; clear all; clc;

alpha = 1; % starting input image number
beta = 1000; % ending input image number

for index = alpha:beta
    
    
    
% image processing first box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FILEnum = num2str(index);
    FILEname = strcat('D:\800nm #29 80mJ 0hr\M\test_01_',FILEnum,'.bmp');
    imdata_raw = imread(FILEname); % read image
    [m_i,n_i,l_i] = size(imdata_raw);
    
    imdata_gray = rgb2gray(imdata_raw);
    imdata_double = im2double(imdata_gray); % 이거 왜 해야되는지 잘 모르겠다.
    imdata_ROI = imdata_double(10:m_i - 70,40:n_i-30);
    [m,n,l] = size(imdata_ROI)
    % figure('Name','imdata ROI'); imshow(imdata_ROI);
    
    imdata_filter = medfilt2(imdata_ROI);
    % figure('Name','imdata filter'); imshow(imdata_filter);
    
    
    % 2차 이미지 프로세싱: 밝기 구간 분할 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imdata_intensity_1 = imdata_filter;
    imdata_intensity_2 = imdata_filter;
    imdata_intensity_3 = imdata_filter;
    
    for ii = 1:m % 2/3 이미지 데이터 군
        for jj = 1:n
            if imdata_filter(ii,jj) < 1/2
                imdata_intensity_1(ii,jj) = 0;
            end
        end
    end
    % figure('Name','imdata intensity 1'); imshow(imdata_intensity_1);
    
    for ii = 1:m % 1/8 이미지 데이터 군
        for jj = 1:n
            if imdata_filter(ii,jj) < 1/8
                imdata_intensity_2(ii,jj) = 0;
            end
        end
    end
    % figure('Name','imdata intensity 2'); imshow(imdata_intensity_2);
    
    for ii = 1:m % 7/8 이미지 데이터 군
        for jj = 1:n
            if imdata_filter(ii,jj) > 1/8
                imdata_intensity_3(ii,jj) = 0;
            end
        end
    end
    % figure('Name','imdata intensity 3'); imshow(imdata_intensity_3);
    
    
    % 3차 이미지 프로세싱: 1차 노이즈 제거 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filter_1 = strel('disk',4);
    imdata_intensity_1 = imdilate(imdata_intensity_1,filter_1);
    % figure('Name','imdata intensity 1'); imshow(imdata_intensity_1);
    
    [p,q] = bwlabel(imdata_intensity_1);
    
    for a = 1:q
        points = find(p==a);
        if(length(points) < 500)
            for jndex = 1:length(points)
                imdata_intensity_1(points(jndex))=0;
            end
        end
    end
    % figure('Name','imdata intensity 1'); imshow(imdata_intensity_1);
    
    for jj = 1:n % 가장 밝은 데이터를 기준으로 노이즈 제거
        point = find(imdata_intensity_1(:,jj) ~= 0);
        if length(point) ~= 0;
            point = min(point);
            imdata_intensity_2(1:point,jj) = 0;
            imdata_intensity_3(1:point,jj) = 0;
        end
    end
    % figure('Name','imdata intensity 2'); imshow(imdata_intensity_2);
    % figure('Name','imdata intensity 3'); imshow(imdata_intensity_3);
    
    
    % 4차 이미지 프로세싱: 2차 노이즈 제거 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filter_2 = strel('disk',3);
    imdata_intensity_2 = imdilate(imdata_intensity_2,filter_2);
    % figure('Name','imdata intensity 2'); imshow(imdata_intensity_2);
    
    [p,q] = bwlabel(imdata_intensity_2);
    
    for a = 1:q
        points = find(p==a);
        if(length(points) < 300)
            for jndex = 1:length(points)
                imdata_intensity_2(points(jndex))=0;
            end
        end
    end
    % figure('Name','imdata intensity 2'); imshow(imdata_intensity_2);
    
    for jj = 1:n % 가장 밝은 데이터를 기준으로 노이즈 제거
        point = find(imdata_intensity_2(:,jj) ~= 0);
        if length(point) ~= 0;
            point = min(point);
            imdata_intensity_3(1:point,jj) = 0;
        end
    end
    % figure('Name','imdata intensity 3'); imshow(imdata_intensity_3);
    
    
    % 5차 이미지 프로세싱: OTSU'S METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    threshold_2 = graythresh(imdata_intensity_2);
    threshold_3 = graythresh(imdata_intensity_3);
    
    imdata_otsu_2 = im2bw(imdata_intensity_2,threshold_2);
    imdata_otsu_3 = im2bw(imdata_intensity_3,threshold_3);
    
    imdata_otsu = imdata_otsu_2 + imdata_otsu_3;
    % figure(); imshow(imdata_otsu);
    
    
    % 6차 이미지 프로세싱: 마스크로 노이즈 제거 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filter_3 = strel('disk',3);
    imdata_mask = imdilate(imdata_otsu,filter_3);
    % figure('Name','imdata otsu'); imshow(imdata_otsu);
    
    for ii = 1:m % 원본데이터를 마스크 씌우기
        for jj = 1:n
            if imdata_mask(ii,jj) == 0
                imdata_ROI(ii,jj) = 0;
            end
        end
    end
%     figure('Name','imdata ROI noise reduced'); imshow(imdata_ROI);
    
    
    
% image processing second box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1차 이미지 프로세싱: OTSU'S METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold = graythresh(imdata_ROI);
    imdata_otsu = im2bw(imdata_ROI,threshold);
%     figure('Name','imdata otsu'); imshow(imdata_otsu);
    
    % 2차 이미지 프로세싱: 거리 구하기 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for jj = 1:n
        point_tmp = find(imdata_otsu(:,jj) == 1);
        if length(point_tmp) ~= 0
           point(jj) = min(point_tmp);
        else
            if jj ~= 1
                point(jj) = point(jj - 1);
            end
        end
    end
    
    x_loc = 1:n;
    y_loc = point;
    
    % 3차 이미지 프로세싱: 트랜드 곡선 구하기 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xx = 1:n;
    
    constant = 4;
    remain = 1:constant - 1;
    
    for kk = 1:floor(n/constant) % 4개 기준으로 샘플링
        first_x_point(kk) = x_loc(constant*(kk-1) + remain(1));
        first_y_point(kk) = y_loc(constant*(kk-1) + remain(1));
    end
    first_y = first_y_point;
    first_x = first_x_point;
    
    yy1 = spline(first_x,first_y,xx);
    
    for kk = 1:floor(n/constant) % 4개 기준으로 샘플링
        second_x_point(kk) = x_loc(constant*(kk-1) + remain(2));
        second_y_point(kk) = y_loc(constant*(kk-1) + remain(2));
    end
    second_y = second_y_point;
    second_x = second_x_point;
    
    yy2 = spline(second_x,second_y,xx);
    
    for kk = 1:floor(n/constant) % 4개 기준으로 샘플링
        third_x_point(kk) = x_loc(constant*(kk-1) + remain(3));
        third_y_point(kk) = y_loc(constant*(kk-1) + remain(3));
    end
    third_y = third_y_point;
    third_x = third_x_point;
    
    yy3 = spline(third_x,third_y,xx);
    
    for kk = 1:floor(n/constant) % 4개 기준으로 샘플링
        fourth_x_point(kk) = x_loc(constant*kk);
        fourth_y_point(kk) = y_loc(constant*kk);
    end
    fourth_y = fourth_y_point;
    fourth_x = fourth_x_point;
    
    yy4 = spline(fourth_x,fourth_y,xx);
    
%     figure();
%     plot(xx,yy1,'y'); hold on;
%     plot(xx,yy2,'k');
%     plot(xx,yy3,'b');
%     plot(xx,yy4,'m');
    
    for jj = 1:n
        sampled_y(jj) = round((yy1(jj) + yy2(jj) + yy3(jj) + yy4(jj))/4);
        if sampled_y(jj) <= 0
            sampled_y(jj) = 1;
        end
    end
    
    for jj = 1:n
        imdata_ROI(sampled_y(jj),jj) = 1;
    end
    figure(); imshow(imdata_ROI);
    
% image processing third box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1차 이미지 프로세싱: 이미지 저장 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name = strcat('result',num2str(index),'.bmp');
    imwrite(imdata_ROI,name); 
    
    close all;
        
    
end
    
