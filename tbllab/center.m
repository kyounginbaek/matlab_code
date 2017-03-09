close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 8; % 받아들인 사진의 개수(계속 바꿔줘야함)

divide = 36; % 360도를 얼마큼 나눌것인지

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\TBL LAB\Sample - label max\1 (',FILEnum,').bmp'); % imfill 되어있는 data를 가져온다
    imdata_raw = imread(FILEname); 
    
    % 세로(행),가로(열) 픽셀 수를 받아들인다
    [row, column] = size(imdata_raw); 
    
    % double 형태로 바꿔준다.
    imdata_double = im2double(imdata_raw);
    
    % 중심점을 찾는다
    [L, n] = bwlabel(imdata_double);
    [r, c] = find(L == 3);
    rbar = mean(r); % 세로(행)의 중심
    cbar = mean(c); % 가로(열)의 중심
    
    % 중심점을 Marker를 이용하여 보여준다
    figure(); imshow(imdata_double);
    hold on
    for k = 1:n
        [r, c] = find(L == k);
        rbar = mean(r);
        cbar = mean(c);
        plot(cbar, rbar, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFacecolor', 'k', 'MarkerSize', 10);
        plot(cbar, rbar, 'Marker', '*', 'MarkerEdgeColor', 'w');
    end
        
    % 중심점을 기준으로 선을 그린다
    imdata_boundary = edge(imdata_raw, 'canny'); % boundary 이미지를 생성한다
    imdata_line = double(imdata_boundary);
    
    % 중심좌표 정수값 변환 (해야 되나?)
    round_rbar = round(mean(r));
    round_cbar = round(mean(c));
    
    figure(); imshow(imdata_double);
    hold on
    for multiple = 1:divide 
        for x = 1:column
            y = tan((2*pi/divide)*multiple) * (x - cbar) + rbar ; 
            plot(x,y,'w'); 
        end
    end

    for x = 1:column % x의 좌표
        for multiple = 1:divide
            y = round ( tan((2*pi/divide)*multiple) * (x - cbar) + rbar ); % y의 좌표
            imdata_line(y,x) = 1; 
        end
    end
    
    figure(); imshow(imdata_line);
    
%     for y = 1:row
%         x = round_cbar;
%         imdata_line(y,x) = 1;
%     end

%     for x = 1:column
%         y = round_rbar;
%         imdata_line(y,x) = 1;
%     end
    

%     hold on
%     for x = 1:column
%         y = -(x-round_cbar) + round_rbar;
%         plot(x,y,'w');
%     end

    % 중심점을 그린 선과 boundary와의 교점을 찾는다
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
