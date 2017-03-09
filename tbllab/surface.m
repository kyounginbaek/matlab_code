close all; clear all; clc;

alpha = 1; % 이미지 시작 넘버
beta = 32; % 받아들인 사진의 개수

distance = 110/beta; % z축 layer간의 간격 (calibration is needed!)
b = zeros(beta,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for index_test = alpha : beta
%     
%     FILEnum_test = num2str(index_test);
%     FILEname_test = strcat('D:\Sample_박사님\andy\4th hole(white) 24h 816-855\1 (',FILEnum_test,').jpg'); % imfill 되어있는 data를 가져온다
%     imdata_raw_test = imread(FILEname_test); 
%     
%     im_check_test = imdata_raw_test(146:608, 559:641,1); % 좌표 추출 (height, width, channel)
%     im_bw_test = im2bw(im_check_test, graythresh(im_check_test)); % boundary를 위한 bw
%     
%     % boundary의 coordinates를 찾는다
%     [b_temp,I_temp] = bwboundaries(im_bw_test, 'holes'); 
%     x_temp = b_temp{1}(:,1);
%     b_size = size(x_temp);
%     b(index_test,1) = b_size(1,1);
% 
% end
% 
% scale = min(b); % boundary 사진에서 얼만큼의 pixel을 받을 것인가 (최소 pixel값을 지향)
scale = 50;

x = zeros(3,scale);
y = zeros(3,scale);
z = zeros(3,scale);

for index = alpha : beta
    
    FILEnum = num2str(index);
    FILEname = strcat('D:\Sample_박사님\andy\2nd hole(green) 48h 521-570\1 (',FILEnum,').jpg'); % imfill 되어있는 data를 가져온다
    imdata_raw = imread(FILEname); 
    
    im_check = imdata_raw(146:608, 559:641,1); % 좌표 추출 (height, width, channel)
    im_bw = im2bw(im_check, graythresh(im_check)); % boundary를 위한 bw
    
    % boundary의 coordinates를 찾아서 분리시킨다
    [b_tmp,I_tmp] = bwboundaries(im_bw, 'holes'); 
    x_tmp = b_tmp{1}(:,1);
    y_tmp = b_tmp{1}(:,2);   
    z_tmp = ones(size(x_tmp), 1) * index * round(distance) ;

    % boundary의 coordinates를 가장 작은 pixel에 맞춰서 줄인다
    step = length(x_tmp) / scale;
    for i = 1:scale
        x(index,i) = x_tmp(round(i*step));
        y(index,i) = y_tmp(round(i*step));
        z(index,i) = z_tmp(round(i*step));
    end
    
    
end

xaxis_calibration = max(max(x)); % 좌표를 똑바로 align 시키기 위함

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     figure(1); plot3(z, y, x,'.');
%     set(gca,'YDir','rev','ZDir','rev')
%     axis equal
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     figure(2); plot3(z,y,x,'b')
%     set(gca,'YDir','rev','ZDir','rev')
%     axis equal
%     hold on
%     for index = alpha : beta
%         plot3(z(index,:),y(index,:),x(index,:),'r')
%     end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     [X,Y] = meshgrid(x,y);
%     figure(3); surf(z,y,x);
%     set(gca,'YDir','rev','ZDir','rev')
%     axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     [X,Y] = meshgrid(x,y);
%     figure(4); surf(z,y,x);
%     set(gca,'YDir','rev','ZDir','rev')
%     axis vis3d tight
% %     axis off
%     axis equal
%     lighting phong
%     shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [X,Y] = meshgrid(x,y);
    figure(5); surf(z, y, xaxis_calibration - x);
    set(gca,'YDir','rev')
    axis vis3d tight
%     axis off
    axis equal
    camlight; lighting phong
    shading interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xlswrite('x_coordinates.xls', x) ;
% xlswrite('y_coordinates.xls', y) ; 
% xlswrite('z_coordinates.xls', z) ;

