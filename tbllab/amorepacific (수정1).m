function amorepacific(name)
close all; clc;

% 초기 이미지 처리 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text=strcat('1 (',num2str(name),').bmp');
imdata_raw = imread(text); % read image
fclose('all');

imdata_gray = rgb2gray(imdata_raw); % convert RGB image to gray
imdata_double = im2double(imdata_gray);

imdata_ROI = imdata_double(5:300,110:650); % cut ROI
[r,c,l] =  size(imdata_ROI);
% figure(); imshow(imdata_ROI);

imdata_intensity_first = zeros(r,c);
imdata_intensity_second = zeros(r,c);
imdata_intensity_third = imdata_ROI;

% 아래부터는 밝기를 나누는 작업 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st, 2nd brightest :  all data
for ii = 1:r % 가장 밝은 이미지 데이터 군
    for jj = 1:c
        if imdata_ROI(ii,jj) > 2/3
            imdata_intensity_first(ii,jj) = 1;
        end
    end
end
% figure('Name','intensity first class'); imshow(imdata_intensity_first);

for ii = 1:r % 중간 밝기의 이미지 데이터 군
    for jj = 1:c
        if imdata_ROI(ii,jj) > 1/3 && imdata_ROI(ii,jj) < 2/3
            imdata_intensity_second(ii,jj) = 1;
        end
    end
end
% figure(); imshow(imdata_intensity_second); title('intensity_second_class');

% 3rd brightest : otsu threshold->get data
for ii = 1:r % 가장 어두운 이미지 데이터 군
    for jj = 1:c
        if imdata_ROI(ii,jj) > 1/3
            imdata_intensity_third(ii,jj) = 0;
        end
    end
end
% figure(); imshow(imdata_intensity_third); title('intensity_third_class');

% 각각의 군을 이미지 처리 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imdata_tmp = imdata_intensity_first;
SE = strel('disk',5);
imdata_dilate = imdilate(imdata_tmp,SE);
% figure(); imshow(imdata_dilate); title('first_class_dilate');

for jj = 1:c % 가장 밝은 데이터를 기준으로 노이즈 제거
    point = find(imdata_dilate(:,jj) ~= 0);
    if length(point) ~= 0;
        point = min(point);
        imdata_intensity_second(1:point,jj) = 0;
        imdata_intensity_third(1:point,jj) = 0;
    end
end

% figure(); imshow(imdata_intensity_second);
% figure(); imshow(imdata_intensity_third);

threshold_3 = graythresh(imdata_intensity_third);

imdata_first_otsu = imdata_intensity_first;
imdata_second_otsu = imdata_intensity_second;
imdata_third_otsu = im2bw(imdata_intensity_third,threshold_3);

% figure(); imshow(imdata_third_otsu); title('3rd otsu');

%% imclose

se_third=strel('disk',1);
imdata_third_otsu=imclose(imdata_third_otsu,se_third);

% figure(); imshow(imdata_third_otsu); title('3rd otsu imclose');


imdata_intensity_otsu = imdata_first_otsu + imdata_second_otsu + imdata_third_otsu; % 밝기를 나눈 otsu 데이터

% figure(); imshow(imdata_first_otsu); title('intensity_first_class');
% figure(); imshow(imdata_second_otsu); title('intensity_second_class');
% figure(); imshow(imdata_third_otsu); title('intensity_third_class');

[label,label_num] = bwlabel(imdata_intensity_otsu);

for epsilon = 1:label_num
    d1 = find(label==epsilon);
    d = size(d1);
    if d(1,1) < 10
        for jndex = 1:d
            imdata_intensity_otsu(d1(jndex)) = 0;
        end
    end
end

% figure(); imshow(imdata_intensity_otsu);

for jj = 1:c
    left_tmp = find(imdata_intensity_otsu(:,jj) == 1);
    left = min(left_tmp);
    if length(left) == 1
        left_column=jj;
        break;
    end    
end
    
for jj = 1:c
    right_tmp = find(imdata_intensity_otsu(:,c - jj + 1) == 1);
    right = min(right_tmp);
    if length(right) == 1
        right_column=c - jj + 1;
        break;
    end
end

% 테두리 라인 긋기
imdata_intensity_otsu(left:r,1:left_column) = 1;
imdata_intensity_otsu(right:r,right_column:c) = 1;
imdata_intensity_otsu(r,:) = 1;

imdata_intensity_otsu = imfill(imdata_intensity_otsu,'holes');
% figure(); imshow(imdata_intensity_otsu);

imdata_medfilt=medfilt2(imdata_intensity_otsu,[3 3]);
% figure(); imshow(imdata_medfilt);

% 그냥 원본에서 노이즈 없는것
for ii = 1:r
    for jj = 1:c
        if imdata_intensity_otsu(ii,jj) == 0
            imdata_ROI(ii,jj) = 0;
        end
    end
end
% figure(); imshow(imdata_ROI);

% threshold = graythresh(imdata_ROI);
% 
% imdata_ROI_otsu = im2bw(imdata_ROI,threshold);
% figure(); imshow(imdata_ROI_otsu); title('otsu imdata_ROI');

%% TO FIND THE UPPER TREND LINE

imdata_discard=label_discard(imdata_intensity_first,18);
% figure(); imshow(imdata_discard); title('label_discard');

imdata_top=top_boundary(imdata_discard);
% figure(); imshow(imdata_top); title('top points');

imdata_top_left=top_left(imdata_top);
imdata_top_right=top_right(imdata_top);
imdata_top_side=imdata_top_left+imdata_top_right;
% figure(); imshow(imdata_top_side); title('imdata top side');

imdata_discard_2=label_discard(imdata_top,5);
% figure(); imshow(imdata_discard_2); title('label discard 2');

imdata_trend_base=imdata_discard_2+imdata_top_side;
% figure(); imshow(imdata_trend_base); title('trend base');

%% get column vector x, y for cubic curve fitting

i=1;

for index=1:c
    for jndex=1:r
        if(imdata_trend_base(jndex,index)==1)
            points_x(i)=index; % column vector : x
            points_y(i)=jndex; % column vector : y
            i=i+1;
            break;
        end
    end
end

%% make interpolation curve

xx=1:1:c;
yy = interp1(points_x,points_y,xx,'cubic');

imdata_spline=zeros(r,c);

for index=1:c
    imdata_spline(floor(yy(index)),index)=1;
end

a=imdata_spline+imdata_ROI;
% figure(); imshow(imdata_spline); title('imdata spline');
figure(); imshow(a); title('a');   


%% 반전 

imdata_intensity_otsu=~imdata_intensity_otsu;
figure(); imshow(imdata_intensity_otsu); title('imdata_intensity_otsu');

% %% median filter
% 
% imdata_medfilt = medfilt2(imdata_intensity_otsu,[3 3]);

%% bwtraceboundary

imdata_bwtraceboundary=zeros(r,c);

left_tmp = find(imdata_intensity_otsu(:,1) == 1);
left=min(left_tmp);

bwtraceboundary_location=bwtraceboundary(imdata_intensity_otsu,[left,1],'E');
bwtraceboundary_location_size=size(bwtraceboundary_location);

for index=1:bwtraceboundary_location_size(1,1)
    imdata_bwtraceboundary(bwtraceboundary_location(index,1),bwtraceboundary_location(index,2))=1;
end

figure(); imshow(imdata_bwtraceboundary); title('imdata bwtraceboundary');

%% (1) two boundary

two_boundary=imdata_spline+imdata_bwtraceboundary;
% figure(); imshow(two_boundary); title('two boundary');

two_boundary(1,:)=0;
two_boundary(r,:)=0;
two_boundary(:,1)=0;
two_boundary(:,c)=0;

two_boundary_fill=imfill(two_boundary,'holes');
% figure(); imshow(two_boundary_fill); title('two boundary fill');

se_2bound=strel('disk',2);
two_boundary_open=imopen(two_boundary_fill,se_2bound);
% figure(); imshow(two_boundary_open); title('two boundary open');

two_boundary_large=label_max(two_boundary_open);
figure(); imshow(two_boundary_large); title('2 boundary large');

% two_boundary_discard=label_discard(two_boundary_open,100);
% figure(); imshow(two_boundary_discard); title('2 boundary discard');


% %% (2) filter out pointy ones
% 
% se2=strel('disk',20);
% erode=imopen(imdata_intensity_otsu,se2);
% 
% erode=imdata_intensity_otsu-erode;
% % figure, imshow(erode)
% 
% imdata_first_medfilt=medfilt2(imdata_first_otsu,[3 3]);
% 
% up=zeros(1,c);
% 
% for jj = 1:c
%     up_tmp = find(imdata_first_medfilt(:,jj) == 1);
%     if(size(up_tmp)~=0)
%     up(1,jj) = min(up_tmp);
%     end
% end
%            
% for index=1:c
%     for jndex=1:r
%     if (size(up(1,index))~=0 & jndex<up(1,index))
%         erode(jndex,index)=0;
%     end
%     end
% end
% 
% figure(),imshow(erode)
% 
% [l,n]=bwlabel(erode);
% 
% for a=1:n
%     points=find(l==a);
%     if(length(points)<20)
%         for index=1:length(points)
%             erode(points(index))=0;
%         end
%     end
% end
% 
% figure(),imshow(erode)
%     
% trace=imdata_ROI+erode;
% 
% % for j = 1:1:c
% %     for i = 1:1:r
% %         if (imdata_first_otsu(i, j) == 1)
% %             trace(i, j) = 0.8;
% %         end
% %     end
% % end
% 
% figure, imshow(trace)
% % % 
% % name=strcat('result_',num2str(name),'.bmp');
% % imwrite(trace,name)
% % 
result=two_boundary_large+imdata_ROI;
figure(); imshow(result);
name=strcat('result_',num2str(name),'.bmp');
imwrite(imdata_ROI,name)


end





