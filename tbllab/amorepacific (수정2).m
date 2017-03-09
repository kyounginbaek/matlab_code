
close all; clc;
clear all;

for name=806:851
%% Pre processing

text=strcat('1 (',num2str(name),').bmp');
imdata_raw = imread(text); % read image
fclose('all');

imdata_gray = rgb2gray(imdata_raw); % convert RGB image to gray
imdata_double = im2double(imdata_gray);

imdata_ROI = imdata_double(5:400,610:695); % cut ROI                        type in the right section
[r,c] =  size(imdata_ROI);
% figure(); imshow(imdata_ROI);
name1=strcat(num2str(name),'_1.bmp');
imwrite(imdata_ROI,name1)
imdata_ROI3=imdata_ROI;

imdata_intensity_first = zeros(r,c); % initialize 3 parts
imdata_intensity_second = zeros(r,c);
imdata_intensity_third = imdata_ROI;

%% 아래부터는 밝기를 나누는 작업

for i = 1:r % class 1
    for j = 1:c
        if imdata_ROI(i,j) > 5/6                                            % input the intensity
            imdata_intensity_first(i,j) = 1;
        end
    end
end
% figure('Name','class1'); imshow(imdata_intensity_first);
% name1=strcat('group1(',num2str(name),').bmp');
% imwrite(imdata_intensity_first,name1)


for ii = 1:r % class 2
    for jj = 1:c
        if imdata_ROI(ii,jj) > 1/2 && imdata_ROI(ii,jj) < 5/6               % input the intensity
            imdata_intensity_second(ii,jj) = 1;
        end
    end
end
% figure('Name','class2'); imshow(imdata_intensity_second);
% name1=strcat('group2(',num2str(name),').bmp');
% imwrite(imdata_intensity_second,name1)


for ii = 1:r % class 3
    for jj = 1:c
        if imdata_ROI(ii,jj) > 1/2
            imdata_intensity_third(ii,jj) = 0;
        end
    end
end
% figure('Name','class3'); imshow(imdata_intensity_third);
% name1=strcat('group3(',num2str(name),').bmp');
% imwrite(imdata_intensity_third,name1)

%% TO FIND THE UPPER TREND LINE

% imdata_discard=label_discard(imdata_intensity_first,23);
imdata_discard=imdata_intensity_first;
% figure(); imshow(imdata_discard); title('label discard');

imdata_top=top_boundary(imdata_discard);
% figure(); imshow(imdata_top); title('top points');
% name1=strcat('trend_top(',num2str(name),').bmp');
% imwrite(imdata_top,name1)

[imdata_top_left,top_l_point]=top_left(imdata_top);
[imdata_top_right,top_r_point]=top_right(imdata_top);
imdata_top_side=imdata_top_left+imdata_top_right;
% figure(); imshow(imdata_top_side); title('imdata top side');
% name1=strcat('trend_side(',num2str(name),').bmp');
% imwrite(imdata_top_side,name1)

imdata_discard_2=label_discard(imdata_top,3);
% figure(); imshow(imdata_discard_2); title('label discard 2');
% name1=strcat('trend_discard(',num2str(name),').bmp');
% imwrite(imdata_discard,name1)

imdata_trend_base=imdata_discard_2+imdata_top_side;
% figure(); imshow(imdata_trend_base); title('trend base');
% name1=strcat('trend_base(',num2str(name),').bmp');
% imwrite(imdata_trend_base,name1)


%% get column vector x, y for cubic curve fitting

[points_y,points_x]=find(imdata_trend_base==1);
% 
% i=1;
% 
% 
% for index=1:c
%     for jndex=1:r
%         if(imdata_trend_base(jndex,index)==1)
%             points_x(i)=index; % column vector : x
%             points_y(i)=jndex; % column vector : y
%             i=i+1;
%             break;
%         end
%     end
% end

%% make interpolation curve

xx=top_l_point(2):1:top_r_point(2);
yy = interp1(points_x,points_y,xx,'cubic');

for i=1:length(xx) % designate distance : trend-bottom -> trend_dist
    j=i+top_l_point(2)-1;
    trend_dist(j)=yy(i); % j : column number / trend_dist(j) : corresponding distance
end

imdata_spline=zeros(r,c);

for index=1:top_r_point(2)-top_l_point(2)+1
    
    if index>1&&abs(floor(yy(index))-floor(yy(index-1)))>1
        if floor(yy(index))<floor(yy(index-1)) % gradient : positive
            imdata_spline(floor(yy(index)):floor(yy(index-1))-1,index+top_l_point(2)-1)=1;
        else                                   % gradient : negative
            imdata_spline(floor(yy(index-1))+1:floor(yy(index)),index+top_l_point(2)-1)=1;
        end
    else
        imdata_spline(floor(yy(index)),index+top_l_point(2)-1)=1;
    end
    
end

% figure('Name','curvefitting'); imshow(imdata_spline);
[spline_r,spline_c]=find(imdata_spline==1);

a=imdata_spline+imdata_ROI;
% figure('Name','trend+ROI'); imshow(a);   
% name1=strcat('trend_interpolation(',num2str(name),').bmp');
% imwrite(a,name1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% imfill lower than trend line
imfill_undercurve=zeros(r,c);
for j=1:c
    for i=1:r
        if(imdata_spline(i,j)==1)
        imfill_undercurve(i:r,j)=1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imdata_intensity_second2=imdata_intensity_second;
imdata_intensity_third2=imdata_intensity_third;
%% remove noise upper than trend line
for index=1:top_r_point(2)-top_l_point(2)+1
    jndex=floor(yy(index));
    imdata_ROI(1:jndex-1,top_l_point(2)+index-1)=0;
    imdata_intensity_second2(1:jndex-1,top_l_point(2)+index-1)=0;
    imdata_intensity_third2(1:jndex-1,top_l_point(2)+index-1)=0;
end

% figure('Name','ROI noise'); imshow(imdata_ROI);
% name1=strcat('noise_ROI1(',num2str(name),').bmp');
% imwrite(imdata_ROI,name1)

% figure('Name','class2 noise'); imshow(imdata_intensity_second2);
% name1=strcat('noise_group2(',num2str(name),').bmp');
% imwrite(imdata_intensity_second2,name1)

t3=graythresh(imdata_intensity_third2);
imdata_intensity_third2=im2bw(imdata_intensity_third2,t3);
% figure('Name','class 3 noise'); imshow(imdata_intensity_third2);
% name1=strcat('noise_group3(',num2str(name),').bmp');
% imwrite(imdata_intensity_third2,name1)

imdata_sum=imdata_intensity_first+imdata_intensity_second2+imdata_intensity_third2;
% figure('Name','sum'); imshow(imdata_sum);
% name1=strcat('noise_all(',num2str(name),').bmp');
% imwrite(imdata_sum,name1)

imdata_sum_label=label_discard(imdata_sum,10);
% figure('Name','sum discard'); imshow(imdata_sum_label);
% name1=strcat('noise_all_label(',num2str(name),').bmp');
% imwrite(imdata_sum_label,name1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imdata_intensity_otsu1(:,:)=imdata_sum_label(:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% PRE HOLE BOUNDARY
% 
% % median filter
% imdata_intensity_otsu=medfilt2(imdata_ROI,[3 3]);
% % imdata_intensity_otsu=imdata_ROI;
% 
% % binary image 
% t=graythresh(imdata_intensity_otsu);
% imdata_intensity_otsu1=im2bw(imdata_intensity_otsu,t);
% figure('Name','imdata intensity otsu1'); imshow(imdata_intensity_otsu1);
% 
% % name1=strcat('ROI_otsu(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu1,name1)

%% mask up
% 
% mask_1=strel('disk',5);
% imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% % imdata_mask_up=imdilate(imdata_mask_up,mask_1);
% % imdata_mask_up=imdilate(imdata_mask_up,mask_1);
% 
% figure('Name','mask up'); imshow(imdata_mask_up);

% name1=strcat('mask_up_ROI(',num2str(name),').bmp');
% imwrite(imdata_mask_up,name1)
% % 
% % mask_1=strel('disk',1);
% % imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% % figure('Name','mask up'); imshow(imdata_mask_up);
% % 
% % mask_2=strel('disk',3);
% % imdata_mask_down=imerode(imdata_mask_up,mask_2);
% % figure('Name','mask down'); imshow(imdata_mask_down);
% % 
% % imdata_intensity_otsu2=imdata_intensity_otsu1+imdata_mask_down;
% % figure('Name','mask in'); imshow(imdata_intensity_otsu1);
% % 
% % name1=strcat('mask_up_ROI1(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu2,name1)
% % 
% % mask_1=strel('disk',2);
% % imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% % figure('Name','mask up'); imshow(imdata_mask_up);
% % 
% % mask_2=strel('disk',4);
% % imdata_mask_down=imerode(imdata_mask_up,mask_2);
% % figure('Name','mask down'); imshow(imdata_mask_down);
% % 
% % imdata_intensity_otsu2=imdata_intensity_otsu1+imdata_mask_down;
% % figure('Name','mask in'); imshow(imdata_intensity_otsu1);
% % 
% % name1=strcat('mask_up_ROI2(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu2,name1)
% % 
% % mask_1=strel('disk',3);
% % imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% % figure('Name','mask up'); imshow(imdata_mask_up);
% % 
% % mask_2=strel('disk',5);
% % imdata_mask_down=imerode(imdata_mask_up,mask_2);
% % figure('Name','mask down'); imshow(imdata_mask_down);
% % 
% % imdata_intensity_otsu2=imdata_intensity_otsu1+imdata_mask_down;
% % figure('Name','mask in'); imshow(imdata_intensity_otsu1);
% % 
% % name1=strcat('mask_up_ROI3(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu2,name1)

mask_1=strel('disk',4);
imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% figure('Name','mask up'); imshow(imdata_mask_up);

mask_2=strel('disk',6);
imdata_mask_down=imerode(imdata_mask_up,mask_2);
% figure('Name','mask down'); imshow(imdata_mask_down);

imdata_intensity_otsu2=imdata_intensity_otsu1+imdata_mask_down;
% figure('Name','mask in'); imshow(imdata_intensity_otsu1);

% % name1=strcat('mask_up_ROI4(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu2,name1)
% % 
% % mask_1=strel('disk',5);
% % imdata_mask_up=imdilate(imdata_intensity_otsu1,mask_1);
% % figure('Name','mask up'); imshow(imdata_mask_up);
% % 
% % mask_2=strel('disk',8);
% % imdata_mask_down=imerode(imdata_mask_up,mask_2);
% % figure('Name','mask down'); imshow(imdata_mask_down);
% % 
% % imdata_intensity_otsu2=imdata_intensity_otsu1+imdata_mask_down;
% % figure('Name','mask in'); imshow(imdata_intensity_otsu1);
% % 
% % name1=strcat('mask_up_ROI5(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu2,name1)

imdata_intensity_otsu1(:,:)=imdata_intensity_otsu2(:,:);
%% mask down
% mask_2=strel('disk',8);
% imdata_mask_down=imerode(imdata_mask_up,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% % imdata_mask_down=imerode(imdata_mask_down,mask_2);
% 
% figure('Name','mask down'); imshow(imdata_mask_down);
% 
% % name1=strcat('mask_down_ROI(',num2str(name),').bmp');
% % imwrite(imdata_mask_down,name1)

%% add to imdata_ROI
% imdata_intensity_otsu1=imdata_intensity_otsu1+imdata_mask_down;
% figure('Name','mask in'); imshow(imdata_intensity_otsu1);
% 
% % name1=strcat('mask_in_ROI(',num2str(name),').bmp');
% % imwrite(imdata_intensity_otsu1,name1)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INNER HOLE
imdata_fill1=imfill(imdata_intensity_otsu1,'hole');
imdata_holes=imdata_fill1-imdata_intensity_otsu1;
% figure('Name','fill'); imshow(imdata_holes);

% name1=strcat('inner_hole_fill(',num2str(name),').bmp');
% imwrite(imdata_holes,name1)

[ha,hb]=label_max(imdata_holes);
imdata_holes_discard=label_discard(imdata_holes,length(hb)*0.68);
% figure('Name','discard holes'); imshow(imdata_holes_discard);

imdata_holes_discard=label_discard(imdata_holes_discard,200);
% figure('Name','discard holes2'); imshow(imdata_holes_discard);

% inner hole sorting
standard=0.5;% n==0

[r,c]=size(imdata_holes_discard);
result=zeros(r,c);
[l,n_hole]=bwlabel(imdata_holes_discard);
c_region=[0,0];

if n_hole>0
    
    if n_hole>1
        min_distance=r;
        for i=1:n_hole
            [r_data,c_data]=find(l==i);
            mean_r=floor(mean(r_data));
            mean_c=floor(mean(c_data));
            
            new_distance=abs(mean_r-trend_dist(mean_c));
            
            if min_distance>new_distance
                min_distance=new_distance;
                
                region=min(abs(min(c_data)-mean_c),abs(max(c_data)-mean_c));
                c_region=[mean_c-region,mean_c+region]; %#ok<NASGU>
                min_inner_hole=i;
     
            end
            
        end
        
        [loca_data]=find(l==min_inner_hole);
        for j=1:length(loca_data)
            result(loca_data(j))=1;
        end
        
%             figure(); imshow(result);
            
        [r_data,c_data]=find(l==min_inner_hole);
        mean_r=floor(mean(r_data));
        mean_c=floor(mean(c_data));
        
        region=min(abs(min(c_data)-mean_c),abs(max(c_data)-mean_c));
        c_region=[mean_c-region,mean_c+region];
        
        if abs(mean_r-trend_dist(mean_c))>abs(r-trend_dist(mean_c))*standard
            result(:,:)=0;
            c_region=[0,0];
        end
        
        
%             figure(); imshow(result);
        
    else
        [r_data,c_data]=find(l==1);
        mean_r=floor(mean(r_data));
        mean_c=floor(mean(c_data));
        
        region=min(abs(min(c_data)-mean_c),abs(max(c_data)-mean_c));
        c_region=[mean_c-region,mean_c+region];
        
        if abs(mean_r-trend_dist(mean_c))>abs(r-trend_dist(mean_c))*standard
            result(:,:)=0;
            c_region=[0,0];
        end
    end
            
end


% name1=strcat('inner_hole_discard(',num2str(name),').bmp');
% imwrite(imdata_holes_discard,name1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 반전
imdata_intensity_otsu=~imdata_intensity_otsu1;
% figure(); imshow(imdata_intensity_otsu); title('imdata intensity otsu');

% name1=strcat('in_out_ROI(',num2str(name),').bmp');
% imwrite(imdata_intensity_otsu,name1)

%% bwtraceboundary
 
imdata_bwtraceboundary=zeros(r,c);

left_tmp = find(imdata_intensity_otsu(:,1) == 1);
left=min(left_tmp);

bwtraceboundary_location=bwtraceboundary(imdata_intensity_otsu,[left,1],'E');
bwtraceboundary_location_size=size(bwtraceboundary_location);

for index=1:bwtraceboundary_location_size(1,1)
    imdata_bwtraceboundary(bwtraceboundary_location(index,1),bwtraceboundary_location(index,2))=1;
end

% figure(); imshow(imdata_bwtraceboundary); title('imdata bwtraceboundary');
% name1=strcat('bwtraceboundary(',num2str(name),').bmp');
% imwrite(imdata_bwtraceboundary,name1)

%% (1) two boundary

two_boundary=imdata_spline+imdata_bwtraceboundary;
% figure(); imshow(two_boundary); title('two boundary');
% name1=strcat('two_boundary(',num2str(name),').bmp');
% imwrite(two_boundary,name1)

two_boundary(1,:)=0;
two_boundary(r,:)=0;
two_boundary(:,1)=0;
two_boundary(:,c)=0;

two_boundary_fill=imfill(two_boundary,'holes');
% figure(); imshow(two_boundary_fill); title('two boundary fill');
% name1=strcat('two_boundary_fill(',num2str(name),').bmp');
% imwrite(two_boundary_fill,name1)

se_2bound=strel('disk',2);
two_boundary_open=imopen(two_boundary_fill,se_2bound);
% figure(); imshow(two_boundary_open); title('two boundary open');
% name1=strcat('two_boundary_fill_section(',num2str(name),').bmp');
% imwrite(two_boundary_open,name1)

[l,n]=bwlabel(two_boundary_open);
if n>=1
    [a2,b2]=label_max(two_boundary_open);
    [a2_min,b2_min]=label_min(two_boundary_open);

    if length(b2)-length(b2_min)>600
        two_boundary_large=a2; % the most largest : hole
    else
        two_boundary_large=two_boundary_open;
        [l,n]=bwlabel(two_boundary_open);
    
        if c_region~=[0,0]
        for i=1:n
            [r_l,c_l]=find(l==i);
            mean_r=mean(r_l);
            mean_c=mean(c_l);
            if mean_c>c_region(2)||mean_c<c_region(1)
                for j=1:length(r_l)
                    two_boundary_large(r_l(j),c_l(j))=0;
                end
            end
        end
        else
            two_boundary_large=zeros(r,c);
        end
    
    end
else
    two_boundary_large(:,:)=two_boundary_open(:,:);
end
    
    % center of each label
    % compare with min_distance_c, min_distance_c_region
    % if in, then okay!
    %%%%%%%%%%%%%%%%%%%%%%%
% %     two_boundary_large=label_discard(two_boundary_open,length(b2)*0.68);
% figure(); imshow(two_boundary_large); title('2 boundary large');
% name1=strcat('two_boundary_discard(',num2str(name),').bmp');
% imwrite(two_boundary_large,name1)

% two_boundary_discard=label_discard(two_boundary_open,100);
% figure(); imshow(two_boundary_discard); title('2 boundary discard');

two_boundary_large3=two_boundary_large+imdata_holes_discard;
% figure(); imshow(two_boundary_large3);
% two_boundary_large2=zeros(r,c,3);
% two_boundary_large2(:,:,1)=two_boundary_large3;
% 
% %% compare size ( outer vs inner )
% max_size(name)=max(length(a2),length(hb));
% sorting_image_data(:,:,name)=two_boundary_large3;

%% show result
imdata_ROI2=zeros(r,c,3);
% imdata_ROI2 = imdata_raw
imdata_ROI2(:,:,1)=imdata_ROI3;
imdata_ROI2(:,:,2)=imdata_ROI3;
imdata_ROI2(:,:,3)=imdata_ROI3;

for index = 1:r
    for jndex = 1:c
        if (two_boundary_large3(index, jndex) == 1)
            imdata_ROI2(index, jndex, 1) = 1;
            imdata_ROI2(index, jndex, 2) = 0;
            imdata_ROI2(index, jndex, 3) = 0;
        end
    end
end
    
    
result=imdata_ROI2;
% figure(); imshow(result);
name1=strcat(num2str(name),'_2.bmp');
imwrite(result,name1)

end
% % 
% % max_max_size=max(max_size); % find max sized pixel이 있는 image number
% % max_image_number=find(max_size==max_max_size);
% % % 
% % % start_image<-max_image_number->end_image
% % 
% % % center 찾기
% % 
% % standard_image=sorting_image_data(max_image_number);
% % [r,c]=find(standard_image==1);
% % center_r=mean(r);
% % center_c=mean(c);
% % 
% % range_c=max(abs(c-mean(c)));
% % 
% % % max pixel 




