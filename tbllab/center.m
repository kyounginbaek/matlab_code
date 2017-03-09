close all; clear all; clc;

alpha = 1; % �̹��� ���� �ѹ�
beta = 8; % �޾Ƶ��� ������ ����(��� �ٲ������)

divide = 36; % 360���� ��ŭ ����������

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for index = alpha : beta
    
    FILEnum = num2str(index);
    FILEname = strcat('C:\Users\xnote\Desktop\TBL LAB\Sample - label max\1 (',FILEnum,').bmp'); % imfill �Ǿ��ִ� data�� �����´�
    imdata_raw = imread(FILEname); 
    
    % ����(��),����(��) �ȼ� ���� �޾Ƶ��δ�
    [row, column] = size(imdata_raw); 
    
    % double ���·� �ٲ��ش�.
    imdata_double = im2double(imdata_raw);
    
    % �߽����� ã�´�
    [L, n] = bwlabel(imdata_double);
    [r, c] = find(L == 3);
    rbar = mean(r); % ����(��)�� �߽�
    cbar = mean(c); % ����(��)�� �߽�
    
    % �߽����� Marker�� �̿��Ͽ� �����ش�
    figure(); imshow(imdata_double);
    hold on
    for k = 1:n
        [r, c] = find(L == k);
        rbar = mean(r);
        cbar = mean(c);
        plot(cbar, rbar, 'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFacecolor', 'k', 'MarkerSize', 10);
        plot(cbar, rbar, 'Marker', '*', 'MarkerEdgeColor', 'w');
    end
        
    % �߽����� �������� ���� �׸���
    imdata_boundary = edge(imdata_raw, 'canny'); % boundary �̹����� �����Ѵ�
    imdata_line = double(imdata_boundary);
    
    % �߽���ǥ ������ ��ȯ (�ؾ� �ǳ�?)
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

    for x = 1:column % x�� ��ǥ
        for multiple = 1:divide
            y = round ( tan((2*pi/divide)*multiple) * (x - cbar) + rbar ); % y�� ��ǥ
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

    % �߽����� �׸� ���� boundary���� ������ ã�´�
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
