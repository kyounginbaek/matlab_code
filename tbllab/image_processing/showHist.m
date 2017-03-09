function[] = showHist(img)
% getting RGB channels
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

% ... and histograms
hR = imhist(R);
hG = imhist(G);
hB = imhist(B);

[imgHist,x] = imhist(rgb2gray(img));

figure('Name','Image histograms');    
    subplot(1,2,1); imshow(img);
    subplot(1,2,2);
        hold on;
        area(x,imgHist,'FaceColor',[0.6 0.6 0.6]);
        plot(x,hR,'r');
        plot(x,hG,'g');
        plot(x,hB,'b');
        xlim([0 255]);
    
return