clear all
close all
clc

load answer_sheet


DataPhy200 = exceldata(2:6:end,3:end); % 200-Hz priming, physical
DataPer200 = exceldata(3:6:end,3:end); % 200-Hz priming, perceived

%% Check the accuracy
err = [];
for k = 1:size(DataPhy200,1)
    err(k,:) = abs(DataPhy200(k,:)-DataPer200(k,:));
end

%% Find the random accuracy and outliers
rnderr = [];
for k = 1:1000
    ANS = 5*floor(floor(200+200*rand(size(DataPhy200)))/5);
    rnderr(k) = sum(sum(abs(DataPhy200 - ANS),2))/(size(ANS,1)*size(ANS,2));
end

figure(1),clf,subplot(121),boxplot(err'),subplot(122),boxplot(rnderr)

outliers = find(mean(err,2) > mean(rnderr));

DataPhy200(outliers,:) = [];
DataPer200(outliers,:) = [];
NumSubj = size(DataPer200,1);

%% Check the normality of responses
Stim = unique(DataPhy200); % stimulus types
NumStim = length(Stim);

RES = reshape(DataPer200,1,size(DataPer200,1)*size(DataPer200,2));
bins = 212.5-25:25:400+12.5;
Ncnt = histc(RES,bins);
figure(2),clf,bar(bins(1:end-1)+12.5,Ncnt(1:end-1))
hold on,
plot(Stim,10*ones(NumStim),'r.-')

[h p] = kstest(RES); % normality test

%% Error distribution for each stimulus frequency
errdist = [];
for s = 1:NumStim
    dif = [];
    for k = 1:NumSubj
        tmp = find(DataPhy200(k,:) == Stim(s));
        dummy = DataPhy200(k,tmp)-DataPer200(k,tmp);
        dif = [dif dummy];
    end
    errdist(s,:) = dif;
end

figure(3),clf,
m = mean(abs(errdist),2); % mean absolute error
ss = std(abs(errdist),[],2); se = ss./sqrt(size(abs(errdist),2));
subplot(121),bar(Stim,m),hold on,plot([Stim Stim]',[m+se m-se]','r-')
m = mean(errdist,2); % mean error
ss = std(errdist,[],2); se = ss./sqrt(size(errdist,2));
subplot(122),bar(Stim,m),hold on,plot([Stim Stim]',[m+se m-se]','r-')

