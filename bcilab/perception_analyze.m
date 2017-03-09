clc;
clear all;
close all;

%% Excel file load

exceldata = xlsread('../answer_sheet.xlsx');
stimdata200 = zeros(10,24);
stimdata300 = zeros(10,24);
stimdata400 = zeros(10,24);

behabdata200 = zeros(10,24);
behabdata300 = zeros(10,24);
behabdata400 = zeros(10,24);

for i = 1:10
    stimdata200(i,:) = exceldata(2+6*(i-1),3:26);
    behabdata200(i,:) = exceldata(3+6*(i-1),3:26);
    
    stimdata300(i,:) = exceldata(4+6*(i-1),3:26);
    behabdata300(i,:) = exceldata(5+6*(i-1),3:26);
    
    stimdata400(i,:) = exceldata(6+6*(i-1),3:26);
    behabdata400(i,:) = exceldata(7+6*(i-1),3:26);
end

%% EEG data load

dirlist = dir('../Recording/');
dirlist = strcat('../Recording/',{dirlist(3:end).name});

nsubject = length(dirlist);
epoc_data = cell(nsubject,1);

for id = 1:length(dirlist)
    EEG = pop_loadbv(dirlist{id}, 'perception.vhdr');
    
    Fs = EEG.srate;
    channel = {EEG.chanlocs.labels};
    nchannel = length(channel); % TP9 31번 채널 제거 할까말까
    
    latency_s49 = [EEG.event(strcmp({EEG.event.type},'S 49')).latency];
    latency_s50 = [EEG.event(strcmp({EEG.event.type},'S 50')).latency];
    nstim = 72;
    %nstim = length(latency_s49);
    
    %data = eegfiltfft_NEW(EEG.data, Fs, LC, HC, EEG.pnts, 0, 0);
    data = EEG.data;
    
    epoc_data49{id} = zeros(nstim,nchannel,Fs*2); % 수정 중
    epoc_data50{id} = zeros(nstim,nchannel,Fs*2);
    
    for k = 1:nstim %length(latency_s49);
        left49 = latency_s49(k) - 249;
        right49 = latency_s49(k) + 750; % 수정 중
        epoc_data49{id}(k,:,:) = data(:,left49:right49);
        
        left50 = latency_s50(k) - 249;
        right50 = latency_s50(k) + 750;
        epoc_data50{id}(k,:,:) = data(:, left50:right50);
    end
    
    epoc_data49{id} = epoc_data49{id} - repmat(mean(epoc_data49{id},2),[1,31,1]); % channel averaging
    epoc_data49{id} = epoc_data49{id} - repmat(mean(epoc_data49{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
    
    epoc_data50{id} = epoc_data50{id} - repmat(mean(epoc_data50{id},2),[1,31,1]); % channel averaging
    epoc_data50{id} = epoc_data50{id} - repmat(mean(epoc_data50{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
    
end

%% Excel behavior data analysis

% close all;
% 
% figure(1);
% X200 = [225 250 275 300 325 350];
% dummy = zeros(1,40);
% Y200 = zeros(1,6);
% for id = 1:6
%     [r,c] = find(stimdata200==X200(1,id));
%     for i = 1:40 
%         dummy(1,i) = behabdata200(r(i),c(i));
%     end
%     Y200(1,id) = trimmean(dummy,10);
% end
% 
% plot(X200, X200,'r-*');
% hold on 
% plot(X200, Y200,'b-*');
% title('Standard frequency = 200Hz'); xlabel('suggested frequency(Hz)'); ylabel('frequency(Hz)');
% legend('suggested frequency','perceived frequency'); ylim([200 400]);
% 
% %%
% 
% figure(2);
% X300 = [225 250 275 325 350 375];
% dummy = zeros(1,40);
% Y300 = zeros(1,6);
% for id = 1:6
%     [r,c] = find(stimdata300==X300(1,id));
%     for i = 1:40 
%         dummy(1,i) = behabdata300(r(i),c(i));
%     end
%     Y300(1,id) = trimmean(dummy,10);
% end
% 
% plot(X300, X300,'r-*');
% hold on 
% plot(X300, Y300,'b-*');
% title('Standard frequency = 300Hz'); xlabel('suggested frequency(Hz)'); ylabel('frequency(Hz)');
% legend('suggested frequency','perceived frequency'); ylim([200 400]);
% 
% %%
% 
% figure(3);
% X400 = [250 275 300 325 350 375];
% dummy = zeros(1,40);
% Y400 = zeros(1,6);
% for id = 1:6
%     [r,c] = find(stimdata400==X400(1,id));
%     for i = 1:40 
%         dummy(1,i) = behabdata400(r(i),c(i));
%     end
%     Y400(1,id) = trimmean(dummy,10);
% end
% 
% plot(X400, X400,'r-*');
% hold on 
% plot(X400, Y400,'b-*');
% title('Standard frequency = 400Hz'); xlabel('suggested frequency(Hz)'); ylabel('frequency(Hz)');
% legend('suggested frequency','perceived frequency'); ylim([200 400]);

%% EEG data analysis

