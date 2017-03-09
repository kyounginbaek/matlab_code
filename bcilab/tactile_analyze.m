clc;
clear all;
close all;

LC = 0.5;
HC = 50;

% EEG = pop_loadbv(pwd,'freqstimuli.vhdr');
% pop_saveset(EEG,'freqstimuli.set');

dirlist = dir('../Recording/');
dirlist = strcat('../Recording/',{dirlist(3:end).name});

nsubject = length(dirlist);

epoc_data = cell(nsubject,1);

for id = 1:length(dirlist)
    EEG = pop_loadbv(dirlist{id}, 'freqstimuli.vhdr');
    
    Fs = EEG.srate;
    channel = {EEG.chanlocs.labels};
    nchannel = length(channel); % TP9 31번 채널 제거 할까말까
    
    latency_s49 = [EEG.event(strcmp({EEG.event.type},'S 50')).latency];
    nstim = 10;
    %nstim = length(latency_s49);
    
    %data = eegfiltfft_NEW(EEG.data, Fs, LC, HC, EEG.pnts, 0, 0);
    data = EEG.data;
    epoc_data{id} = zeros(nstim,nchannel,Fs*2); % 수정 중
    
    for k = 1:nstim %length(latency_s49);
        left = latency_s49(k) - 249;
        right = latency_s49(k) + 750; % 수정 중
        epoc_data{id}(k,:,:) = data(:,left:right);
    end
    
    epoc_data{id} = epoc_data{id} - repmat(mean(epoc_data{id},2),[1,31,1]); % channel averaging
    epoc_data{id} = epoc_data{id} - repmat(mean(epoc_data{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
end

epoc_data = cell2mat(epoc_data);
nstim = 100;
%nstim = size(epoc_data,1);

nwindow = 200;
nshift = 20;
nfft = 256;

epoc_spec = zeros([nstim,nchannel,size(spectrogram(epoc_data(1,1,:),hamming(nwindow),nwindow-nshift,nfft,Fs))]);

for st = 1:nstim
    for ch = 1:nchannel
        [S,F,T] = spectrogram(epoc_data(st,ch,:),hamming(nwindow),nwindow-nshift,nfft,Fs);
        epoc_spec(st,ch,:,:) = abs(S);
    end
end

epoc_spec = epoc_spec - repmat(mean(epoc_spec(:,:,:,T<0.5),4),[1, 1, 1, size(epoc_spec,4)]);

%% C3, C4 theta, alpha, beta graph plot
% 
% channel = 7;
% 
% alpha = and(F>8,F<12);
% figure();
% plot(T-0.5, squeeze(mean(trimmean(epoc_spec(:,channel,alpha,:),10,1),3)),'r-o');
% hold on; ylim([-40,15]); 
% xlabel('time(s)'); ylabel('power');
% 
% low_beta = and(F>13,F<20);
% plot(T-0.5, squeeze(mean(trimmean(epoc_spec(:,channel,low_beta,:),10,1),3)),'g-*');
% 
% high_beta = and(F>21,F<30);
% plot(T-0.5, squeeze(mean(trimmean(epoc_spec(:,channel,high_beta,:),10,1),3)),'b-s');
% legend('alpha(8-12Hz)','low beta(13-20Hz)','high beta(21-30Hz)');
% 
% Line = line([0,0],[-50,50],'color','k'); set(Line,'linewidth',1);


%% Spectrogram outlier removed, only C3, C4

% FIDX = and(F>5,F<40);
% figure(8);
% pcolor(T-0.5,F(FIDX),squeeze(trimmean(epoc_spec(:,7,FIDX,:),10,1))); shading interp; colorbar; % 7번 = C3. 23번 = C4, 31번 = TP9(안 쓸 예정, 고장)
% title('C3'); xlabel('time(s)'); ylabel('frequency(Hz)'); % caxis([-20, 20]); 
% hold on; Line = line([0,0],[0,50],'color','k'); set(Line,'linewidth',2);
% set(gca,'FontSize',11)
% set(findall(gcf,'type','text'),'FontSize',11)
% caxis([-20, 20]);
% 
% FIDX = and(F>5,F<40);
% figure(9);
% pcolor(T-0.5,F(FIDX),squeeze(trimmean(epoc_spec(:,23,FIDX,:),10,1))); shading interp; colorbar; % 7번 = C3. 23번 = C4, 31번 = TP9(안 쓸 예정, 고장)
% title('C4'); xlabel('time(s)'); ylabel('frequency(Hz)'); % caxis([-20, 20]); 
% hold on; Line = line([0,0],[0,50],'color','k'); set(Line,'linewidth',2);
% set(gca,'FontSize',11)
% set(findall(gcf,'type','text'),'FontSize',11)
% caxis([-20, 20]);

%% Spectrogram outlier removed, individual loop

% FIDX = and(F>5,F<40);
% figure();
% for id =1:10
%     subplot(2,5,id)
%     pcolor(T-0.5,F(FIDX),squeeze(mean(epoc_spec(1+10*(id-1):10+10*(id-1),24,FIDX,:),1))); shading interp; colorbar; % 7번 = C3. 23번 = C4, 31번 = TP9(안 쓸 예정, 고장)
%     title('C3'); xlabel('time(s)'); ylabel('frequency(Hz)'); %caxis([-20, 20]); 
%     hold on; Line = line([0,0],[0,50],'color','k'); set(Line,'linewidth',2);
% end

%% Topography
% for i = 1:7
%     figure(i); 
%     range = (5*i-2):(5*i+4);
%     EPOC = squeeze(mean(mean(trimmean(epoc_spec(:,1:30,11:17,range),10,1),4),3)); % [stim, channel, Freqeuncy, time] 순서
%     % -0.2~0초(3:9) 0초~0.2초(8:14) 0.2초~0.4초(13:19) 0.4초~0.6초(18:24)
%     % 0.6초~0.8초(23:29) 0.8초~1.0초(28:34) 1.0초~1.2초(33:39) 1.2초~1.4초(38:44)
%     % 1.4초~1.6초(43:49) 1.6초~1.8초(48:54) 1.8초~2.0초(53:59) 2.0초~2.2초(58:64)
%     topoplot(EPOC,'eeglab_chan32_BKI.locs'); colorbar; caxis([-20, 20]);
% end

%% Topography outlier(T8, T7) deletion

% epoc_spec = epoc_spec(:,setdiff(1:size(epoc_spec,2), [8,24]),:,:);
% FIDX = and(F>21, F<30);
% 
% for i = 1:7
%     figure(i); 
%     range = (5*i-2):(5*i+4);
% %     EPOC = squeeze(mean(mean(trimmean(epoc_spec(:,1:28,11:17,range),10,1),4),3)); % [stim, channel, Freqeuncy, time] 순서
%     EPOC = squeeze(trimmean(mean(mean(epoc_spec(:,1:28,FIDX,range),1),4),10,3));
% 
%     % -0.2~0초(3:9) 0초~0.2초(8:14) 0.2초~0.4초(13:19) 0.4초~0.6초(18:24)
%     % 0.6초~0.8초(23:29) 0.8초~1.0초(28:34) 1.0초~1.2초(33:39) 1.2초~1.4초(38:44)
%     % 1.4초~1.6초(43:49) 1.6초~1.8초(48:54) 1.8초~2.0초(53:59) 2.0초~2.2초(58:64)
%     topoplot(EPOC,'eeglab_chan32_BKI_sample.locs'); colorbar; 
% %     caxis([-10, 10]);
% end

%% hist check
% figure();
% i = 5;
% range = (5*i-2):(5*i+4);
% hist(mean(mean(epoc_spec(:,24,11:17,range),3),4),100)