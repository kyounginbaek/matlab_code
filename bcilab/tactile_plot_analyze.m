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

epoc_data_f = cell(nsubject,1);
epoc_data_p = cell(nsubject,1);
epoc_data_b = cell(nsubject,1);

for id = 1:length(dirlist)-1
    EEG_f = pop_loadbv(dirlist{id}, 'freqstimuli.vhdr');
    EEG_p = pop_loadbv(dirlist{id}, 'pressure.vhdr');
    
    Fs = EEG_f.srate;
    channel = {EEG_f.chanlocs.labels};
    nchannel = length(channel); % TP9 31번 채널 제거 할까말까
    
    latency_s50_f = [EEG_f.event(strcmp({EEG_f.event.type},'S 50')).latency];
    latency_s52_f = [EEG_f.event(strcmp({EEG_f.event.type},'S 52')).latency];
    
    latency_s49_p = [EEG_p.event(strcmp({EEG_p.event.type},'S 49')).latency];
    latency_s52_p = [EEG_p.event(strcmp({EEG_p.event.type},'S 52')).latency];
    
    nstim_f = 10;
    nstim_p = 8;
    nstim_b = 2;
    %nstim = length(latency_s49);
    
    data_f = eegfiltfft_NEW(EEG_f.data, Fs, LC, HC, EEG_f.pnts, 0, 0);
    data_p = eegfiltfft_NEW(EEG_p.data, Fs, LC, HC, EEG_p.pnts, 0, 0);
    
    epoc_data_f{id} = zeros(nstim_f,nchannel,Fs*2); % 수정 중
    epoc_data_p{id} = zeros(nstim_p,nchannel,Fs*2);
    epoc_data_b{id} = zeros(nstim_b,nchannel,Fs*2);
    
    for k = 1:nstim_f %length(latency_s49);
        left = latency_s50_f(k) - (Fs*0.5-1);
        right = latency_s50_f(k) + (Fs*1.5); % 수정 중
        
        epoc_data_f{id}(k,:,:) = data_f(:,left:right);
    end
    
    for k = 1:nstim_p %length(latency_s49);
        left = latency_s49_p(k) - (Fs*0.5-1);
        right = latency_s49_p(k) + (Fs*1.5); % 수정 중
        
        epoc_data_p{id}(k,:,:) = data_p(:,left:right);
    end
    
    %% baseline epoching
%     left_b = latency_s52_f(1);
%     right_b = latency_s52_f(1) + 999; % 수정 중
%     epoc_datab{id}(1,:,:) = data(:,leftb:rightb);

%     left_b = latency_s52_p(2);
%     right_b = latency_s52_p(2) + 999;
%     epoc_data_b{id}(2,:,:) = data_b(:,left2b:right2b);
    
    %% modulating
    
    epoc_data_f{id} = epoc_data_f{id} - repmat(mean(epoc_data_f{id},2),[1,31,1]); % channel averaging
%     epoc_data_f{id} = epoc_data_f{id} - repmat(mean(epoc_data_f{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
    
    epoc_data_p{id} = epoc_data_p{id} - repmat(mean(epoc_data_p{id},2),[1,31,1]); % channel averaging
%     epoc_data_p{id} = epoc_data_p{id} - repmat(mean(epoc_data2{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
    
%     epoc_data_b{id} = epoc_datab{id} - repmat(mean(epoc_datab{id},2),[1,31,1]); % channel averaging
%     epoc_data_b{id} = epoc_datab{id} - repmat(mean(epoc_datab{id}(:,:,1:250),3),[1,1,1000]); % baseline correction
    
end

epoc_data_f = cell2mat(epoc_data_f);
epoc_data_p = cell2mat(epoc_data_p);
% epoc_data_b = cell2mat(epoc_datab);
nstim_f = 90;
nstim_p = 72;
% nstim_b = 20;
%nstim = size(epoc_data,1);

nwindow = 400;
nshift = 40; % !! oscillation 줄이기
nfft = 256;

epoc_spec_f = zeros([nstim_f,nchannel,size(spectrogram(epoc_data_f(1,1,:),hamming(nwindow),nwindow-nshift,nfft,Fs))]);
epoc_spec_p = zeros([nstim_p,nchannel,size(spectrogram(epoc_data_p(1,1,:),hamming(nwindow),nwindow-nshift,nfft,Fs))]);
% epoc_spec_b = zeros([nstim_b,nchannel,size(spectrogram(epoc_data_b(1,1,:),hamming(nwindow),nwindow-nshift,nfft,Fs))]);

for st = 1:nstim_f
    for ch = 1:nchannel
        [Sf,F,T,Pf] = spectrogram(epoc_data_f(st,ch,:),hamming(nwindow),nwindow-nshift,nfft,Fs);
        epoc_spec_f(st,ch,:,:) = Pf;
    end
end

for st = 1:nstim_p
    for ch = 1:nchannel
        [Sp,F,T,Pp] = spectrogram(epoc_data_p(st,ch,:),hamming(nwindow),nwindow-nshift,nfft,Fs);
        epoc_spec_p(st,ch,:,:) = Pp;
    end
end

% for st = 1:nstimb
%     for ch = 1:nchannel
%         [Sb,Fb,Tb] = spectrogram(epoc_datab(st,ch,:),hamming(nwindow),nwindow-nshift,nfft,Fs);
%        
%         epoc_specb(st,ch,:,:) = abs(Sb);
%     end
% end

channel = 7;
alpha = and(F>8,F<12);
low_beta = and(F>13,F<20);
high_beta = and(F>21,F<30);

epoc_spec_f = epoc_spec_f - repmat(mean(epoc_spec_f(:,:,:,T<0.5),4),[1,1,1,size(epoc_spec_f,4)]);
% figure, imagesc(squeeze(mean(epoc_spec_f(:,channel,high_beta,:),3)))
% for k = 1:100, figure(k), plot(squeeze(mean(epoc_spec_f(k,7,alpha,:),3))), pause, end

epoc_spec_p = epoc_spec_p - repmat(mean(epoc_spec_p(:,:,:,T<0.5),4),[1, 1, 1, size(epoc_spec_p,4)]);
% figure, imagesc(squeeze(mean(epoc_spec_p(:,channel,high_beta,:),3)))
% for k = 1:100, figure(k), plot(squeeze(mean(epoc_spec_p(k,7,alpha,:),3))), pause, end


%% Spectrogram outlier removed, only C3, C4

% target = epoc_spec_f;
% 
% FIDX = and(F>5,F<40);
% figure();
% pcolor(T-0.5,F(FIDX),squeeze(trimmean(target(:,7,FIDX,:),10,1))); shading interp; colorbar; % 7번 = C3. 23번 = C4, 31번 = TP9(안 쓸 예정, 고장)
% title('C3'); xlabel('time(s)'); ylabel('frequency(Hz)'); % caxis([-20, 20]); 
% hold on; Line = line([0,0],[0,50],'color','k'); set(Line,'linewidth',2);
% set(gca,'FontSize',11)
% set(findall(gcf,'type','text'),'FontSize',11)
% caxis([-0.15, 0.15]);
% 
% FIDX = and(F>5,F<40);
% figure();
% pcolor(T-0.5,F(FIDX),squeeze(trimmean(target(:,23,FIDX,:),10,1))); shading interp; colorbar; % 7번 = C3. 23번 = C4, 31번 = TP9(안 쓸 예정, 고장)
% title('C4'); xlabel('time(s)'); ylabel('frequency(Hz)'); % caxis([-20, 20]); 
% hold on; Line = line([0,0],[0,50],'color','k'); set(Line,'linewidth',2);
% set(gca,'FontSize',11)
% set(findall(gcf,'type','text'),'FontSize',11)
% caxis([-0.15, 0.15]);

%% Topography

% target = epoc_spec_f;
% band = high_beta;
% 
% for i = 1:8
%     figure(i); 
%     range = (2*i-1):(2*i);
%     EPOC = squeeze(trimmean(mean(mean(target(:,1:30,band,range),3),4),10,1)); % [stim, channel, Freqeuncy, time] 순서
%     topoplot(EPOC,'eeglab_chan32_BKI.locs'); colorbar; 
% %     caxis([-0.6, 0.6]);
% end

%% Topography outlier(T8, T7) deletion

epoc_spec_f = epoc_spec_f(:,setdiff(1:size(epoc_spec_f,2), [8,24]),:,:);
epoc_spec_p = epoc_spec_p(:,setdiff(1:size(epoc_spec_p,2), [8,24]),:,:);

band = alpha;

% for i = 1:8
%     figure(i); 
%     range = (2*i-1):(2*i);
%     EPOC = squeeze(mean(mean(trimmean(epoc_spec(:,1:28,11:17,range),10,1),4),3)); % [stim, channel, Freqeuncy, time] 순서
    figure();
    EPOC = squeeze(trimmean(mean(epoc_spec_f(:,1:28,band,8),3),10,1));
    topoplot(EPOC,'eeglab_chan32_BKI_sample.locs'); colorbar; 
    caxis([-0.70, 0.70]);
    
    figure();
    EPOC = squeeze(trimmean(mean(epoc_spec_p(:,1:28,band,8),3),10,1));
    topoplot(EPOC,'eeglab_chan32_BKI_sample.locs'); colorbar; 
    caxis([-0.70, 0.70]);
% end

%% C3, C4 theta, alpha, beta graph plot

close all;
channel = 23;

figure();
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_f(:,channel,alpha,:),3),10,1)),'b-*');
hold on
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_p(:,channel,alpha,:),3),10,1)),'r-*');
% plot(Tb-0.5, squeeze(mean(trimmean(epoc_specb(:,channel,alphab,:),10,1),3)),'k-');
ylim([-0.25,0.15]);
Line = line([0,0],[-0.5,0.5],'color','k'); set(Line,'linewidth',1);
xlabel('time(s)'); ylabel('power'); title('alpha(8-12Hz)');
legend('frequency(250Hz)','pressure');
set(gca,'FontSize',11)
set(findall(gcf,'type','text'),'FontSize',11)

figure();
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_f(:,channel,low_beta,:),3),10,1)),'b-*');
hold on
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_p(:,channel,low_beta,:),3),10,1)),'r-*');
% plot(T-0.5, squeeze(mean(trimmean(epoc_specb(:,channel,low_betab,:),10,1),3)),'k-');
ylim([-0.06,0.04]);
Line = line([0,0],[-0.5,0.5],'color','k'); set(Line,'linewidth',1);
xlabel('time(s)'); ylabel('power'); title('low beta(13-20Hz)');
legend('frequency(250Hz)','pressure');
set(gca,'FontSize',11)
set(findall(gcf,'type','text'),'FontSize',11)

figure();
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_f(:,channel,high_beta,:),3),10,1)),'b-*');
hold on
plot(T-0.5, squeeze(trimmean(mean(epoc_spec_p(:,channel,high_beta,:),3),10,1)),'r-*');
% plot(T-0.5, squeeze(mean(trimmean(epoc_specb(:,channel,high_betab,:),10,1),3)),'k-');
ylim([-0.06,0.06]);
Line = line([0,0],[-0.5,0.5],'color','k'); set(Line,'linewidth',1);
xlabel('time(s)'); ylabel('power'); title('high beta(21-30Hz)');
legend('frequency(250Hz)','pressure');
set(gca,'FontSize',11)
set(findall(gcf,'type','text'),'FontSize',11)

%% C3, C4 band power & phase plot

figure(4);
channel = 23; % 7 = C3, 23 = C4
k = 5; % trimmean 대신 outlier 제거를 위함. 

alpha_f_data = sort(squeeze(mean(mean(epoc_spec(:,channel,alpha,9:21),3),4)));
alpha_f = mean(alpha_f_data(k+1:end-k));
alpha_f_std = std(alpha_f_data(k+1:end-k));
low_beta_f_data = sort(squeeze(mean(mean(epoc_spec(:,channel,low_beta,9:21),3),4)));
low_beta_f = mean(low_beta_f_data(k+1:end-k));
low_beta_f_std = std(low_beta_f_data(k+1:end-k));
high_beta_f_data = sort(squeeze(mean(mean(epoc_spec(:,channel,high_beta,9:21),3),4)));
high_beta_f = mean(high_beta_f_data(k+1:end-k));
high_beta_f_std = std(high_beta_f_data(k+1:end-k));

alpha_p_data = sort(squeeze(mean(mean(epoc_spec2(:,channel,alpha2,9:21),3),4)));
alpha_p = mean(alpha_p_data(k+1:end-k));
alpha_p_std = std(alpha_p_data(k+1:end-k));
low_beta_p_data = sort(squeeze(mean(mean(epoc_spec2(:,channel,low_beta2,9:21),3),4)));
low_beta_p = mean(low_beta_p_data(k+1:end-k));
low_beta_p_std = std(low_beta_p_data(k+1:end-k));
high_beta_p_data = sort(squeeze(mean(mean(epoc_spec2(:,channel,high_beta2,9:21),3),4)));
high_beta_p = mean(high_beta_p_data(k+1:end-k));
high_beta_p_std = std(high_beta_p_data(k+1:end-k));

err = [alpha_f_std alpha_p_std; low_beta_f_std low_beta_p_std; high_beta_f_std high_beta_p_std];

X = {'alpha(8-12Hz)'; 'low_beta(13-20Hz)'; 'high_beta(21-30Hz)'};
Y = [alpha_f, alpha_p; low_beta_f, low_beta_p; high_beta_f, high_beta_p];
hArray = barwitherr(err,Y);
set(gca, 'XTickLabel', X);
set(gca,'YGrid','on')
set(gca,'GridLineStyle','-')
set(hArray(1), 'FaceColor', 'b');
set(hArray(2), 'FaceColor', 'r');
ylim([-40 40]);
legend('frequency', 'pressure');
title('Freqeuncy & pressure difference');
set(gca,'FontSize',11)
set(findall(gcf,'type','text'),'FontSize',11)
% ylabel('mean value'); xlabel('frequency band');

%% 개인별 C3, C4 band power plot

