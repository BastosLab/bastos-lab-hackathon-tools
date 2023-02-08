%Good one: Lucky-Search-05182018-001

cfg = [];
cfg.viewmode = 'vertical';
% cfg.channel = 1:16;
% cfg.channel = 17:32;
cfg.channel = 33:48;
ft_databrowser(cfg, data)

U1 = 1:16; %Area PFC
U2 = 17:32; %Area 7A 
U3 = 33:48; %Area V4 Foveal


%interpolate bad channels

% badchan1 = 7;
% badchan2 = 9;
% badchan3 = [];
% 
% lfp(:,U1(badchan1),:) = (lfp(:,U1(badchan1)+1,:) + lfp(:,U1(badchan1)-1,:))./2;
% lfp(:,U2(badchan2),:) = (lfp(:,U2(badchan2)+1,:) + lfp(:,U2(badchan2)-1,:))./2;
% lfp(:,U3(badchan3),:) = (lfp(:,U3(badchan3)+1,:) + lfp(:,U3(badchan3)-1,:))./2;

figure; 
subplot(1,3,1); imagesc(squeeze(mean(lfp(:,U1,:),3))'); caxis([-0.1 .1])
subplot(1,3,2); imagesc(squeeze(mean(lfp(:,U2,:),3))'); caxis([-0.1 .1])
title('LFP ERP')
subplot(1,3,3); imagesc(squeeze(mean(lfp(:,U3,:),3))'); caxis([-0.1 .1])


data =[];
for t = 1:size(lfp, 3)
    data.time{t,1} = -1.5:0.001:3.0;
    data.trial{t,1} = zeros(size(electrodeInfo,1), 4501); %num channel x num time
    data.trial{t,1} = squeeze(lfp(:,:,t))';
end
data.fsample = 1000;
data.label = {};
%make each channel unique
for c = 1:size(electrodeInfo,1)
    data.label{c} = [electrodeInfo.area{c} int2str(electrodeInfo.channel(c))];
end
data.trialinfo = zeros(size(lfp, 3), 1);
data.trialinfo(:,1) = 1:size(lfp, 3);

cfg = [];
cfg.latency = [-1 1];
data = ft_selectdata(cfg, data)

cfg = [];
% cfg.trials = 1:100;
data = ft_selectdata(cfg, data)


cfg = [];
cfg.method = 'mtmfft';
% cfg.tapsmofrq = 5;
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.keeptrials = 'no';
cfg.trials = 1:100;
cfg.foi = [0:150];
mpow = ft_freqanalysis(cfg, data)


% mpow = ft_freqdescriptives([], pow);


relpow1 = mpow.powspctrm(U1,:) ./ repmat(max(mpow.powspctrm(U1,:)), [16 ,1]);
relpow2 = mpow.powspctrm(U2,:) ./ repmat(max(mpow.powspctrm(U2,:)), [16 ,1]);
relpow3 = mpow.powspctrm(U3,:) ./ repmat(max(mpow.powspctrm(U3,:)), [16,1]);

figure; 
subplot(1,3,1); imagesc(relpow1)
subplot(1,3,2); imagesc(relpow2)
title('LFP relative power')
subplot(1,3,3); imagesc(relpow3)


b1 = nearest(mpow.freq, 8); b2 = nearest(mpow.freq, 30);
t1 = nearest(mpow.freq, 1); t2 = nearest(mpow.freq, 4);
g1 = nearest(mpow.freq, 50); g2 = nearest(mpow.freq, 150);

figure; 
subplot(1,3,1);
plot(mean(relpow1(:,t1:t2),2), 1:16, 'k'); hold on;
plot(mean(relpow1(:,b1:b2),2), 1:16, 'r'); hold on;
plot(mean(relpow1(:,g1:g2),2), 1:16, 'g'); hold on;
set(gca, 'ydir', 'reverse')
title('PFC')
subplot(1,3,2);
plot(mean(relpow2(:,t1:t2),2), 1:16, 'k'); hold on;
plot(mean(relpow2(:,b1:b2),2), 1:16, 'r'); hold on;
plot(mean(relpow2(:,g1:g2),2), 1:16, 'g'); hold on;
set(gca, 'ydir', 'reverse')
title('7A')
subplot(1,3,3);
plot(mean(relpow3(:,t1:t2),2), 1:16, 'k'); hold on;
plot(mean(relpow3(:,b1:b2),2), 1:16, 'r'); hold on;
plot(mean(relpow3(:,g1:g2),2), 1:16, 'g'); hold on;
set(gca, 'ydir', 'reverse')
title('V4')

epoched_U1 = lfp(:,U1,:);
epoched_U2 = lfp(:,U2,:);
epoched_U3 = lfp(:,U3,:);

NSampl=4501; %number of time samples per trial
NCh=16;
NTrl=size(lfp,3);

sig=0.4; %cortical conductivity
dis=200.*10.^-3; %inter-site distance in micrometers
sep = 2; %separation between channels
CSD1=zeros(NSampl,NCh-sep*2,NTrl);
CSD2=zeros(NSampl,NCh-sep*2,NTrl);
CSD3=zeros(NSampl,NCh-sep*2,NTrl);
for ch = sep+1:NCh-sep
    CSD1(:,ch-sep,:) = -sig*(epoched_U1(:,ch-sep,:)-(2.*epoched_U1(:,ch,:))+epoched_U1(:,ch+sep,:))/(((dis).*sep).^2);
    CSD2(:,ch-sep,:) = -sig*(epoched_U2(:,ch-sep,:)-(2.*epoched_U2(:,ch,:))+epoched_U2(:,ch+sep,:))/(((dis).*sep).^2);
    CSD3(:,ch-sep,:) = -sig*(epoched_U3(:,ch-sep,:)-(2.*epoched_U3(:,ch,:))+epoched_U3(:,ch+sep,:))/(((dis).*sep).^2);
end

mCSD1 = squeeze(mean(CSD1(:,:,:),3)); %1000:2500
mCSD2 = squeeze(mean(CSD2(:,:,:),3));
mCSD3 = squeeze(mean(CSD3(:,:,:),3));

figure; 
subplot(1,3,1); imagesc(mCSD1');  caxis([-0.05 .05])
subplot(1,3,2); imagesc(mCSD2');  caxis([-0.025 .025])
title('CSD')
subplot(1,3,3); imagesc(mCSD3');  caxis([-0.15 .15])

figure; 
subplot(1,3,1); imagesc(mCSD1'); xlim([1400 2500]); caxis([-0.05 .05])
subplot(1,3,2); imagesc(mCSD2'); xlim([1400 2500]); caxis([-0.025 .025])
title('CSD stimulus')
subplot(1,3,3); imagesc(mCSD3'); xlim([1400 2500]); caxis([-0.15 .15])


figure; 
subplot(1,3,1); imagesc(mCSD1'); xlim([3600 4000]); caxis([-0.05 .05])
subplot(1,3,2); imagesc(mCSD2'); xlim([3600 4000]); caxis([-0.025 .025])
title('CSD test on')
subplot(1,3,3); imagesc(mCSD3'); xlim([3600 4000]); caxis([-0.15 .15])

figure; 
subplot(1,3,1); imagesc(mCSD1'); xlim([3800 4200]); caxis([-0.05 .05])
subplot(1,3,2); imagesc(mCSD2'); xlim([3800 4200]); caxis([-0.025 .025])
title('CSD response locked')
subplot(1,3,3); imagesc(mCSD3'); xlim([3800 4200]); caxis([-0.15 .15])



