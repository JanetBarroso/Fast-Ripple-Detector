% This is for pseudo online analysis to test slice viability.
% This code could expand to search SPW and Ripples in several channels


%The fist part of the code convert .brw files in to voltage data according
%to MAE/3brain formula. ::::...IMPORTANT NOTE...::: If you import data as a
% *.brw file comment ": If is a MATLAB file :" section below.
%If you import data as a *.mat file then comment *.BRW section

%* I will fix this so the program ask what kind of file are you using so
% you don't have to comment lines. Janet Barroso

clear all
close all
tic

%
% ...:::: *.BRW files ...:::
% 
% Load Files from file.brw importing all folders 
FilePath=('~/CompartidaWindows/021215/');
FileName=('m021215s01_20min_CA1.brw');
%load([FilePath FileName])
Data=h5read([FilePath FileName],'/3BData/Raw');
Data=single(Data);

%Data structure: Rows=Numbers of channels, colums=time

SelectedChannels=h5read([FilePath FileName],'/3BRecInfo/3BMeaStreams/Raw/Chs');
BitDepth = cast(h5read([FilePath FileName], '/3BRecInfo/3BRecVars/BitDepth'),'uint32');
MaxVolt = h5read([FilePath FileName], '/3BRecInfo/3BRecVars/MaxVolt'); %uVolts
MinVolt = h5read([FilePath FileName], '/3BRecInfo/3BRecVars/MinVolt'); %uVolts
NRecFrames = cast(h5read([FilePath FileName], '/3BRecInfo/3BRecVars/NRecFrames'),'uint32');
SamplingRate = h5read([FilePath FileName], '/3BRecInfo/3BRecVars/SamplingRate');
SignalInversion = h5read([FilePath FileName], '/3BRecInfo/3BRecVars/SignalInversion');


% Converting data to voltage according to MAE/3Brain formula:
% HDFfileslayout.pdf*(pp26)

ADCCountsToMV=SignalInversion.*(MaxVolt-MinVolt)./cast(2^BitDepth,'double');
MVOffset=SignalInversion*MinVolt;
%AnalogValue = MVOffset+DigitalValue*ADCCountsToMV;
Data=MVOffset+Data*ADCCountsToMV;                                                                                                                                                                                                                                                                                                                            
data= cast(Data, 'double');                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
L=length(data);
t=(0:L-1)./SamplingRate;

figure(1)
scatter(SelectedChannels.Col, SelectedChannels.Row, 'sb')
%set(gca,'YDir','Reverse')
title('Selected Channels')
xlim([1 64])
ylim([1 64])


%Removing artefacts

[r c]= size(data);

for i=1:r
    arts=find(data(i,:)>=500 | data(i,:)<=-500);
    data(i,arts)=0;
    %subplot(211); plot(t,data(1,:))
    %subplot(211); plot(t,data(randi([1, r]),:))
    %title('Raw Data')
end

% Apply filter

% LFP filter (<500 Hz)
Nyquist=SamplingRate/2;
LowCutOff=1/Nyquist;
HighCutOff=500/Nyquist;
OrderFilter=5;
BandPass=[LowCutOff HighCutOff];
[Bc Ac]=butter(OrderFilter, BandPass);

     LFP= filtfilt(Bc, Ac, data);
%plot(t,LFP(1,:))
%title('LFP')

% Ripple filter (150-300 Hz)
Nyquist=SamplingRate/2;
LowCutOffR=150/Nyquist;
HighCutOffR=300/Nyquist;
OrderFilter=5;
BandPassR=[LowCutOff HighCutOff];
[BcR AcR]=butter(OrderFilter, BandPassR);

RippleBand= filtfilt(BcR, AcR, LFP);
% plot(t,RippleBand(1,:))
% title('Ripple Band LFP')

% SPW filter (1-60 Hz)
Nyquist=SamplingRate/2;
LowCutOffSPW=1/Nyquist;
HighCutOffSPW=60/Nyquist;
OrderFilter=5;
BandPassSPW=[LowCutOff HighCutOff];
[BcSPW AcSPW]=butter(OrderFilter, BandPassSPW);

SPWBand=filtfilt(BcSPW, AcSPW, LFP);
% plot(t,SPWBand(1,:))
% title('SPW Band LFP')

%downSampling x4 
downsampled=downsample(LFP',4)';
%plot(downsampled)

%wavelet de cada 10 seg
tic
for i=56:90
    i
    figure(i+1)
    % Looking for Ripples (200 hz)
    subplot(323)
    coeffRipp=cwt(downsampled(i,:), 100:300, 'morl', 'plot'); colormap jet;
    title('Wavelet Analysis looking for Ripples')
    % subplot(223)
    % plot(coeffRipp(100,:), 'r-')
    % title('200 hz Coeff')
    % Looking for SPW (50 Hz)
    subplot(324)
    coeffSPW=cwt(downsampled(i,:), 30:60, 'morl', 'plot'); colormap jet;
    title('Wavelet Analysis looking for SPW')
    % subplot(224)
    % plot(coeffSPW(26,:), 'r-')
    % title('50 hz Coeff')


% Now we search index where possible Ripples are in a 200 Hz band (200Hz = coeffRipp(100,:)
% since coeffRipp(1,:)=100 hz and coeffRipp(end,:)=300 Hz, 200Hz = coeffRipp(100,:)
%ie, coeffRipp(1,:)= min-frec of the wavelet and coeffRipp(end,:) = max-frec of the wavelet)

    IndRippDownsampled=find(coeffRipp(100,:)>=4.*std(coeffRipp(100,:)));
    TimesPossibleRipples=4*IndRippDownsampled;
    Window= 1500;
    
    [r1 c1]=size(TimesPossibleRipples);
    IntervalsOfInterestRipp=cell(c1,1);
    for j=1:c1
        region=[];
        down=TimesPossibleRipples(j)-Window;
        up=TimesPossibleRipples(j)+Window;
        if down<=0
            down=1;
            up=2.*Window;
        else
            down=down;
            up=up;
        end    
        region=[down:up];
        IntervalsOfInterestRipp{j,1}=region;
    end

% Now we search index where possible SPW are in a 50 Hz band (50Hz = coeffSPW(26,:)
% since coeffSPW(1,:)=30 hz and coeffSPW(end,:)=60 Hz, 50Hz = coeffSPW(26,:)
%ie, coeffSPW(1,:)= min-frec of the wavelet and coeffSPW(end,:) = max-frec of the wavelet)

    IndSPWDownsampled=find(coeffSPW(26,:)>=4.*std(coeffSPW(26,:)));
    TimesPossibleSPW=4*IndSPWDownsampled;
    %L_Window= 1500;
    %R_Window=1500;
    [r c2]=size(TimesPossibleSPW);
    IntervalsOfInterestSPW=cell(c2,1);
    for j=1:c2
        region=[];
        down=TimesPossibleSPW(j)-Window;
        up=TimesPossibleSPW(j)+Window;
        if down<=0
            down=1;
            up=2.*Window;
        else
            down=down;
            up=up;
        end    
        region=[down:up];
        IntervalsOfInterestSPW{j,1}=region;
    end


    subplot(3,2,[1 2]); plot(t,LFP(i,:))
    title('Possible Ripples(red) ans SPW(green) on LFP')
    hold on
    for h= 1:c1
        plot(t(IntervalsOfInterestRipp{h,1}(1, 1500)), -300:300, '-.r')
        hold on
        %waitbar(i/c)
    end
    for h= 1:c2
        plot(t(IntervalsOfInterestSPW{h,1}(1, 1500)), -200:200, '-.g', 'linewidth', 2)
        hold on
        %waitbar(i/c)
    end
    hold off

        subplot(3,2,5); plot(t,RippleBand(i,:))
        title('Possible Ripples(red)' )
    hold on
    for k= 1:c1
        plot(t(IntervalsOfInterestRipp{k,1}(1, 1500)), -300:300, '-.r')
        hold on
        %waitbar(i/c)
    end
    hold off
        subplot(3,2,6); plot(t,SPWBand(i,:))
        title('Possible SPW(green)')
    hold on
    
    for k= 1:c2
        plot(t(IntervalsOfInterestSPW{k,1}(1, 1500)), -200:200, '-.g', 'linewidth', 2)
        hold on
        %waitbar(i/c)
    end
    hold off
    
    
end
toc