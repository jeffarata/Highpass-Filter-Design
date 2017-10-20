% Jeff Arata
% 10/19/17

% This project and the associated files were provided by Joe Hoffbeck and
% are found in his paper "Enhance Your DSP Course With These Interesting
% Projects.pdf"

% With this project, we listen to the recording of some birds chirping that
% has some jet engine noise throughout it. We attempt to filter out as much
% of the jet engine noise as possible.

clear;
clc;

[x, fs] = audioread('birds_jet_noise.wav');
soundsc(x, fs);

figure(1)
spectrogram(x, 512, [], 512, fs, 'yaxis');
[S, F, T, P] = spectrogram(x, 512, [], 512, fs, 'yaxis');  % Spectrogram
logP = 10*log10(abs(P));
[PSD, F] = pwelch(x, [], [], [], fs);
logPSD = 10*log10(PSD);                     % Get PSD for plotting

figure(2)           % Plotting is just nice to help see the fundamental 
plot(F, logPSD)     % frequencies and pick a tolerance to find those values.

% Here, the PSD plot doesn't help us very much in finding the frequencies,
% but the spectrogram is very clear. Since the jet engine is mostly noise,
% as opposed to a pitched sound, it is hard to distinguish individual
% frequencies from each other. We see that the jet engine noise stops in
% the spectrogram at around 1500 Hz. In the PSD plot, the dB at this
% frequency is about -102. So bins finds where most frequencies throughout
% all time are less than this value. The first nonzero value in rows is
% then the index of the frequency where the jet engine noise doesn't exist
% for the most part anymore.

bins = (sum(logP, 2) < -102*size(logP, 2));
rows = 1:size(logP, 1);
rows = rows(bins);
cutoff_idx = rows(1);
cutoff = rows(1) * (fs/2) / size(logP,1);

% We'll use that frequency plus a buffer to ensure our filter cuts out as
% much jet engine noise as possible. The value of the buffer in this case
% should be about 150-200Hz, so as not to cut out any bird chirps which the
% spectrogram shows begin at about 1900 Hz.

buffer = 175;
cutoff_plus = cutoff + buffer;

b = fir1(200, cutoff_plus/(fs/2), 'high');    % Coefficients of filter. 'high'
                                              % makes it a highpass filter
            
% Let's see if we can find a more effective filter by another order. Here 
% we will try to minimize the sum of all PSD values from frequencies 0Hz to
% the found cutoff frequency. If a filter of a new order as smaller sum
% here (recall we're working in dB values of the PSD, which are negative),
% then we will accept that as the order of our filter. We keep the order
% within reason, only looping up to an order of 300.

%cutoff_idx = round(size(logP, 1) * cutoff / (fs/2));

db_low_freq = 10;       % Arbitrary initialization as long as it's positive

for jj = 2:2:300        % Increment by 2 for even ordered highpass filters
    b = fir1(jj, cutoff_plus/(fs/2), 'high');
    y = filter(b, 1, x);
    % Optional visualization as the succesive filters change the
    % spectrogram
    %figure(1)
    %spectrogram(y, 512, [], 512, fs, 'yaxis')
    [S, F, T, P] = spectrogram(y, 512, [], 512, fs, 'yaxis');
    logP = 10*log10(abs(P));
    if sum(sum(logP(1:cutoff_idx, :))) < db_low_freq
        order = jj;
        db_low_freq = sum(sum(logP(1:cutoff_idx, :)));
    end
end
                                         
b = fir1(order, cutoff_plus/(fs/2), 'high');
    
figure(3)                       % Plot zeros and poles
zplane(b, 1)

figure(4)                       % Magnitude and Phase of Frequency Response
freqz(b, 1, 512)    

y = filter( b, 1, x );          % Filter the sound
%pause(10)                      % Separates playback if not looping for
                                % most effective filter order.
soundsc(y, fs);                 % Play the sound