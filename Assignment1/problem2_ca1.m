% Setting random number generator seed
rng(1234);
% Setting variables

% a 
N = 1024;
%%
% Determining level crossing and fade duration 

nIter = 50;
rholin = [1, 0.1, 0.01];
fdArray = [10, 100];

nvFinal = [];
tauFinal = [];

for i = 1:length(fdArray)
    fd = fdArray(i);
    nvMatrix = [];
    tauMatrix = [];
    
    for j = 1:nIter
        [waveform, sampRate, duration] = ClarkeGansModel(N, fd);
        envelopeMatrix = repmat(waveform, 1, length(rholin));
        % Performing Level crossing detection
        levelDetect = envelopeMatrix > rholin;
        levelCrossShift = [zeros(1, length(rholin)); levelDetect(1:end-1,:)];
        levelCrossDetect = (levelDetect - levelCrossShift) > 0;
        nvIteration = sum(levelCrossDetect, 1);    
        nvMatrix = [nvMatrix; nvIteration];
        
        % Determine the amount of time the signal spends under a given rho
        % value
        levelDetect = envelopeMatrix < rholin;
        nSamplesLevels = sum(levelDetect, 1);
        tauIteration = (1/(sampRate*duration))*nSamplesLevels;
        tauMatrix = [tauMatrix; tauIteration];

    end
    
    % normalization for Nv
    nv = mean(nvMatrix, 1);
    nvNormalized = nv/fd/duration;
    nvFinal = [nvFinal; nvNormalized];
    
    % Normalization for tau_bar
    tau = mean(tauMatrix, 1);
    tau = tau./(nvNormalized);
    tauFinal = [tauFinal; tau]
end
% computing the theoretical values
%%
nvtheoretical = sqrt(2*pi)*rholin.*exp(-rholin.^2);
tauideal = (exp(rholin.^2) - 1)./(sqrt(2*pi)*rholin);
disp("fd values");
disp(fdArray');
disp("Simulated results");
disp(nvFinal);
disp("Theoretical values")
disp(nvtheoretical);
disp("fd values");
disp(fdArray');
disp("Simulated results");
disp(tauFinal);
disp("Theoretical values");
disp(tauideal);
%% Functions

%{
sezArray = SEZ(f,fd)

Inputs:
    f (array): array of freqeuencies at which the slope of the spectrum
    must be computed, all in Hz
    fd (float): maximum doppler frequency shift in Hz
Outputs:
    sezArray (array): array of values of the given SEZ function at the
    specified frequencies
%}
%%
function sezArray = SEZ(f,fd)
    sezArray = 1.5./(pi*fd*sqrt(1 - (f/fd).^2));
    % Dealing with the infinity values at the ends of the filter
    
    slopePrev = SEZslope(f(end-1), fd);
    newValue = sezArray(end - 1) + slopePrev*(f(end) - f(end - 1));
    sezArray(end) = newValue;
    sezArray(1) = newValue;
end

%{
slope = SEZslope(f,fd)

Inputs:
    f (array): array of freqeuencies at which the slope of the spectrum
    must be computed, all in Hz
    fd (float): maximum doppler frequency shift in Hz
Outputs:
    slope (array): array of slopes of the given SEZ function at the
    specified frequencies
%}
function slope = SEZslope(f,fd)
    slope = (1.5.*f.*(1 - (f/fd).^2).^(-3/2))/(pi*fd^3);
    
end


%{
function waveformRayleigh = ClarkGansModel(N,fd)

Inputs 
    N (int): Number of points in the line spectrum
    fd (float): maximum doppler frequency shift in Hz
Outputs
    waveformRayleigh (array): array of size N x 1 withe time domain values
        of the rayleigh fading waveform
    sampRate (float): sampling rate of the fading waveform
    T (float): duration of the waveform in s
%}
function [waveformRayleigh, sampRate, T] = ClarkeGansModel(N, fd)
    
    % b
    deltaF = 2*fd/(N - 1);
    T = 1/deltaF;
    K = 256*N;
    Ts = T/K;
    sampRate = 1/Ts;
    
    
    % c
    mu = 0;
    sigma = 1/sqrt(2);
    
    % generating noise source 1
    inPhasePosFreq = normrnd(mu, sigma, [1, N/2]);
    quadraturePosFreq = normrnd(mu, sigma, [1, N/2]);
%     posFreqValues = inPhasePosFreq + 1j*quadraturePosFreq;
%     
    % d
%     negFreqValues = inPhasePosFreq - 1j*quadraturePosFreq;
    % the higher the index in negFreqValues, the more negative the frequency
    inPhase = [flip(inPhasePosFreq), inPhasePosFreq];
    quadrature = [-1*flip(quadraturePosFreq), quadraturePosFreq];
    noiseSource1 = inPhase + 1j*quadrature;

    
    % generating noise source 2
    inPhasePosFreq = normrnd(mu, sigma, [1, N/2]);
    quadraturePosFreq = normrnd(mu, sigma, [1, N/2]);
%     posFreqValues = inPhasePosFreq + 1j*quadraturePosFreq;
%     
    % d
%     negFreqValues = inPhasePosFreq - 1j*quadraturePosFreq;
    % the higher the index in negFreqValues, the more negative the frequency
    inPhase = [flip(inPhasePosFreq), inPhasePosFreq];
    quadrature = [-1*flip(quadraturePosFreq), quadraturePosFreq];
    noiseSource2 = inPhase + 1j*quadrature;

    
    % e 
    freqs = -fd:deltaF:fd;
    
    filter = sqrt(SEZ(freqs, fd));
    
    
    noiseSourceFiltered1 = filter.*noiseSource1;
    noiseSourceFiltered2 = filter.*noiseSource2;
    

    %padding
    ns1 = [zeros(1, K/2 - N/2), noiseSourceFiltered1, zeros(1, K/2- N/2)];
    ns2 = [zeros(1, K/2 - N/2), noiseSourceFiltered2, zeros(1, K/2- N/2)];
    
    % f
    inPhaseShifted = ifftshift(ns1); %circshift(noiseSourceFiltered1, N/2);
    quadratureShifted = ifftshift(ns2); %circshift(noiseSourceFiltered2, N/2);
    
    
    
    inPhaseTD = ifft( [inPhaseShifted]);
    quadratureTD = ifft([ quadratureShifted]);
    
    % g
    
    signalTD = sqrt(abs(inPhaseTD).^2 + abs(quadratureTD).^2);
    signalNorm = signalTD/(sqrt(mean(signalTD.^2))); % normalize envelope to 0dB
    waveformRayleigh = signalNorm';
end