%% ESE 217 Case Study 2 Part 2
% Chuan Shuo Chen and Simon He

clear;
close all;

S = load("noisyhandel.mat");    % loads data

% setup audio signal and metdata
noisy = S.Vsound;
Fs    = S.Fs;  % 44100
T     = 1/Fs;
L     = length(noisy);

%playSound(noisy, Fs);

% figure 1: noisyhandle soundwaves
time = 0:T:(L - 1) * T;
plot(time, noisy);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

% fourier analysis on original signal to isolate 'noisy' frequency
noisyf = abs(fft(noisy));   % only care about absolute value for amp.
X_s    = Fs/L*(0:L-1);      % setup x-axis

% figure 2: noisyhandel fft
figure;
plot(X_s, noisyf);

xlabel("Frequency (Hz)");
ylabel("fft output");
grid on;

% sort noisyf data in descending order
noisyf_map = [noisyf, X_s.'];                   % each amplitude is assigned to its frequency
noisyf_map = sortrows(noisyf_map, 'descend');   % key-value pairings preserved

% dominant frequencies from fft are: 60 Hz, 44040 Hz
% design a circuit to isolate the above frequencies ie.
% vout = a1sin(2pi60t) + a2sin(2pi44040t)

% define circuit constants and initial conditions
R2 = 1000;      % Resistance in Ohms for R2
R4 = 16;        % Resistance in Ohms for R4 **CONST**
C1 = 1e-6;      % Capacitance in Farads for C1
C3 = 1e-6;      % Capacitance in Farads for C3
V0 = [0; 0];    % Initial voltages across C1 and C3

tspan = [0, max(time)];

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t, V] = ode45(@(t, V) cascadedRCODE(t, noisy, V, R2, R4, C1, C3), tspan, V0);
V_out = V(:,1) - V(:,2);

function dVdt = cascadedRCODE(~, Vin, V, R2, R4, C1, C3)
    dV1dt = (Vin - V(1)) / (R2*C1); % dV1/dt for C1
    dV3dt = (V(1) - V(2)) / (R4*C3); % dV3/dt for C3
    dVdt = [dV1dt; dV3dt]; % return derivatives as a column vector
end

%% function playSound(y, Fs)
% scales the audio signal in y to range from -1.0 to 1.0, then plays the
% entire sequence. Control is returned once the sequence is finished
% playing. 
%
% y - audio time sequence 
% Fs - sampling rate for y (Hz)
%
% Matthew Lew 10/25/2018

function playSound(y, Fs)
obj = audioplayer(y/max(abs(y)), Fs);
playblocking(obj);
end