%% ESE 217 Case Study 2 Part 2
% Chuan Shuo Chen and Simon He

clear;
close all;

S = load("noisyhandel.mat");    % loads data

%% Setup audio signal source and metadata
noisy = S.Vsound;
Fs    = S.Fs;  % 44100
T     = 1/Fs;
L     = length(noisy);

tspan  = 0:T:(L - 1) * T;
dnoisy = gradient(noisy(:)) ./ gradient(tspan(:));

%playSound(noisy, Fs);

% figure 1: noisyhandle soundwaves
tspan = 0:T:(L - 1) * T;
plot(tspan, noisy);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

% fourier analysis on original signal to isolate 'noisy' frequency
noisyf = abs(fft(noisy));   % only care about absolute value for amp.
X_s    = Fs/L*(0:L-1);      % setup x-axis

% figure 2: noisyhandel fft
figure;
plot(X_s(1:L/2), noisyf(1:L/2));

xlabel("Frequency (Hz)");
ylabel("fft output");
grid on;

% sort noisyf data in descending order
noisyf_map = [noisyf, X_s.'];                   % each amplitude is assigned to its frequency
noisyf_map = sortrows(noisyf_map, 'descend');   % key-value pairings preserved

% dominant frequencies from fft are: 60 Hz (likely from electrical grid)
% the above frequencies are 'noisy' so we want to filter them out as much
% as possible (fft output on these peaks should be far lower after
% filtering)
 
%% 2nd Order High-pass Circuit (similar config to fig. 3)
% define circuit constants and initial conditions
R2 = 1.6;       % Resistance in Ohms for R2 -> lower impedence on first circuit (R4 = R2 * 10)
R4 = 16;        % Resistance in Ohms for R4 **CONST**
C1 = 9.947e-5;  % Capacitance in Farads for C1 -> C1 from cutoff frequency eqn: f_c = 100 Hz (sufficiently above grid freq)
C3 = 9.947e-6;  % Capacitance in Farads for C3 -> higher impedence on second circuit (C3 = C1 / 10)
V0 = [0; 0];    % Initial voltages across C1 and C3
% note: cutoff freq is f_c = 1/(2pi(rc))

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t1, V] = ode45(@(t, V) cascadedRCODE_A(t, dnoisy, V, R2, R4, C1, C3, round(t / T) + 1), tspan, V0);
V_out = V(:,2);

% fourier analysis on output signal to verify captured frequencies
V_outf = abs(fft(V_out));

% figure 3: v_out fft
figure;
hold on;
plot(X_s(1:L/2), V_outf(1:L/2));
hold off;

xlabel("Frequency (Hz)");
ylabel("fft output")
grid on;

% figure 4: v_out
figure;
plot(tspan, V_out);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out, Fs);

%% Band-pass filter (similar config to fig. 4)
% Define changed circuit elements
R1 = R2;        % Resistance in Ohms for R1 -> low impedence on first stage (low-pass)
C2 = 1.856e-5;  % Capacitance in Farads for C2 -> low-pass stage: cutoff freq is 5358.15 Hz*
C3 = C1;        % Update C3 value based on C1 since this is now a single stage high-pass

% * from band-pass centre frequency formula: f = sqrt(f_c_high * f_c_low)
%   high-pass stage cutoff freq is same as before at 100 Hz
%   from noisyf_map: most dominant freq that is audible and not 60 Hz (grid)
%   is 567 Hz so we set that as centre frequency for the above

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t2, V2] = ode45(@(t, V2) cascadedRCODE_B(t, noisy, V2, R1, R4, C2, C3, round(t / T) + 1), tspan, V0);
V_out2 = V2(:,2);

% Fourier analysis on alternate circuit 
V_out2f = abs(fft(V_out2));

% figure 5: v_out2 fft
figure;
hold on;
plot(X_s(1:L/2), V_out2f(1:L/2));
hold off;

xlabel("Frequency (Hz)");
ylabel("fft output")
grid on;

% figure 6: v_out2
figure;
plot(t2, V_out2);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out2, Fs);

%% Band-Pass Filter with High-Pass Filter input

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t3, V3] = ode45(@(t, V3) cascadedRCODE_B(t, V_out, V2, R1, R4, C2, C3, round(t / T) + 1), tspan, V0);
V_out3 = V3(:,2);

% Fourier analysis on alternate circuit 
V_out3f = abs(fft(V_out3));

% figure 5: v_out2 fft
figure;
hold on;
plot(X_s(1:L/2), V_out3f(1:L/2));
hold off;

xlabel("Frequency (Hz)");
ylabel("fft output")
grid on;

ylim([0, 3e4])

% figure 6: v_out2
figure;
plot(t3, V_out3);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out3, Fs);

%% Results:
% The high-pass filter does a much better job at filtering out the 60 Hz
% hum from the electrical grid than the band-pass filter. This can be
% clearly seen by comparing the two fft plots for the two filters, with a
% noticeable spike at 60 Hz still appearing in the band-pass filter output
% fft but no noticeable spike on the high-pass filter fft. This is likely
% because the second-order nature of the high-pass filter is able to much
% more effectively remove the 60 Hz sound waves compared to the band-pass
% filter's high-pass stage, which is only a first-order circuit. However,
% the band-pass filter does a much better job at removing higher frequency
% hissing as well as it has a lower and upper cutoff frequency as opposed
% to the high-pass filter's single lower cutoff. This also means that the
% clarity of the music is better preserved in the band-pass filter.
%
% By passing the output of the high-pass filter through the band-pass
% filter, we do actually get a signal that has the hum mostly filtered out
% and music with a comparatively high degree of clarity. However, this
% output signal is far more quiet compared to either of the previous filter
% outputs and also has a very noticeable 'crackling' noise is the
% background.

%% DEs:
function dVdt = cascadedRCODE_A(~, dVin, V, R2, R4, C1, C3, counter)
    dV1dt = ((R2 * C1 * dVin(counter)) - V(1))/(R2 * C1);
    dV3dt = ((R4 * C3 * dV1dt) - V(2))/(R4 * C3);
    dVdt  = [dV1dt; dV3dt];                     % return derivatives as a column vector
end

function dVdt = cascadedRCODE_B(~, Vin, V, R1, R4, C2, C3, counter)
    dV1dt = ((Vin(counter) - V(1)))/(R1 * C2);
    dV3dt = ((R4 * C3 * dV1dt) - V(2))/(R4 * C3);
    dVdt  = [dV1dt; dV3dt];
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