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
plot(X_s, noisyf);

xlabel("Frequency (Hz)");
ylabel("fft output");
grid on;

% sort noisyf data in descending order
noisyf_map = [noisyf, X_s.'];                   % each amplitude is assigned to its frequency
noisyf_map = sortrows(noisyf_map, 'descend');   % key-value pairings preserved

% dominant frequencies from fft are: 60 Hz (likely from electrical grid), 44040 Hz
% design a circuit to isolate the above frequencies ie.

% define circuit constants and initial conditions
R2 = 16;        % Resistance in Ohms for R2
R4 = 16;        % Resistance in Ohms for R4 **CONST**
C1 = 1.658e-4;  % Capacitance in Farads for C1
C3 = 1.658e-4;  % Capacitance in Farads for C3 -> C1 and C3 from cutoff frequency eqn where f_c = 60 Hz (grid freq)
V0 = [0; 0];    % Initial voltages across C1 and C3
% note: cutoff freq is f_c = 1/(2pi(rc))

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t, V] = ode45(@(t, V) cascadedRCODE_A(t, noisy, V, R2, R4, C1, C3, round(t / T) + 1), tspan, V0);
V_out = V(:,1) - V(:,2);

% figure 3: output voltage
figure;
plot(tspan, V_out);

xlabel("Time (s)");
ylabel("Output Signal Strength (V)");
grid on

% fourier analysis on output signal to verify captured frequencies
V_outf = abs(fft(V_out));

% figure 4: v_out fft
figure;
hold on;
xline(44040, '--magenta');  % target freq
plot(X_s, V_outf);
hold off;

xlabel("Time (s)");
ylabel("fft output")

% figure 5: v_out
figure;
plot(tspan, V_out);

xlabel("Time (s)");
ylabel("Signal Strength (V)");

playSound(V_out, Fs);

% Alternate circuit design
% Define changed circuit elements
R1 = 16;
C2 = 1.658e-4;

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t, V2] = ode45(@(t, V2) cascadedRCODE_B(t, noisy, V2, R1, R4, C2, C3, round(t / T) + 1), tspan, V0);
V_out2 = V2(:,1) - V2(:,2);

% Fourier analysis on alternate circuit 
V_out2f = abs(fft(V_out2));

% figure 6: v_out2 fft
figure;
hold on;
xline(44040, '--magenta');  % target freq
plot(X_s, V_out2f);
hold off;

xlabel("Time (s)");
ylabel("fft output")

% figure 7: v_out2
figure;
plot(tspan, V_out2);

xlabel("Time (s)");
ylabel("Signal Strength (V)");

playSound(V_out2, Fs);

% We would prefer the original circuit (modeled by cascadedRCODE_A) since
% it seems to have a better signal strength, both circuits seem to filter
% out similar frequencies as shown in their respective ffts so quality is
% approximately similar.

function dVdt = cascadedRCODE_A(~, Vin, V, R2, R4, C1, C3, counter)
    dV1dt = (Vin(counter) - V(1)) / (R2 * C1);
    dV3dt = (V(1) - V(2)) / (R4 * C3);
    dVdt  = [dV1dt; dV3dt];                     % return derivatives as a column vector
end

function dVdt = cascadedRCODE_B(~, Vin, V, R1, R4, C2, C3, counter)
    dV1dt = (Vin(counter) - V(1) / (R1 * C2));
    dV3dt = (V(1) - V(2)) / (R4 * C3);
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