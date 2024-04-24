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

playSound(noisy, Fs);

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
% the above frequencies are 'noisy' so we want to filter them out as much
% as possible (fft output on these peaks should be far lower after
% filtering)
 
%% 2nd Order High-pass Circuit (similar config to fig. 3)
% define circuit constants and initial conditions
R2 = 1.6;       % Resistance in Ohms for R2 -> lower impedence on first circuit (R4 = R2 * 10)
R4 = 16;        % Resistance in Ohms for R4 **CONST**
C1 = 9.947e-5;  % Capacitance in Farads for C1 -> C1 from cutoff frequency eqn: f_c = 100 Hz (sufficiently above 60 Hz grid)
C3 = 9.947e-6;  % Capacitance in Farads for C3 -> higher impedence on second circuit (C3 = C1 / 10)
V0 = [0; 0];    % Initial voltages across C1 and C3
% note: cutoff freq is f_c = 1/(2pi(rc))

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t1, V] = ode45(@(t, V) cascadedRCODE_A(t, noisy, V, R2, R4, C1, C3, round(t / T) + 1), tspan, V0);
V_out = V(:,1) - V(:,2);

% figure 3: output voltage
figure;
plot(t1, V_out);

xlabel("Time (s)");
ylabel("Output Signal Strength (V)");
grid on;

% fourier analysis on output signal to verify captured frequencies
V_outf = abs(fft(V_out));

% figure 4: v_out fft
figure;
hold on;
plot(X_s, V_outf);
hold off;

xlabel("Time (s)");
ylabel("fft output")
grid on;

% figure 5: v_out
figure;
plot(tspan, V_out);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out, Fs);

%% Band-pass filter (similar config to fig. 4)
% Define changed circuit elements
R1 = R2;        % Resistance in Ohms for R1 -> low impedence on first stage (low-pass)
C2 = 3.094e-5;  % Capacitance in Farads for C2 -> low-pass stage: cutoff freq is 3214.89 Hz*

% * from band-pass centre frequency formula: f = sqrt(f_c_high * f_c_low)
%   high-pass stage cutoff freq is same as before at 100 Hz
%   from noisyf_map: most dominant freq that is audible and not 60 Hz (grid)
%   is 567 Hz so we set that as centre frequency for the above

% Solve the system of circuit ODEs
% V(1) = V_1, V(2) = V_C3
[t2, V2] = ode45(@(t, V2) cascadedRCODE_B(t, noisy, V2, R1, R4, C2, C3, round(t / T) + 1), tspan, V0);
V_out2 = V2(:,1) - V2(:,2);

% Fourier analysis on alternate circuit 
V_out2f = abs(fft(V_out2));

% figure 6: v_out2 fft
figure;
hold on;
plot(X_s, V_out2f);
hold off;

xlabel("Time (s)");
ylabel("fft output")
grid on;

% figure 7: v_out2
figure;
plot(t2, V_out2);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out2, Fs);

%% Band-pass filter (similar config to fig. 4) but with 2nd Order High-Pass Input
% This circuit is equivalent to cascading a 2nd order high-pass filter
% described in fig. 3 with ther first order band-pass filter described in
% fig. 4

% Our target isolating frequency ranges are the same here so no circuit
% parameters are changed

[t3, V3] = ode45(@(t, V3) cascadedRCODE_B(t, V_out, V2, R1, R4, C2, C3, round(t / T) + 1), tspan, V0);
V_out3 = V3(:,1) - V3(:,2);

% Fourier analysis   
V_out3f = abs(fft(V_out3));

% figure 8: v_out3 fft
figure;
hold on;
plot(X_s, V_out3f);
hold off;

xlabel("Time (s)");
ylabel("fft output")
grid on;

% figure 9: v_out3
figure;
plot(t3, V_out3);

xlabel("Time (s)");
ylabel("Signal Strength (V)");
grid on;

playSound(V_out3, Fs);


%% Results:
% We would prefer the 2nd order high-pass circuit over the band-pass
% circuit in this case as though the band-pass circuit is theoretically
% more versatile in that there is an upper and lower cutoff frequency, the
% single order setup in this case likely does a worse job at filtering out 
% unwanted frequencies. Furthermore, while both circuits do a good job at
% filtering out the 60 Hz hum, the high-pass circuit actually does a better
% job at filtering out the hissing. This could suggest that the hiss is
% made up of a significant amount of lower frequency elements not covered
% by the low-pass portion of the band-pass filter, thus diminishing the
% usefulness of the band-pass filter's extra functionality. Finally, we
% also experimented on running the 2nd order high-pass circuit's output
% through the first order band-pass circuit. Though this method seems to
% have gotten rid of the hissing present in the other two circuit's
% outputs, the signal making up the music is also heavily diminished and is
% quite hard to make out as its signal is now quite close in strength to
% the reduced hum. Furthermore, this 'cascaded' circuit seems to introduce
% a high pitched noise in lieu of the hissing as well so the audio signal
% outputted here is of lower quality by most metrics compared to the 2nd
% order high-pass circuit by itself.

%% DEs:
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