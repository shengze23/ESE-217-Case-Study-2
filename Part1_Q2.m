close all;
clear all;

% Initialization
R = 1000;  % Resistance (Ohm)
C = 1e-6;  % Capacitance (F)
omega = 50;  % Frequency (Hz)

initialI = 0;  % Initial current
t_span = [0 0.05];  % Time span
% Use ode45
[t, I] = ode45(@(t, I) rcCircuitODE(t, I, R, C, omega), t_span, initialI);
Fs = length(I)/0.05;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(I); % Length of signal
V_out = (cumtrapz(t, I)) / C;
Y = fft(V_out);
X_s = Fs/L*(0:L-1);
figure;
subplot(2, 1, 1);
plot(t,V_out);
xlabel("time")
ylabel("V_out")
grid on;
subplot(2, 1, 2);
plot(X_s(1:L/2),abs(Y(1:L/2)));
xlabel("f (Hz)")
ylabel("Amplitude for Voltage at frequancy 50")
grid on;


omega_1 = 1000;
[t, I] = ode45(@(t, I) rcCircuitODE(t, I, R, C, omega_1), t_span, initialI);
Fs = length(I)/0.05;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(I);             % Length of signal
V_out_1 = (cumtrapz(t, I)) / C;
Y = fft(V_out_1);
X_s = Fs/L*(0:L-1);
figure;
subplot(2, 1, 1);
plot(t, V_out_1);;
xlabel("time(s)");
ylabel("Out put Voltage");
grid on;
subplot(2, 1, 2);
plot(X_s(1:L/2), abs(Y(1:L/2)));
xlabel("f (Hz)")
ylabel("Amplitude for Voltage at frequancy 1000")
grid on;

function dIdt = rcCircuitODE(t, I, R, C, omega)
    % rcCircuitODE Calculates the derivative of I for a given t, I, R, C, and omega.
    % This function returns the value of the derivative dIdt at the time t.

    % Calculate the time constant
    tau = R * C;

    % Differential equation
    dIdt = (10 * pi * omega / R) * cos(2 * pi * omega * t) - (1 / tau) * I;
end