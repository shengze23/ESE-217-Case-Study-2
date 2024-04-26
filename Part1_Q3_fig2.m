clear;
close all;

% Constants and initial conditions
R1 = 1000; % Resistance in Ohms for R2
R4 = 1000; % Resistance in Ohms for R4
C2 = 1e-6; % Capacitance in Farads for C1
C3 = 1e-6; % Capacitance in Farads for C3
a1 = 1; % Amplitude for the first sinusoidal component
a2 = 1; % Amplitude for the second sinusoidal component
initialV = [0; 0]; % Initial voltages across C1 and C3

% Time span for the simulation
t_span = [0, 0.1]; % 0 to 0.1 seconds

% Solve the system of ODEs
[t, V] = ode45(@(t, V) cascadedRCODE(t, V, R1, R4, C2, C3, a1, a2), t_span, initialV);
V_out = V(:,2);

Fs = length(V_out)/0.1;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(V_out);             % Length of signal
Y = fft(V_out);
X_s = Fs/L*(0:L-1);

% Plot the results
figure;
plot(t, V(:,1), 'b', t, V(:,2), 'r', t, V_out, 'g');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_{C1}', 'V_{C3}','V_{out}', 'Location', 'best');
title('Response of Cascaded RC Circuits to Compound Signal');
grid on;

figure;
subplot(2, 1, 1);
plot(t, V_out);
xlabel("time(s)");
ylabel("Out put Voltage");
grid on;
subplot(2, 1, 2);
plot(X_s(1:L/2), abs(Y(1:L/2)));
xlabel("f (Hz)")
ylabel("fft graph")
grid on;

function dVdt = cascadedRCODE(t, V, R1, R4, C2, C3, a1, a2)
    V_in = a1*sin(2*pi*50*t) + a2*sin(2*pi*(10^5)*t);% compound input signal
    V_in_der = 100 * pi * a1 * cos(2*pi*50*t) + 2*pi*(10^5) * a2 * cos(2*pi*(10^5)*t);
    dV1dt = (V_in - V(1))/(R1 * C2);
    dV3dt = ((R4 * C3 * dV1dt) - V(2))/(R4 * C3); % dV3/dt for C3
    dVdt = [dV1dt; dV3dt]; % return derivatives as a column vector
end
