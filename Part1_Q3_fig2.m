clear;
close all;
% Constants and initial conditions
R1 = 1000; % Resistance in Ohms for R2
R4 = 1000; % Resistance in Ohms for R4
C2 = 1e-6; % Capacitance in Farads for C1
C3 = 1e-6; % Capacitance in Farads for C3
a1 = 1; % Amplitude for the first sinusoidal component
a2 = 1; % Amplitude for the second sinusoidal component
Vin = @(t)(a1 * sin(2 * pi * 50 * t) + a2 * sin(2 * pi * 10^5 * t));
initialV = [0; 0]; % Initial voltages across C1 and C3
% Time span for the simulation
t_span = [0, 0.1]; % 0 to 0.1 seconds
% Solve the system of ODEs
[t, V] = ode45(@(t, V) cascadedRCODE(t, V, R1, R4, C2, C3, Vin), t_span, initialV);
V_out = V(:,1) - V(:,2);

Fs = length(V_out)/0.1;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(V_out);             % Length of signal
Y = fft(V_out);
X_s = Fs/L*(0:L-1);

figure;
plot(t, V(:,1), 'b', t, V(:,2), 'r', t, V_out, 'g');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend('V_{C2}', 'V_{C3}', 'v_{out} ','Location', 'best');
title('Response of Cascaded RC Circuits to Compound Signal');
grid on;

figure;
subplot(2, 1, 1);
plot(t, V_out);;
xlabel("time(s)");
ylabel("Out put Voltage");
grid on;
subplot(2, 1, 2);
plot(X_s, abs(Y));
xlabel("f (Hz)")
ylabel("fft graph")
grid on;

function dVdt = cascadedRCODE(t, V, R1, R4, C2, C3, Vin)
    % Input voltage function handle
    Vin = Vin(t); % Evaluate Vin at time t
    
    % Voltage across C2 and C3
    V_C2 = V(1); 
    V_C3 = V(2); 
    
    % Current through R1
    i_R1 = (Vin - V_C2) / R1; 
    
    % Current through R4 is the same as current through C3 due to series connection
    i_R4 = (V_C2 - V_C3) / R4;
    
    % Differential equations for the voltages across capacitors C2 and C3
    dV_C2dt = i_R1 / C2; % Charging of C2 by current i_R1
    dV_C3dt = i_R4 / C3; % Charging of C3 by current i_R4 (also i3)
    
    % Output derivative vector
    dVdt = [dV_C2dt; dV_C3dt];
end