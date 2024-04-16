close all;
clear all;

% they are in series
R = 1000; % Resistance (Ohm)
C = 1e-6; % Capacitance (F)
V = 5; % Constant voltage (V_tot)

% set the time
t_final = 0.01; % Total time
dt = 1e-5; % Time step
time = 0:dt:t_final; % array of time

% Initialization
Vc = zeros(size(time)); % Voltage across the capacitor
I = zeros(size(time)); % Current through the circuit

% Initially the voltage across capaciator will be 0v
Vc(1) = 0; % set the first element in array be zero

% Using eruler's method
for k = 1:length(time)-1
    I(k) = (V - Vc(k)) / R; % Ohm's law: V = IR, so I = (V - Vc)/R
    dVc = I(k) / C; % Change in capacitor voltage
    Vc(k+1) = Vc(k) + dVc * dt; % Update capacitor voltage
end

% Calculate current for the last step
I(end) = (V - Vc(end)) / R;

% Plot
figure;
subplot(2,1,1);
plot(time, Vc);
title('Capacitor Voltage (Vc) vs. Time');
xlabel('Time (s)');
ylabel('Voltage (V)');

subplot(2,1,2);
plot(time, I);
title('Circuit Current (I) vs. Time');
xlabel('Time (s)');
ylabel('Current (A)');


% plot the equation solving by hand
I = (V / R) * exp(-time / (R * C)); % Current through the resistor
Vc = V * (1 - exp(-time / (R * C))); % Voltage across the capacitor

% Plotting the results
figure;
subplot(2,1,2); % Subplot for current
plot(time, I, 'b', 'LineWidth', 2);
title('Circuit Current I(t)');
xlabel('Time (s)');
ylabel('Current (A)');
grid on;

subplot(2,1,1); % Subplot for voltage
plot(time, Vc, 'r', 'LineWidth', 2);
title('Capacitor Voltage V_C(t)');
xlabel('Time (s)');
ylabel('Voltage (V)');
grid on;


