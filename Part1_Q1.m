clear all;
close all;
% Constants
R = 1000; % Resistance in Ohms
C = 1e-6; % Capacitance in Farads

% Time settings
t_final = 0.01; % Total time in seconds
dt = 0.001; % Time step in seconds
time = 0:dt:t_final; % Time array

% Initialize arrays for current and voltage across the capacitor
I = zeros(size(time));
V_resis = zeros(size(time));

% Compute current and voltage at each time step
for k = 1:length(time) % start from the second index to preserve zero initial condition
    I(k) = 0.005 * exp(-time(k)/(C * R));
    V_resis(k) = R * 0.005 * exp(-time(k)/(C * R));
end

% Plotting the current and voltage across the capacitor
figure;
subplot(2, 1, 1); % Current plot
plot(time, I, 'b');
title('Current through the Resistor');
xlabel('Time (s)');
ylabel('Current (A)');
grid on;

subplot(2, 1, 2); % Voltage across the capacitor plot
plot(time, V_resis, 'r');
title('Voltage across the Capacitor');
xlabel('Time (s)');
ylabel('Voltage (V)');
grid on;
