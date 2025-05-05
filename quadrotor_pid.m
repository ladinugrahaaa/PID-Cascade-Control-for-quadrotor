clear; clc;

%% MODEL PARAMETERS
g = 9.81;        
m = 1.343;       
Ix = 0.01864;    
Iy = 0.037;      
Iz = 0.00554;  
kt = 0.1939; 
Jr = 0.007588; 
Omega = 1;       

%% TIME
dt = 0.01; 
t = 0:dt:100; 
steps = length(t);

%% STATE INITIALIZATION
X = zeros(12, steps);
X(:,1) = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

%% DESIRED TRAJECTORY
length_traj = 20;
width      = 20;
z_height   = 10;
side_time  = 20;
total_time = 4 * side_time;

xsp = zeros(3, steps);

for k = 1:steps
    current_time = t(k);
    if current_time <= side_time
        xsp(:, k) = [length_traj * current_time / side_time; 0; z_height];
    elseif current_time <= 2 * side_time
        xsp(:, k) = [length_traj; width * (current_time - side_time) / side_time; z_height];
    elseif current_time <= 3 * side_time
        xsp(:, k) = [length_traj - length_traj * (current_time - 2 * side_time) / side_time; width; z_height];
    elseif current_time <= 4 * side_time
        xsp(:, k) = [0; width - width * (current_time - 3 * side_time) / side_time; z_height];
    else
        xsp(:, k) = [0; 0; z_height];
    end
end

%% PID GAIN PARAMETERS
% Outer Loop
Kp_x = 6.2590;  Ki_x = 9.9998;  Kd_x = 7.8902;
Kp_y = 5.3639;  Ki_y = 5.5670;  Kd_y = 7.3942;
Kp_z = 10.0000;  Ki_z = 0.0001;  Kd_z = 3.4794;

% Inner Loop
Kp_phi = 8;  Ki_phi = 0.0;  Kd_phi = 3.0;
Kp_theta = 8; Ki_theta = 0.0; Kd_theta = 3.0;
Kp_psi = 4;   Ki_psi = 0.0;  Kd_psi = 1.0;

%% error integral (outer & inner loop)
e_int_x = 0;  e_int_y = 0;  e_int_z = 0;
e_int_phi = 0;  e_int_theta = 0;  e_int_psi = 0;

%% Array kontrol input
U = zeros(4, steps);

%% Simulation Loop
figure(1);
for k = 1:steps-1
    %% Outer Loop
    x = X(1,k);   xdot = X(2,k);
    y = X(3,k);   ydot = X(4,k);
    z = X(5,k);   zdot = X(6,k);
    
    x_des = xsp(1,k);
    y_des = xsp(2,k);
    z_des = xsp(3,k);
    
    error_x = x_des - x;
    error_y = y_des - y;
    error_z = z_des - z;
    
    error_dot_x = - xdot;
    error_dot_y = - ydot;
    error_dot_z = - zdot;
    
    e_int_x = e_int_x + error_x * dt;
    e_int_y = e_int_y + error_y * dt;
    e_int_z = e_int_z + error_z * dt;
    
    a_x_des = Kp_x * error_x + Ki_x * e_int_x + Kd_x * error_dot_x;
    a_y_des = Kp_y * error_y + Ki_y * e_int_y + Kd_y * error_dot_y;
    a_z_des = Kp_z * error_z + Ki_z * e_int_z + Kd_z * error_dot_z;
    
    theta_des = - a_x_des / g;
    phi_des   =   a_y_des / g;
    psi_des = 0;
    
    U1 = m * (g - a_z_des);
    
    %% Inner Loop
    phi = X(7,k);    phidot = X(8,k);
    theta = X(9,k);  thetadot = X(10,k);
    psi = X(11,k);   psidot = X(12,k);
    
    error_phi = phi_des - phi;
    error_theta = theta_des - theta;
    error_psi = psi_des - psi;

    error_dot_phi = - phidot;
    error_dot_theta = - thetadot;
    error_dot_psi = - psidot;

    e_int_phi = e_int_phi + error_phi * dt;
    e_int_theta = e_int_theta + error_theta * dt;
    e_int_psi = e_int_psi + error_psi * dt;

    U2 = Kp_phi   * error_phi   + Ki_phi   * e_int_phi   + Kd_phi   * error_dot_phi;
    U3 = Kp_theta * error_theta + Ki_theta * e_int_theta + Kd_theta * error_dot_theta;
    U4 = Kp_psi   * error_psi   + Ki_psi   * e_int_psi   + Kd_psi   * error_dot_psi;

    U(:,k) = [U1; U2; U3; U4];
    
    %% Update state menggunakan model dinamika nonlinear
    Xdot = computeStateDerivative(X(:,k), U(:,k), g, m, Ix, Iy, Iz, kt, Jr);
    X(:,k+1) = X(:,k) + dt * Xdot;
    
    %% Plotting Trajectory setiap beberapa iterasi
    if mod(k, 10) == 0
        fprintf('Step %d: x=%.2f, y=%.2f, z=%.2f, phi=%.2f°, theta=%.2f°, psi=%.2f°\n', ...
            k, X(1,k), X(3,k), X(5,k), rad2deg(X(7,k)), rad2deg(X(9,k)), rad2deg(X(11,k)));
        clf;
        hold on;
        % Plot trajectory referensi dan lintasan aktual
        plot3(xsp(1,1:k), xsp(2,1:k), xsp(3,1:k), 'k--', 'LineWidth', 1.5);
        plot3(X(1,1:k), X(3,1:k), X(5,1:k), 'b', 'LineWidth', 1.5);
    
        % Gambar quadrotor dengan bentuk yang lebih nyata
        L = 1; % panjang setengah lengan quadrotor (sesuaikan dengan skala)
        drawQuadrotor(X(:,k), L);
    
        % Plot titik-titik sudut lintasan (misalnya, jika diperlukan)
        scatter3([0, length_traj, length_traj, 0], [0, 0, width, width], ...
             [z_height, z_height, z_height, z_height], 50, 'm', 'filled', 's');
        axis equal;
        title('3D Trajectory of the Drone with Quadrotor Shape');
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        legend('Reference Trajectory', 'Actual Trajectory', 'Quadrotor');
        grid on;
        view(3);
        drawnow;
    end
    if k * dt >= total_time
        break;
    end
end

%% PLOT STATES DAN CONTROL INPUTS
figure(2);
for i = 1:12
    subplot(3,4,i);
    plot(t(1:k), X(i,1:k), 'LineWidth', 1.5);
    title(['State ', num2str(i)]);
    xlabel('Time [s]');
    ylabel(['X', num2str(i)]);
    grid on;
end

figure(3);
for i = 1:4
    subplot(2,2,i);
    plot(t(1:k), U(i,1:k), 'LineWidth', 1.5);
    title(['Control Input U', num2str(i)]);
    xlabel('Time [s]');
    ylabel(['U', num2str(i)]);
    grid on;
end

%% NONLINEAR STATE
function Xdot = computeStateDerivative(X, U, g, m, Ix, Iy, Iz, kt, Jr)
    xdot   = X(2);
    ydot   = X(4);
    zdot   = X(6);
    phi    = X(7);
    phidot = X(8);
    theta  = X(9);
    thetadot = X(10);
    psi    = X(11);
    psidot = X(12);
    
    c_phi = cos(phi);   s_phi = sin(phi);
    c_theta = cos(theta); s_theta = sin(theta);
    c_psi = cos(psi);   s_psi = sin(psi);
    
    Xdot = zeros(12,1);
    Xdot(1) = xdot;
    Xdot(2) = -U(1)/m * (c_phi*s_theta*c_psi + s_phi*s_psi) + (kt/m)*xdot;
    Xdot(3) = ydot;
    Xdot(4) = -U(1)/m * (c_phi*s_theta*s_psi - s_phi*c_psi) + (kt/m)*ydot;
    Xdot(5) = zdot;
    Xdot(6) = g - U(1)/m * (c_phi*c_theta) + (kt/m)*zdot;
    Xdot(7) = phidot;
    Xdot(8) = (0.225/Ix)*U(2) + (Iy - Iz)/Ix * thetadot*psidot - kt/Ix*phidot - Jr/Ix*thetadot;
    Xdot(9) = thetadot;
    Xdot(10) = (0.225/Iy)*U(3) + (Iz - Ix)/Iy * phidot*psidot - kt/Iy*thetadot - Jr/Iy*phidot;
    Xdot(11) = psidot;
    Xdot(12) = (0.225/Iz)*U(4) + (Ix - Iy)/Iz * phidot*thetadot - kt/Iz*psidot;
end

%% Quadrotor
function drawQuadrotor(state, L)
    % Fungsi untuk menggambar quadrotor berdasarkan state
    % state: [x; xdot; y; ydot; z; zdot; phi; phidot; theta; thetadot; psi; psidot]
    % L: panjang setengah lengan quadrotor (misal, L = 1)
    
    % Ekstrak posisi dan orientasi (ambil yaw, pitch, roll)
    x = state(1);
    y = state(3);
    z = state(5);
    phi   = state(7);  % roll
    theta = state(9);  % pitch
    psi   = state(11); % yaw
    
    % Buat matriks rotasi menggunakan konvensi ZYX
    % MATLAB menyediakan fungsi eul2rotm dengan format [yaw, pitch, roll]
    R = eul2rotm([psi, theta, phi]);
    
    % Definisikan titik-titik ujung lengan di kerangka badan (body frame)
    % Misalnya, lengan horizontal (sumbu x) dan lengan vertikal (sumbu y)
    p1 = [ L; 0; 0];   % ujung kanan
    p2 = [-L; 0; 0];   % ujung kiri
    p3 = [0;  L; 0];   % ujung depan
    p4 = [0; -L; 0];   % ujung belakang
    
    % Transformasi titik-titik tersebut ke kerangka inersial
    center = [x; y; z];
    P1 = center + R * p1;
    P2 = center + R * p2;
    P3 = center + R * p3;
    P4 = center + R * p4;
    
    % Gambar lengan quadrotor sebagai garis
    plot3([P1(1), P2(1)], [P1(2), P2(2)], [P1(3), P2(3)], 'r-', 'LineWidth', 2);
    plot3([P3(1), P4(1)], [P3(2), P4(2)], [P3(3), P4(3)], 'r-', 'LineWidth', 2);
    
    % Gambar pusat quadrotor
    plot3(x, y, z, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
end
