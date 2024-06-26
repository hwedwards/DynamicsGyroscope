close all; 
clear all; 
clc; 
% Symbolic Variables:
syms t alpha(t) beta(t) gamma(t)
syms h d1 d2 L L2
syms t1 t2
syms r1 r2 r3 r4
syms m1 m2 m3 m4

% L - length of gyroscope axle
% L2 - length of rotating axle
% d1 - length from A to start of big disc (B)
% d2 - length from A to start of counterweight
% t1 - thickness of disc
% t2 - thickness of counterweight
% m1 - mass of large disc
% m2 - mass of counterweight
% m3 - mass of axle
% m4 - mass of rotating axle
% r1 - radius of big disc
% r2 - radius of small disc
% r3 - radius of gyroscope axle
% r4 - radius of rotating axle

% Positive rotation around z axis.
R01 = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
R10 = R01';

% Negative rotation around y axis.
R12 = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
R21 = R12.';

% Positive rotation around x axis.
R23 = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
R32 = R23.';

R20 = R21 * R10;
R02 = R20.'; 
R30 = R32 * R20;
R03 = R30.'; 

w1_1 = [0; 0; diff(alpha, t)];

w21_2 = [0; -diff(beta, t); 0];
w2_2 = w21_2 + R21*w1_1;

w32_3 = [diff(gamma, t); 0; 0];
w3_3 = w32_3 + R32*w21_2 + R32*R21*w1_1;

rOC3_1 = [0; 0; L2/2]; % Base to COM of rotating axle
rAO_1 = [0; 0; -L2/2]; % A to COM of rotating axle

rACaxle_2 = [-(L/2 - d1); 0; 0]; % A to COM of gyroscope axle
rAD1_2 = [d1 + t1/2; 0; 0]; % A to COM of big disc
rAD2_2 = [-(d2 + t2/2); 0; 0]; % A to COM of counterweight
rAC2_2 = (m3*rACaxle_2 + m2*rAD2_2)/(m3 + m2); % combined COM of axle + counterweight

rOC3_1_dot = diff(rOC3_1, t) + cross(w1_1, rOC3_1);
rOC3_1_ddot = diff(rOC3_1_dot, t) + cross(w1_1, rOC3_1_dot);

rAC2_2_dot = diff(rAC2_2, t) + cross(w2_2, rAC2_2);
rAC2_2_ddot = diff(rAC2_2_dot, t) + cross(w2_2, rAC2_2_dot);

rAD1_2_dot = diff(rAD1_2, t) + cross(w2_2, rAD1_2); 
rAD1_2_ddot = diff(rAD1_2_dot, t) + cross(w2_2, rAD1_2_dot);

%%%% Rotating Rod stuff
% Inertia tensor of axle about its own COM
xx = (1/12)*m4*(L2^2) + (1/4)*m4*(r3^2);
yy = xx;
zz = (1/2)*m4*(r4^2);
IRAxle_1 = [xx 0 0; 0 yy 0; 0 0 zz];

% derivative of linear(p) and angular(h) momentum of rotating axle in frame
% 2 about rotating axle COM
p_rotating_C3_dot = m4*rOC3_1_ddot;

h_rotating_C3 = IRAxle_1*w1_1;
h_rotating_C3_dot = diff(h_rotating_C3, t) + cross(w1_1, h_rotating_C3);


%%%% Gyroscope Rod stuff
% Inertia Tensor of counterweight about its own COM
xx = (1/2)*m2*(r2^2);
yy = (1/12)*m2*(t2^2) + (1/4)*m2*(r2^2);
zz = yy;
ID2_2 = [xx 0 0; 0 yy 0; 0 0 zz];

% Intertia tensor of just the axle about its own COM
xx = (1/2)*m3*(r3^2);
yy = (1/12)*m3*(L^2) + (1/4)*m3*(r3^2);
zz = yy;
IAxle_2 = [xx 0 0; 0 yy 0; 0 0 zz];

% Inertia tensor of counterweight about the COM of the combined gyroscope axle (C2) in frame {2}
delta_x = (rAD2_2(1) - rAC2_2(1));
delta_y = (rAD2_2(2) - rAC2_2(2));
delta_z = (rAD2_2(3) - rAC2_2(3));
PAT = [delta_y^2 + delta_z^2 -delta_x*delta_y -delta_x*delta_z; -delta_x*delta_y delta_x^2 + delta_z^2 -delta_y*delta_z; -delta_x*delta_z -delta_y*delta_z delta_x^2 + delta_y^2];
ID2_C2_2 = ID2_2 + m2*PAT;

% Inertia tensor of gyro axle about the COM of the combined gyroscope axle (C2) in frame {2}
delta_x = (rACaxle_2(1) - rAC2_2(1));
delta_y = (rACaxle_2(2) - rAC2_2(2));
delta_z = (rACaxle_2(3) - rAC2_2(3));
PAT = [delta_y^2 + delta_z^2 -delta_x*delta_y -delta_x*delta_z; -delta_x*delta_y delta_x^2 + delta_z^2 -delta_y*delta_z; -delta_x*delta_z -delta_y*delta_z delta_x^2 + delta_y^2];
IAxle_C2_2 = IAxle_2 + m3*PAT;

% Inertia tensor of combined counterweight and gyro axle about the COM of combined gyroscope axle (C2) in frame {2}
ICombAxle_C2_2 = IAxle_C2_2 + ID2_C2_2;

% derivative of linear(p) and angular(h) momentum of gyroscope axle in frame
% 2 about gyroscope axle COM
p_Axle_C2_2_dot = (m2 + m3)*rAC2_2_ddot;

h_Axle_C2 = ICombAxle_C2_2*w2_2;
h_Axle_C2_dot = diff(h_Axle_C2, t) + cross(w2_2, h_Axle_C2);


%%%% Disc stuff
% Inertia tensor of disc about its own COM
xx = (1/2)*m1*(r1^2);
yy = (1/12)*m1*(t1^2) + (1/4)*m1*(r1^2);
zz = yy;
ID1_3 = [xx 0 0; 0 yy 0; 0 0 zz];

% derivative of linear and angualr momentum of rotating disc in frame 3
p_D1_3_dot = m1*R32*rAD1_2_ddot;
h_D1_3 = ID1_3*w3_3;
h_D1_3_dot = diff(h_D1_3, t) + cross(w3_3, h_D1_3);

Fg3_0 = [0; 0; m1*9.8];
F3_3 = p_D1_3_dot - R30*Fg3_0;

MG3_3 = h_D1_3_dot;

Fg2_0 = [0; 0; (m2+m3)*9.81];
F2_2 = p_Axle_C2_2_dot - R20*Fg2_0 + R23*F3_3;

MG2_2 = h_Axle_C2_dot + R23*MG3_3 - cross(rAC2_2, F2_2) + cross(rAD1_2, R23*F3_3);

Fg1_0 = [0; 0; m4*9.81];
F1_1 = p_rotating_C3_dot - R10*Fg1_0 + R12*F2_2;

MG1_1 = h_rotating_C3_dot + R12*MG2_2 - cross(rOC3_1, F1_1) + cross(rAO_1, F2_2);


syms F3x F3y F3z tau_3 M3y M3z
[F3x, F3y, F3z] = solve([F3x;F3y;F3z] == F3_3, [F3x, F3y, F3z])

[tau_3, M3y, M3z] = solve([tau_3; M3y; M3z] == MG3_3, [tau_3, M3y, M3z]); 
tau_3

syms F2x F2y F2z tau_2 M2x M2z
[F2x,F2y,F2z]=solve([F2x;F2y;F2z]==F2_2, [F2x,F2y,F2z]); 


[M2x, tau_2, M2z] = solve([M2x; tau_2; M2z] == MG2_2, [M2x, tau_2, M2z]); 
tau_2

syms F1x F1y F1z tau_1 M1x M1y
[F1x,F1y,F1z]=solve([F1x;F1y;F1z]==F1_1, [F1x,F1y,F1z]);

[M1x,M1y,tau_1]=solve([M1x;M1y;tau_1]==MG1_1,[M1x,M1y,tau_1]); 
tau_1

tau1 = 0; 
tau2 = 0; 
tau3 = 0; 
tau_1 = tau1-tau_1 ==0; 
tau_2 = tau2 -tau_2 ==0; 
tau_3 = tau3 - tau_3 == 0; 
eqns= [tau_1, tau_2, tau_3]; 

syms alpha_new beta_new gamma_new alpha_dot alpha_ddot beta_dot beta_ddot gamma_dot gamma_ddot 
old_syms = [diff(alpha, t, 2), diff(alpha, t), alpha, diff(beta,t, 2), diff(beta, t), beta, diff(gamma, t, 2), diff(gamma, t) , gamma] 
new_syms = [alpha_ddot, alpha_dot, alpha_new, beta_ddot, beta_dot, beta_new, gamma_ddot, gamma_dot, gamma_new] 
eqns_new = subs(eqns, old_syms, new_syms)

[A ,b] = equationsToMatrix(eqns_new, [alpha_ddot, beta_ddot, gamma_ddot]) 
X_ddotNew = A\b; 
X_ddotNew = simplify(X_ddotNew)

L_num =  0.48; %- length of gyroscope axle
L2_num = 0.195; %- length of rotating axle
d1_num = 0.128; %- length from A to start of big disc (B)
d2_num = 0.2175;  %- length from A to start of counterweight
t1_num = 0.0221;  %- thickness of disc
t2_num = 0.0313; %- thickness of counterweight
m1_num = 1.747;  %- mass of large disc
m2_num = 0.9; % - mass of counterweight
m3_num = 0.15; % - mass of axle
m4_num = 0.15; %- mass of rotating axle
r1_num = 0.254/2;  %- radius of big disc
r2_num = 0.07/2;% - radius of counterweight
r3_num = 12.1/2/1000;  %- radius of gyroscope axle
r4_num = 12.1/2/1000; % - radius of rotating axle

sym_measurements = [L L2 d1 d2 t1 t2 m1 m2 m3 m4 r1 r2 r3 r4]; 
numerical_measurements =[L_num L2_num d1_num d2_num t1_num t2_num m1_num m2_num m3_num m4_num r1_num r2_num r3_num r4_num]; 
%I should be subbing in for the new equations
% Let's try something that makes a little more sense
% Should I just subinto Xdot? That would make a lot more sense

X_ddotNew = subs(X_ddotNew, sym_measurements, numerical_measurements) 
filename = "X_ddotNew_data.mat"; 
save(filename, "X_ddotNew"); 

tau_1=subs(tau_1, sym_measurements, numerical_measurements);
tau_2=subs(tau_2, sym_measurements, numerical_measurements);
tau_3=subs(tau_3, sym_measurements, numerical_measurements);

[X_dot,S] = odeToVectorField(tau_1, tau_2, tau_3);
S
eom = matlabFunction(X_dot,'Vars',{'t','Y'})

options = odeset('RelTol',1e-7,'AbsTol',1e-7'); % solver options
alpha0=0;
alpha_dot0=-0.8;
beta0=0.122;
beta_dot0=-0.5;
gamma0 = 0; 
gamma_dot0 = 34.27; 
x_init = [alpha0; alpha_dot0; beta0; beta_dot0;  gamma0; gamma_dot0];   % initial condition

tspan=[0 20]; 
EOM_Solved = ode45(eom,tspan,x_init,options);

%% EVAULATE THE SOLUTION
dt = 0.05;                                  % set time step     
t = tspan(1):dt:tspan(2);                   % creat time vector
X = deval(EOM_Solved,t);                           % deval


%% PLOT THE STATES
% Plot for alpha and alpha_dot
subplot(3, 2, 1);
plot(t, X(1,:), 'b', t, X(2,:), 'r');
xlabel('Time');
ylabel('Alpha States');
legend('\alpha', '\dot{\alpha}');
title('Alpha States');

% Plot for gamma and gamma_dot
subplot(3, 2, 3);
plot(t, X(5,:), 'b', t, X(6,:), 'r');
xlabel('Time');
ylabel('Gamma States');
legend('\gamma', '\dot{\gamma}');
title('Gamma States');

% Plot for beta and beta_dot
subplot(3, 2, 5);
plot(t, X(3,:), 'b', t, X(4,:), 'r');
xlabel('Time');
ylabel('Beta States');
legend('\beta', '\dot{\beta}');
title('Beta States');

% Overall plot settings
sgtitle('State Variables over Time');


%% ANIMATION
%% CREATE CYLINDERS
VIDEO = 1; 
%Gyroscope Axel - this is in frame 2
[zg,yg,xg] = cylinder(r3_num,100);
xg = xg*L_num - (L_num-d1_num); 


%Gyroscope CounterWeight - this is in frame 2
[zc, yc, xc] = cylinder(r2_num,100);
xc = xc* t2_num - t2_num - d2_num;

%Rotating Disk - this is in frame 3
[yd, zd, xd] = cylinder(r1_num, 100);
xd = xd*t1_num + d1_num; 

%% SETUP VIDEO IF REQUIRED
if VIDEO
    fps = 1/dt;
    MyVideo = VideoWriter('gyro_animation','MPEG-4');
    MyVideo.FrameRate = fps;
    open(MyVideo);
end

%% CREATE ANIMATION
handle = figure;
hold on % ; grid on
for i = 1:length(t)
    cla 

    alpha = X(1,i)/40;
    beta = X(3, i)/500; 
    gamma = X(5, i); 
   
    % Positive rotation around z axis.
    R01 = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
    R10 = R01';

    % Negative rotation around y axis.
    R12 = [cos(beta) 0 -sin(beta); 0 1 0; sin(beta) 0 cos(beta)];
    R21 = R12.';

    % Positive rotation around x axis.
    R23 = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
    R32 = R23.';

    R20 = R21 * R10;
    R02 = R20.'; 
    R30 = R32 * R20;
    R03 = R30.'; 
   
    %Rotate the GyroScope axel

    [xg, yg, zg] = rotation(xg, yg, zg, R02); 

    %Rotate the Counterweight
    [xc, yc, zc] = rotation(xc, yc, zc, R02); 

    %Rotate the Rotating Disk
    [xd, yd, zd] = rotation(xd, yd, zd, R02); 
    %Surf colours the surface of a 3 dimensional plot
    surf(xg, yg, zg, xg, 'EdgeColor', 'none'); 
    surf(xc, yc, zc, xc, 'EdgeColor', 'none'); 
    surf(xd, yd, zd, xd, 'EdgeColor', 'none');
    
    axis square
    view(3)
    axis(1*[-0.5 0.5 -0.5 0.5 -0.5 0.5])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    if VIDEO
        writeVideo(MyVideo,getframe(handle));
    else
        pause(dt)
    end
end

if VIDEO
close(MyVideo)
end
function [Xf,Yf,Zf]=rotation(Xi,Yi,Zi,R)

I=size(Xi,1);
J=size(Xi,2);

Xf=zeros(I,J);
Yf=zeros(I,J);
Zf=zeros(I,J);

for ii=1:I
    for jj=1:J
        vector=[Xi(ii,jj);Yi(ii,jj);Zi(ii,jj)];
        vector=R*vector;
            Xf(ii,jj)=vector(1);
            Yf(ii,jj)=vector(2);
            Zf(ii,jj)=vector(3);
    end
end

end

%}
%{
function xdot = eom(t, x)
    %{
    syms alpha_new beta_new gamma_new alpha_dot beta_dot gamma_dot 
    sym_param = [alpha_dot alpha_new beta_dot beta_new gamma_dot gamma_new];
    num_param = [x(4) x(1) x(5) x(2) x(6) x(5)];
    %}
    alpha_new = x(1); 
    beta_new = x(5); 
    gamma_new = x(3); 
    alpha_dot = x(2); 
    beta_dot = x(6); 
    gamma_dot = x(4); 
    % Load data and handle potential errors
    try
        alpha_ddot = -(beta_dot*((760788801*gamma_dot*cos(beta_new))/5000000000 - (589361080059*alpha_dot*sin(2*beta_new))/1000000000000))/((665439960159*cos(beta_new)^2)/1000000000000 + 11907/2000000);
        %alpha_ddot = subs(alpha_ddot, sym_param, num_param); 
        alpha_ddot = double(alpha_ddot);
        beta_ddot = (10116922404000*cos(beta_new))/74599273351 - (57031355551*alpha_dot^2*sin(2*beta_new))/149198546702 + (16906417800*alpha_dot*gamma_dot*cos(beta_new))/74599273351;
        %beta_ddot = subs(beta_ddot, sym_param, num_param);
        beta_ddot = double(beta_ddot); 
        gamma_ddot = (beta_dot*cos(beta_new)*((760788801*gamma_dot*sin(beta_new))/5000000000 - (592337830059*alpha_dot)/500000000000 + (513282199959*alpha_dot*cos(beta_new)^2)/1000000000000))/((665439960159*cos(beta_new)^2)/1000000000000 + 11907/2000000); 
        %gamma_ddot = subs(gamma_ddot, sym_param, num_param);
        gamma_ddot = double(gamma_ddot); 
    catch
        error('Error loading or processing data from file.');
    end
    
    % Display loaded data for debugging
    disp('Loaded acceleration data:');
    disp(alpha_ddot);
    disp(beta_ddot);
    disp(gamma_ddot);
    % These are still being recorded as sym expressions
    % To fix, I gotta use subs! And sub the values of alpha_new, 
    % beta_new, gamma_new, and there derivatives for X(1), X(2)
    % X(3) X(4) X(5) X(6) the above I thought I was doing is but not really
    % Assign derivatives to xdot
    xdot = [alpha_dot;
            beta_dot;
            gamma_dot;
            alpha_ddot;
            beta_ddot;
            gamma_ddot];
end
%}