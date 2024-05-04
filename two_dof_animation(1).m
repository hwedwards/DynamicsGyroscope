function two_dof_animation()
    clear all
    close all
    clc
    %% DO I WANT TO RECORD THE VIDEO
    VIDEO = 1;
    
    %% SETUP THE PROBLEM
    X_init = [0;pi/3;10;0];                         % initial conditions
    tspan = [0 8];                                 % start and finish times
    options = odeset('RelTol',1e-7,'AbsTol',1e-7'); % solver options
    sol = ode45(@eom,tspan,X_init,options);         % SOLVE the eoms

    %% EVAULATE THE SOLUTION
    dt = 0.03;                                  % set time step                        
    t = tspan(1):dt:tspan(2);                   % creat time vector
    X = deval(sol,t);                           % deval

    %% PLOT THE STATES
    plot(t,X)
    xlabel('time')
    ylabel('states')
    h = legend('$\theta$','$\alpha$','$\dot{\theta}$','$\dot{\alpha}$');
    set(h,'Interpreter','latex')

    %% CREATE CYLINDERS
    l1 = 0.3;
    l2 = 0.9;
    [x1,y1,z1] = cylinder(0.1,10);
    z1 = z1*l1;
    [x2,y2,z2] = cylinder(0.1,100);
    z2 = z2*l2;
    z2 = z2-l2;

    %% SETUP VIDEO IF REQUIRED
    if VIDEO
        fps = 1/dt;
        MyVideo = VideoWriter('two_dof_animation','MPEG-4');
        MyVideo.FrameRate = fps;
        open(MyVideo);
    end

    %% CREATE ANIMATION
    handle = figure;
    hold on % ; grid on
    for i = 1:length(t)
        cla

        theta = X(1,i);
        alpha = X(2,i);

        R1to0 = [cos(theta),-sin(theta),0;
                 sin(theta),cos(theta),0;
                 0,0,1];
        R2to1 = [cos(alpha),0,sin(alpha);
                 0 1 0;
                 -sin(alpha),0,cos(alpha)];

        [x1_rotated,y1_rotated,z1_rotated] = rotation(x1,y1,z1,R1to0);
        [x2_rotated,y2_rotated,z2_rotated] = rotation(x2,y2,z2,R1to0*R2to1);
        surf(x1_rotated,y1_rotated,z1_rotated,x1)%'EdgeColor','none'
        surf(x2_rotated,y2_rotated,z2_rotated,x2,'EdgeColor','none')
        axis square
        view(3)
        axis(1*[-1 1 -1 1 -1 1])
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

end

%%
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

%%
function xdot = eom(t,x)

    tau1 = 10*sin(t);
    tau2 = 5*cos(t);

    theta = x(1);
    alpha = x(2);
    theta_d = x(3);
    alpha_d = x(4);
    
    g = 9.8;
       
    l2 = 0.9;
    m2 = 1;
    
    Ix2 = 0.5;
    Iy2 = Ix2;
    Iz2 = 1;
    
    Iz1 = 1;
        
    f1_num = tau1 + alpha_d*theta_d*sin(alpha)*cos(alpha)*(2*(Iz2-Ix2)-m2*l2^2/4);
    f1_den = Iz1 + Iz2*cos(alpha)^2+(Ix2+m2*l2^2/4)*sin(alpha)^2;
    
    f1 = f1_num/f1_den;
    
    f2_num = tau2+theta_d^2*sin(alpha)*cos(alpha)*(Ix2-Iz2+m2*l2^2/4)-m2*l2*g*sin(alpha)/2;
    f2_den = Iy2 + m2*l2^2/4;
    
    f2 = f2_num/f2_den;
    
    xdot = [theta_d;
            alpha_d;
            f1
            f2]; 
end

