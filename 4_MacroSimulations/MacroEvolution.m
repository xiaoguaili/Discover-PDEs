% This script is used to compute evolution.

clear all
clc

constraint = 0; % 0: read fit result from independent fit
                % 1: read fit result from fit with conservation constraint
              
boundaryCondition = 'periodic'; % 'periodic': periodic boundary condition
                                 % 'Dirichlet': Dirichlet boundary condition

%% coarse-graining discretization

Ngamma = 40;
dx = 1/Ngamma;
dt = dx^2/1000;                  % time discretization
T = 0.01;                       % total simulation time
timesteps = T/dt;               % total simulation steps
%Domain
xfinal = 1.0;
if strcmp(boundaryCondition,'periodic')
    x = [0:dx:xfinal-dx];
elseif strcmp(boundaryCondition,'Dirichlet')
    x = [0:dx:xfinal];
else
    print('Uncompatible boundary condition!')
end
N = length(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial profile
% frequency = 3;
% alpha = 5;
% beta = 7;
frequency = 2;
alpha = 3;
beta = 7;
density_0 = -alpha*cos(frequency*pi*x/xfinal)+beta; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Video setting
dump = 100; % number of frames
mkdir('movies')
v = VideoWriter(['movies/Movie3'], 'MPEG-4') % video name
v.FrameRate = 5; % frame speed
open(v)
fig = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output snapshots
mkdir('snapshots')
snapshots=[1 3 6 10]; % specify time to output snapshots; real time is snapshots/1000
s = 1;  % number of snapshots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mass matrix
if strcmp(boundaryCondition,'periodic')
    M_0 = gallery('tridiag',N,1/(6),2/(3),1/(6));
    M_0(1,N) = 1/(6);
    M_0(N,1) = 1/(6);
    M = M_0 / Ngamma;
elseif strcmp(boundaryCondition,'Dirichlet')
    % mass matrix partition: first and last entries in rho are
    % fixed, and the other entries in rho are unknown. Therefore, M
    % could be organized and be separated as four sections:
    % [Muu, Muk; Mku, Mkk]
    M_0 = gallery('tridiag',N,1/(6),2/(3),1/(6));
    M = M_0 / Ngamma;
    Muu=M(2:N-1,2:N-1);  
else
    print('Uncompatible boundary condition!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nv = 0; % number of frame
t = 0;

%% import particle data
if constraint == 1
    load(['../3b_FitWithConstraint/fitParameter.mat']);
else
    fitID = 22;    
    load(['../3a_IndependentFit/dia_fit_',num2str(fitID),'.mat']);
    fit00 = fitresult.p00;
    fit10 = fitresult.p10;
    fit01 = fitresult.p01;
    fit11 = fitresult.p11;
    fit20 = fitresult.p20;
    fit02 = fitresult.p02;
    clear fitresult
    
    load(['../3a_IndependentFit/off_b_bPlus1_fit_',num2str(fitID),'.mat']);
    fitPlus00 = fitresult.p00;
    fitPlus10 = fitresult.p10;
    fitPlus01 = fitresult.p01;
    fitPlus11 = fitresult.p11;
    fitPlus20 = fitresult.p20;  
    fitPlus02 = fitresult.p02;
    clear fitresult
    
    load(['../3a_IndependentFit/off_b_bMinus1_fit_',num2str(fitID),'.mat']);
    fitMinus00 = fitresult.p00;
    fitMinus10 = fitresult.p10;
    fitMinus01 = fitresult.p01;
    fitMinus11 = fitresult.p11;
    fitMinus20 = fitresult.p20;
    fitMinus02 = fitresult.p02;
    clear fitresult
end

%% output first snapshot

rho_oldAnal = density_0'; % initialize analytic solution
rho_oldPart = density_0'; % initialize particle-based solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(boundaryCondition,'periodic')
    plot([x xfinal],[rho_oldAnal; rho_oldAnal(1)],'-','Color',[228/255,158/255,37/255],'LineWidth',3,'DisplayName','analytic solution');
    hold on
    plot([x xfinal],[rho_oldPart; rho_oldPart(1)],'--','Color',[16/255,115/255,176/255],'LineWidth',3,'DisplayName','particle-based solution');
    hold off  
elseif strcmp(boundaryCondition,'Dirichlet')
    plot(x,rho_oldAnal,'-','Color',[228/255,158/255,37/255],'LineWidth',3,'DisplayName','analytic solution');
    hold on
    plot(x,rho_oldPart,'--','Color',[16/255,115/255,176/255],'LineWidth',3,'DisplayName','particle-based solution');
    hold off
end
axis([0, xfinal, beta-alpha-1, beta+alpha+1]);
set(gca,'FontSize',20);
legend('Location','south')
box on
title(['time = ',num2str(t)])
xlabel('$x$','Interpreter','latex','FontSize',20)
ylabel('$\rho$','Interpreter','latex','FontSize',20)
set(0,'defaultfigurecolor',[1 1 1])
frame = getframe(fig);      
writeVideo(v, frame);
% output data of first snapshot
save('snapshots/0.mat','x','rho_oldAnal','rho_oldPart')

%% prepare analytic m(rho) at different rho for interpolation

[r, m, dmdr] = prepare_analytic_m(beta-alpha, beta+alpha); % input: lower and upper bound of m

%% compute evolution

while (t <= T)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize K and F
    % particle-based
    KPart = zeros(N,N);
    F_oldPart = zeros(N,1);
    if strcmp(boundaryCondition,'periodic')
        drhoPart = grad(rho_oldPart)/dx; 
    elseif strcmp(boundaryCondition,'Dirichlet')
        drhoPart = gradient(rho_oldPart)/dx;
    else
        print('Uncompatible boundary condition!')
    end
    % analytic
    KAnal = zeros(N,N);
    F_oldAnal = zeros(N,1); 
    if strcmp(boundaryCondition,'periodic')
        drhoAnal = grad(rho_oldAnal)/dx;
    elseif strcmp(boundaryCondition,'Dirichlet')
        drhoAnal = gradient(rho_oldAnal)/dx;
    else
        print('Uncompatible boundary condition!')
    end

    % analytic
    mAnal=interp1(r,m,real(rho_oldAnal),'spline');
    dmAnal = interp1(r,dmdr,real(rho_oldAnal),'spline');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:N
        % compute analytic F = -log(2*m(\rho))
        % particle-based
        F_oldPart(i) = - log(2*interp1(r,m,real(rho_oldPart(i)),'spline')); 
        % analytic
        F_oldAnal(i) = - log(2*interp1(r,m,real(rho_oldAnal(i)),'spline'));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % define indices according to periodic boundary condition
        if i==1
            ip=2;
            im=N;
        elseif i==N
            ip=1;
            im=N-1;
        else
            ip=i+1; % next
            im=i-1; % previous
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define K
        % particle-based
        if constraint == 1
            %Note: K_{b,a} entry: <K r_a, r_b> expanded on a.
            %K_{i,i} entry: <K r_i, r_i> expanded on i (a=i, b=i)
            KPart(i,i) = a(1)*rho_oldPart(i)^2 + a(2)*rho_oldPart(i)*drhoPart(i) +...
                a(3)*drhoPart(i)^2 + a(4)*rho_oldPart(i) + a(5)*drhoPart(i) +a(6);
            %K_{i+1,i} entry: <K r_i, r_{i+1}> expanded on i (a=b-1=i, b=i+1)
            KPart(ip,i) = a(7)*rho_oldPart(i)^2 + a(8)*rho_oldPart(i)*drhoPart(i) +...
                a(9)*drhoPart(i)^2 + a(10)*rho_oldPart(i) + a(11)*drhoPart(i) +a(12);
            %K_{i-1,i} entry: <K r_i, r_{i-1}> expanded on i (a=b+1=i, b=i-1)
            KPart(im,i) = (-a(1)-a(7))*rho_oldPart(i)^2 + (-a(2)-a(8))*rho_oldPart(i)*drhoPart(i) +...
                (-a(3)-a(9))*drhoPart(i)^2 + (-a(4)-a(10))*rho_oldPart(i) +...
                (-a(5)-a(11))*drhoPart(i) + (-a(6)-a(12));
        else
            %K_{i,i} entry: <K r_i, r_i> expanded on i (a=i, b=i)
            KPart(i,i) = (fit00 + fit10*rho_oldPart(i) + fit01*drhoPart(i) +...
                     fit20*rho_oldPart(i)^2 + fit11*rho_oldPart(i)*drhoPart(i) +...
                     fit02*drhoPart(i)^2);
            %K_{i+1,i} entry: <K r_i, r_{i+1}> expanded on i (a=b-1=i, b=i+1)
            KPart(ip,i) = (fitPlus00 + fitPlus10*rho_oldPart(i) + fitPlus01*drhoPart(i) +...
                     fitPlus20*rho_oldPart(i)^2 + fitPlus11*rho_oldPart(i)*drhoPart(i) +...
                     fitPlus02*drhoPart(i)^2);
            %K_{i-1,i} entry: <K r_i, r_{i-1}> expanded on i (a=b+1=i, b=i-1)
            KPart(im,i) = (fitMinus00 + fitMinus10*rho_oldPart(i) + fitMinus01*drhoPart(i) +...
                     fitMinus20*rho_oldPart(i)^2 + fitMinus11*rho_oldPart(i)*drhoPart(i) +...
                     fitMinus02*drhoPart(i)^2);
        end
        %analytic
        %K_{i,i} entry: <K r_i, r_i> expanded on i (a=i, b=i)
        KAnal(i,i) = mAnal(i) * 2 * Ngamma;
        %K_{i+1,i} entry: <K r_i, r_{i+1}> expanded on i (a=b-1=i, b=i+1)
        KAnal(ip,i) = -mAnal(i)*Ngamma - dmAnal(i)*drhoAnal(i)/2;
        %K_{i-1,i} entry: <K r_i, r_{i-1}> expanded on i (a=b+1=i, b=i-1)
        KAnal(im,i) = -mAnal(i)*Ngamma + dmAnal(i)*drhoAnal(i)/2;        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(boundaryCondition,'periodic')
        % particle-based
        rho_newPart=rho_oldPart+M\(dt*KPart*F_oldPart );
        rho_oldPart=rho_newPart;    
        %analytic
        rho_newAnal=rho_oldAnal+M\(dt*KAnal*F_oldAnal );
        rho_oldAnal=rho_newAnal;
    elseif strcmp(boundaryCondition,'Dirichlet')
        % correct nonperiodic stiffness matrix
        KPart(1,N)=0;
        KPart(N,1)=0;
        KAnal(1,N)=0;
        KAnal(N,1)=0;
    
        % stiffness matrix partition: first and last entries in rho are
        % fixed, and the other entries in rho are unknown. Therefore, K
        % could be organized and be separated as four sections:
        % [Kuu, Kuk; Kku, Kkk], and F could be organized and be separated
        % as two sections: [Fu; Fk];
        KuuPart=KPart(2:N-1,2:N-1);
        KukPart=[(KPart(2:N-1,1))'; (KPart(2:N-1,N))'];
        KukPart=KukPart';
        KuuAnal=KAnal(2:N-1,2:N-1);
        KukAnal=[(KAnal(2:N-1,1))'; (KAnal(2:N-1,N))'];
        KukAnal=KukAnal';
    
        Fu_oldPart=F_oldPart(2:N-1);
        Fk_oldPart=[F_oldPart(1);F_oldPart(N)] ;
        Fu_oldAnal=F_oldAnal(2:N-1);
        Fk_oldAnal=[F_oldAnal(1);F_oldAnal(N)] ;
    
        % coarse_graining
        rho_newPart(2:N-1,1)=rho_oldPart(2:N-1,1)+(Muu\(dt*KukPart*Fk_oldPart + dt*KuuPart*Fu_oldPart));
        rho_newPart(1)=density_0(1);
        rho_newPart(N)=density_0(end);
        rho_oldPart=rho_newPart;  
        %analytic
        rho_newAnal(2:N-1,1)=rho_oldAnal(2:N-1,1)+ (Muu\(dt*KukAnal*Fk_oldAnal + dt*KuuAnal*Fu_oldAnal));
        rho_newAnal(1)=density_0(1);
        rho_newAnal(N)=density_0(end);
        rho_oldAnal=rho_newAnal;
    else
        print('Uncompatible boundary condition!')
    end
    
    t=t+dt;
    Nv=Nv+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot figures
    if(rem(Nv,round(timesteps/dump))==0)
        if strcmp(boundaryCondition,'periodic')
            plot([x xfinal],[rho_newAnal; rho_newAnal(1)],'-','Color',[228/255,158/255,37/255],'LineWidth',3,'DisplayName','analytic solution');
            hold on
            plot([x xfinal],[rho_newPart; rho_newPart(1)],'--','Color',[16/255,115/255,176/255],'LineWidth',3,'DisplayName','particle-based solution');
            hold off  
        elseif strcmp(boundaryCondition,'Dirichlet')
            plot(x,rho_newAnal,'-','Color',[228/255,158/255,37/255],'LineWidth',3,'DisplayName','analytic solution');
            hold on
            plot(x,rho_newPart,'--','Color',[16/255,115/255,176/255],'LineWidth',3,'DisplayName','particle-based solution');
            hold off
        end
        axis([0,xfinal,beta-alpha-1,beta+alpha+1]);
        set(gca,'FontSize',20);
        legend('Location','south')
        title(['time = ',num2str(t)])
        xlabel('$x$','Interpreter','latex','FontSize',20)
        ylabel('$\rho$','Interpreter','latex','FontSize',20)
        box on
        set(0,'defaultfigurecolor',[1 1 1])
        frame=getframe(fig);
        writeVideo(v,frame);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output snapshots
    if abs(t-snapshots(s)/1000)<dt
        save(['snapshots/',num2str(snapshots(s)),'.mat'],'x','rho_newAnal','rho_newPart')
        if s<length(snapshots)
            s=s+1;
        end
    end
end
close(v);

