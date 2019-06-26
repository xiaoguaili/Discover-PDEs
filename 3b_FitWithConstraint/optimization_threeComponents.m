%This script is used to fit three components in matrix K
%with conservation constraint. Practically, sum of the three components is
%zero.
clear all
clc

load('Data.mat')
%Load data of three components from different profiles
%There are 9 variables in Data.mat, which come from Summary.m.
%They are rho_b_b, drho_b_b, diagonal_b_b for diagonal component;
%rho_b_bPlus1, drho_b_bPlus1, off_b_bPlus1 for <K \gamma_b, \gamma_{b+1}>
%expanded on \rho_b and \nabla \rho_b;
%rho_b_bMinus1, drho_b_bMinus1, off_b_bMinus1 for <K \gamma_b, \gamma_{b-1}>
%expanded on \rho_b and \nabla \rho_b.

%% Define fit model
%Least square error
lsq = @(v) sum((v(1)*rho_b_b.^2 + v(2)*rho_b_b.*drho_b_b + v(3)*drho_b_b.^2 + v(4)*rho_b_b + v(5)*drho_b_b + v(6) - diagonal_b_b).^2)+...
   sum((v(7)*rho_b_bPlus1.^2 + v(8)*rho_b_bPlus1.*drho_b_bPlus1 + v(9)*drho_b_bPlus1.^2 + v(10)*rho_b_bPlus1 + v(11)*drho_b_bPlus1 + v(12) - off_b_bPlus1).^2)+...
   sum(((-v(1)-v(7))*rho_b_bMinus1.^2 + (-v(2)-v(8))*rho_b_bMinus1.*drho_b_bMinus1 + (-v(3)-v(9))*drho_b_bMinus1.^2 + (-v(4)-v(10))*rho_b_bMinus1 + (-v(5)-v(11))*drho_b_bMinus1 + (-v(6)-v(12)) - off_b_bMinus1).^2);

a0 = ones(1,12); %initial values
options.TolX = 1e-22; %tolerance
options.TolFun = 1e-20; %tolerance
options.MaxFunEvals = 10000; %maximum iteration
[a, residue] = fminunc(lsq,a0,options)
clear a0 lsq options
save(['fitParameter.mat'],'a')

%% compute analytic m
rhoNew = [4:0.1:10];
drhoNew = [-19:0.1:19];
Ngamma=40;
% find analytic m for each rhoNew and drhoNew
for i = 1:length(rhoNew)
    targetrho = rhoNew(i);
    error = 10;
    tolerance = 2e-15;
    lbound = 0;
    rbound = 100;
    while(error > tolerance && lbound < rbound)
        mid = (lbound + rbound) / 2;
        rho=sqrt(2*mid)*besseli(1,2*sqrt(2*mid))*besseli(0,2*sqrt(2*mid))^(-1);
        error = abs(rho - targetrho);
        if(rho > targetrho && lbound < rbound)
            rbound = mid;
        elseif(rho < targetrho && lbound < rbound) 
            lbound = mid;
        end
    end
    m = mid;
    dmdrho=(2-2*besseli(1,2*sqrt(2*m))^2.*besseli(0,2*sqrt(2*m))^(-2))^(-1);
    for j = 1:length(drhoNew)
        %Analytically, <K \gamma_b, \gamma_b> = 2*m(\rho_b)*Ngamma + O(1/Ngamma)
        anal_dia(j,i) = 2*m*Ngamma;
        %Analytically, <K \gamma_b, \gamma_{b+1}> = -m(\rho_b)*Ngamma - 1/2*\nabla m(\rho_b) + O(1/Ngamma)
        anal_off_b_bPlus1(j,i) = -m*Ngamma - 1/2*dmdrho*drhoNew(j);
        %Analytically, <K \gamma_b, \gamma_{b-1}> = -m(\rho_b)*Ngamma + 1/2*\nabla m(\rho_b) + O(1/Ngamma)
        anal_off_b_bMnus1(j,i) = -m*Ngamma + 1/2*dmdrho*drhoNew(j);
    end
    clear m 
end
clear error 

%% compute fit values
for i = 1:length(rhoNew)
    for j = 1:length(drhoNew)
        fit_dia(j,i)=a(1)*rhoNew(i)^2 + a(2)*rhoNew(i)*drhoNew(j) +...
            a(3)*drhoNew(j)^2 + a(4)*rhoNew(i) + a(5)*drhoNew(j) +a(6);
        fit_off_b_bPlus1(j,i)=a(7)*rhoNew(i)^2 + a(8)*rhoNew(i)*drhoNew(j) +...
            a(9)*drhoNew(j)^2 + a(10)*rhoNew(i) + a(11)*drhoNew(j) +a(12);
        fit_off_b_bMnus1(j,i)=(-a(1)-a(7))*rhoNew(i)^2 + (-a(2)-a(8))*rhoNew(i)*drhoNew(j) +...
            (-a(3)-a(9))*drhoNew(j)^2 + (-a(4)-a(10))*rhoNew(i) +...
            (-a(5)-a(11))*drhoNew(j) + (-a(6)-a(12));
        error_dia(j,i) = (fit_dia(j,i)-anal_dia(j,i))/anal_dia(j,i) * 100;
        error_off_b_bPlus1(j,i) = (fit_off_b_bPlus1(j,i)-anal_off_b_bPlus1(j,i))/anal_off_b_bPlus1(j,i) * 100;
        error_off_b_bMnus1(j,i) = (fit_off_b_bMnus1(j,i)-anal_off_b_bMnus1(j,i))/anal_off_b_bMnus1(j,i) * 100;
    end
end
clear tolerance lbound rbound targetrho mid i j

%% plot errors

figure
[X,Y] = meshgrid(rhoNew,drhoNew);
h = surf(X, Y, error_dia);
ylabel('$\nabla \rho$','Interpreter','latex','FontSize',20)
xlabel('$\rho$','Interpreter','latex','FontSize',20)
zlabel('error of $\langle K_{\rho_{b}+\nabla \rho_{b}} \gamma_b, \gamma_{b} \rangle\ [\%]$','Interpreter','latex','FontSize',20);
set(gca,'FontSize',20);
h.EdgeColor = 'none';
box on
clear h X Y
savefig(['error_dia.fig'])

figure
[X,Y] = meshgrid(rhoNew,drhoNew);
h = surf(X, Y, error_off_b_bPlus1);
ylabel('$\nabla \rho$','Interpreter','latex','FontSize',20)
xlabel('$\rho$','Interpreter','latex','FontSize',20)
zlabel('error of $\langle K_{\rho_{b}+\nabla \rho_{b}} \gamma_b, \gamma_{b+1} \rangle\ [\%]$','Interpreter','latex','FontSize',20);
set(gca,'FontSize',20);
h.EdgeColor = 'none';
box on
clear h X Y
savefig(['error_off_b_bPlus1.fig'])
% 
figure
[X,Y] = meshgrid(rhoNew,drhoNew);
h = surf(X, Y, error_off_b_bMnus1);
ylabel('$\nabla \rho$','Interpreter','latex','FontSize',20)
xlabel('$\rho$','Interpreter','latex','FontSize',20)
zlabel('error of $\langle K_{\rho_{b}+\nabla \rho_{b}} \gamma_b, \gamma_{b-1} \rangle\ [\%]$','Interpreter','latex','FontSize',20);
set(gca,'FontSize',20);
h.EdgeColor = 'none';
box on
clear h X Y
savefig(['error_off_b_bMnus1.fig'])


