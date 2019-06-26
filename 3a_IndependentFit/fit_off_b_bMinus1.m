%This group of scripts is used to fit three components in matrix K
%independently. This script is to fit one off-diagonal component.
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
%In this script, only rho_b_bMinus1, drho_b_bMinus1, off_b_bMinus1 are used.

fitID=22; %fit options

%% fit offdiagonal term
[xData, yData, zData] = prepareSurfaceData(rho_b_bMinus1, drho_b_bMinus1, off_b_bMinus1);

% Set up fittype and options.
if fitID==22
    ft = fittype( 'poly22' );
elseif fitID==21
    ft = fittype( 'poly21' );
end

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft );

% Plot fit with data.
figure
h = plot( fitresult, [xData, yData], zData );

% Label axes
ylabel('$\nabla \rho$','Interpreter','latex','FontSize',20)
xlabel('$\rho$','Interpreter','latex','FontSize',20)
zlabel('$\langle K_{\rho_{b}+\nabla \rho_{b}(x-x_{b})} \gamma_{b}, \gamma_{b-1} \rangle$','Interpreter','latex','FontSize',20);
mkdir('figures');
savefig('figures/off_b_bMinus1.fig')
save(['off_b_bMinus1_fit_',num2str(fitID),'.mat'],'fitresult')

%%  compute analytical value
rhoNew = [4:0.1:10];
drhoNew = [-19:0.1:19];

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
        %Analytically, <K \gamma_b, \gamma_{b-1}> = -m(\rho_b)*Ngamma + 1/2*\nabla m(\rho_b) + O(1/Ngamma)
        Zanalytical(j,i) = -m*Ngamma + 1/2*dmdrho*drhoNew(j);
    end
    clear m 
end
clear tolerance lbound rbound targetrho mid

%% compute fit value
fitGradTerm00 = fitresult.p00;
fitGradTerm10 = fitresult.p10;
fitGradTerm01 = fitresult.p01;
fitGradTerm11 = fitresult.p11;
fitGradTerm20 = fitresult.p20;
if fitID==21
    fitGradTerm02 = 0;
elseif fitID==22
    fitGradTerm02 = fitresult.p02;
end

for i = 1:length(rhoNew)
    for j = 1:length(drhoNew)
        Zfit(j,i)=fitGradTerm00 + fitGradTerm10*rhoNew(i) + fitGradTerm01*drhoNew(j) +...
            fitGradTerm20*rhoNew(i)^2 + fitGradTerm11*rhoNew(i)*drhoNew(j) + fitGradTerm02*drhoNew(j)^2 ;
        error(j,i) = (Zfit(j,i)-Zanalytical(j,i))/Zanalytical(j,i) * 100;
    end
end
clear fitGradTerm00 fitGradTerm10 fitGradTerm01 fitGradTerm11 fitGradTerm20 fitGradTerm02

%% plot error
figure
[X,Y] = meshgrid(rhoNew,drhoNew);
h = surf(X, Y, error);
ylabel('$\nabla \rho$','Interpreter','latex','FontSize',20)
xlabel('$\rho$','Interpreter','latex','FontSize',20)
zlabel('error of $\langle K_{\rho_{b}+\nabla \rho_{b}(x-x_{b})} \gamma_{b-1}, \gamma_{b} \rangle\ [\%]$','Interpreter','latex','FontSize',20);
set(gca,'FontSize',20);
h.EdgeColor = 'none';
box on
clear h 
mkdir('errorPlot')
savefig(['errorPlot/error_off_bMinus1_b_',num2str(fitID),'.fig'])