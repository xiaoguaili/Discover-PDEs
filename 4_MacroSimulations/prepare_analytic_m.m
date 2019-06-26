function [r, m, dmdr] = prepare_analytic_m(ml, mu)

r=ml:0.001:mu;

for i = 1:length(r)
    targetrho = r(i);
    error = 10;
    tolerance = 1e-14;
    lbound = 0;
    rbound = 150;
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
    m(i) = mid;
    dmdr(i)=(2-2*besseli(1,2*sqrt(2*m(i)))^2.*besseli(0,2*sqrt(2*m(i)))^(-2))^(-1);

end
clear error i lbound mid rbound targetrho tolerance rho

% figure
% plot(r,m,'LineWidth',3)
% ylabel('$m$','Interpreter','latex','FontSize',20)
% xlabel('$\rho$','Interpreter','latex','FontSize',20)
% set(gca,'FontSize',20); 
% save('analytic_m_dmdrho.mat','r','m','dmdr')
end