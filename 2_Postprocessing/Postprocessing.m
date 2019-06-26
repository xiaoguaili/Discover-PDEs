%This matlab code reads in the data from zero range process, and compute
%each component in K. The output of this script is saved in file 'Summary.m'

clear all
clc

feval('Summary'); % read in Summary file, including all information about particle simulation
dx=1/Nbin;

diagonal=zeros(Ngamma,1);
offdiagonal=zeros(Ngamma,1);
rho=zeros(Ngamma,1);

for j=2:2  % 1 is at boundary.
    rho_gamma_initial=[];
    rho_gamma_final=[];
    for req=1:R1
        req
        run(['Data_',num2str(req-1)]);
        rho_gamma_initial=[rho_gamma_initial;rho_gamma_initial_final(1:R2,2*(j-1)+1)];
        rho_gamma_final=[rho_gamma_final;rho_gamma_initial_final(1:R2,2*(j-1)+2)];
    end    
    rho_g_ave_initial=mean(rho_gamma_initial);
    rho_g_ave_final=mean(rho_gamma_final);
    %Define basis function
    gamma=max(0,1.0-Ngamma*abs(x-x_basis(j)));  
    gamma_int=sum(gamma)*dx;
    rho(j)=0.5*(rho_g_ave_initial + rho_g_ave_final)/gamma_int;
    %Compute Y_gamma(t_0) and  Y_gamma(t_0+h)
    Y_gamma_initial1=sqrt(Nbin)*(rho_gamma_initial-rho_g_ave_initial);
    Y_gamma_final1=sqrt(Nbin)*(rho_gamma_final-rho_g_ave_final);
    %Compute diagonal component
    diagonal(j)=mean((Y_gamma_final1-Y_gamma_initial1).^2)/h/2;
    diagonal_N=(Y_gamma_final1-Y_gamma_initial1).^2/h/2;   
    diagonal_SD(j)=std(diagonal_N); %Compute standard deviation
    clear diagonal_N rho_gamma_initial rho_gamma_final rho_g_ave_initial rho_g_ave_final
end

for j=3:Ngamma 
    rho_gamma_initial=[];
    rho_gamma_final=[];
    for req=1:R1
        run(['Data_',num2str(req-1)]);
        rho_gamma_initial=[rho_gamma_initial;rho_gamma_initial_final(1:R2,2*(j-1)+1)];
        rho_gamma_final=[rho_gamma_final;rho_gamma_initial_final(1:R2,2*(j-1)+2)];
    end    
    rho_g_ave_initial=mean(rho_gamma_initial);
    rho_g_ave_final=mean(rho_gamma_final);
    %Define basis function
    gamma=max(0,1.0-Ngamma*abs(x-x_basis(j)));  
    gamma_int=sum(gamma)*dx;
    rho(j)=0.5*(rho_g_ave_initial + rho_g_ave_final)/gamma_int; 
    %Compute Y_gamma(t_0) and  Y_gamma(t_0+h)
    Y_gamma_initial2=sqrt(Nbin)*(rho_gamma_initial-rho_g_ave_initial);
    Y_gamma_final2=sqrt(Nbin)*(rho_gamma_final-rho_g_ave_final);
    %Compute diagonal component and off-diagonal component
    diagonal(j)=mean((Y_gamma_final2-Y_gamma_initial2).^2)/h/2;
    diagonal_N=(Y_gamma_final2-Y_gamma_initial2).^2/h/2;   
    diagonal_SD(j)=std(diagonal_N);
    
    offdiagonal(j)=mean((Y_gamma_final2-Y_gamma_initial2).*(Y_gamma_final1-Y_gamma_initial1))/h/2;
    offdiagonal_N=(Y_gamma_final2-Y_gamma_initial2).*(Y_gamma_final1-Y_gamma_initial1)/h/2;   
    offdiagonal_SD(j)=std(offdiagonal_N);
   
    Y_gamma_final1 = Y_gamma_final2;
    Y_gamma_initial1 = Y_gamma_initial2;
    clear diagonal_N offdiagonal_N rho_gamma_initial rho_gamma_final rho_g_ave_initial rho_g_ave_final
end
%% output data
myfile = fopen('Summary.m','a');

% diagonal component
fprintf(myfile,'%% Expected value of \\rho, \\nabla rho, <K \\gamma_b,\\gamma_b> associated to each gamma\n');
%%%%%%%%%%%%%%
fprintf(myfile,'rho_b_b = [%6.6f',rho(2));
for j=3:Ngamma
	fprintf(myfile,', %6.6f',rho(j));
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'drho_b_b = [%6.6f',slope);
for j=3:Ngamma
	fprintf(myfile,', %6.6f',slope);
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'diagonal_b_b = [%6.6f',diagonal(2));
for j=3:Ngamma
	fprintf(myfile,', %6.6f',diagonal(j));
end
fprintf(myfile,'];\n');


% off-diagonal component: <K \gamma_b, \gamma_{b+1}> expanded on \rho_b and \nabla \rho_b
fprintf(myfile,'%% Expected value of \\rho, \\nabla rho, <K \\gamma_b, \\gamma_{b+1}> expanded on \\rho_b and \\nabla \\rho_b\n');
%%%%%%%%%%%%%%
fprintf(myfile,'rho_b_bPlus1 = [%6.6f',rho(2));
for j=3:Ngamma-1
	fprintf(myfile,', %6.6f',rho(j));
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'drho_b_bPlus1 = [%6.6f',slope);
for j=3:Ngamma-1
	fprintf(myfile,', %6.6f',slope);
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'off_b_bPlus1 = [%6.6f',offdiagonal(3));
for j=4:Ngamma
	fprintf(myfile,', %6.6f',offdiagonal(j));
end
fprintf(myfile,'];\n');

% off-diagonal component: <K \gamma_{b}, \gamma_{b-1}> expanded on \rho_b and \nabla \rho_b
fprintf(myfile,'%% Expected value of \\rho, \\nabla rho, <K \\gamma_{b}, \\gamma_{b-1}> expanded on \\rho_b and \\nabla \\rho_b\n');
%%%%%%%%%%%%%%
fprintf(myfile,'rho_b_bMinus1 = [%6.6f',rho(3));
for j=4:Ngamma
	fprintf(myfile,', %6.6f',rho(j));
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'drho_b_bMinus1 = [%6.6f',slope);
for j=4:Ngamma
	fprintf(myfile,', %6.6f',slope);
end
fprintf(myfile,'];\n');
%%%%%%%%%%%%%%
fprintf(myfile,'off_b_bMinus1 = [%6.6f',offdiagonal(3));
for j=4:Ngamma
	fprintf(myfile,', %6.6f',offdiagonal(j));
end
fprintf(myfile,'];\n');


fclose(myfile);

