clc; clear;

% discretize the spatial part
N = 100; dx = 1/(N+1);
x = (1:N)'*dx;
epsilon = 0.01;

% Laplacian matrix
e = ones(N,1);
L = epsilon*(diag(-2*e)+diag(e(1:end-1),1)+diag(e(1:end-1),-1))/dx^2 + eye(N);

% initial value
u0 = 0.1*sin(pi*x);
h = 0.01;
tfinal = 1;
Nt = round(tfinal/h);

% nonlinearly function
Nfun = @(u) -u.^3;

% ETDRK4
h_ref = 0.001; Nt_ref = round(tfinal/h_ref);
u_ref = u0;
E = expm(h_ref*L);
for n = 1:Nt_ref
    k1 = Nfun(u_ref);
    k2 = Nfun(u_ref + 0.5*h_ref*k1);
    k3 = Nfun(u_ref + 0.5*h_ref*k2);
    k4 = Nfun(u_ref + h_ref*k3);
    u_ref = E*u_ref + h_ref/6*(k1+2*k2+2*k3+k4);
end

% compare
methods = {'ETD1','ETD2','ETDRK4'};
errors = zeros(size(methods));

for m = 1:length(methods)
    u = u0;
    for n = 1:Nt
        switch methods{m}
            case 'ETD1'
                u = expm(h*L)*u + h*expm(h*L)*Nfun(u);
            case 'ETD2'
                % Second Taylor
                u = expm(h*L)*u + h*expm(h*L)*(Nfun(u) + 0.5*h* (Nfun(u)-Nfun(u)));
            case 'ETDRK4'
                % 四階 RK
                E = expm(h*L);
                k1 = Nfun(u);
                k2 = Nfun(u + 0.5*h*k1);
                k3 = Nfun(u + 0.5*h*k2);
                k4 = Nfun(u + h*k3);
                u = E*u + h/6*(k1+2*k2+2*k3+k4);
        end
    end
    errors(m) = norm(u-u_ref,inf);
end

% plot
figure; bar(errors);
set(gca,'XTickLabel',methods);
ylabel('Max norm error');
title('Error vs method');
grid on;
