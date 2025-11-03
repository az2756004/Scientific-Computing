clear
clc

% %barycentric
% m = 10;
% z_j = randn(m,1); % set of real or complex distinct support points.
% f_j = zeros(m,1); % the set of real or complex data values.
% w_j = randn(m,1); % the set of real or complex weights.
% n_j = zeros(m,1);
% d_j = zeros(m,1);
% 
% % Consider function F：
% 
% for j =1:m
%     f(j) = F(z_j(j));
% end
% 
% % The AAA algorithm

Z = exp(linspace(-.5,.5+.15i*pi,1000));
F = @(z) tan(pi*z/2);

[r,pol,res,zer,z,f,w,err] = aaa(F,Z);


