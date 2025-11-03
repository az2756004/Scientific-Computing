clear
clc
%% AAA rational approximation (flattened version)
% Make sure F and Z are defined before running this script.

% Example setup (you can modify or comment this out)
% F = @(z) tan(pi*z/2);
% Z = exp(linspace(-.5, .5 + .15i*pi, 1000));

tol  = 1e-13;   % tolerance
mmax = 100;     % maximum iterations

% --- initialization ---
Z = exp(linspace(-.5,.5+.15i*pi,1000));
M  = length(Z);
F = @(x) tan(pi*x/2);
if isa(F, 'function_handle')
    F = F(Z); % evaluate if function handle
end
Z  = Z(:); 
F  = F(:);

SF = spdiags(F, 0, M, M); % left scaling matrix
J  = 1:M; 
z  = []; 
f  = []; 
C  = [];
errvec = []; 
R  = mean(F);

% --- main AAA iteration loop ---
for m = 1:mmax
    [~, j] = max(abs(F - R));   % select next support point
    z = [z; Z(j)];              % add new support point
    f = [f; F(j)];              % corresponding function value
    J(J == j) = [];             % remove from remaining set

    C = [C, 1 ./ (Z - Z(j))];   % add next Cauchy column
    Sf = diag(f);               % right scaling matrix
    A = SF * C - C * Sf;        % Loewner matrix

    [~, ~, V] = svd(A(J,:), 0); % SVD on reduced matrix
    w = V(:, m);                % weight vector = min sing. vector

    N = C * (w .* f);           % numerator
    D = C * w;                  % denominator

    R = F;
    R(J) = N(J) ./ D(J);        % rational approximation at sample pts

    err = norm(F - R, inf);
    errvec = [errvec; err];

    % --- stopping condition ---
    if err <= tol * norm(F, inf)
        break;
    end
end

%% --- Compute poles, residues, zeros (flattened PRZ section) ---
m = length(w);
B = eye(m + 1); 
B(1,1) = 0;
E = [0, w.'; ones(m,1), diag(z)];

pol = eig(E, B);
pol = pol(~isinf(pol));     % poles

% --- define evaluation function for AAA approximant ---
r_eval = @(zz) ( (1./(zz(:) - z.')) * (w.*f) ) ./ ( (1./(zz(:) - z.')) * w );

% --- compute residues ---
dz = 1e-5 * exp(2i * pi * (1:4) / 4);
res = zeros(size(pol));
for k = 1:length(pol)
    res(k) = r_eval(pol(k) + dz(k)) .* (dz(k).') / 4;
end



% --- compute zeros ---
E = [0, (w .* f).'; ones(m,1), diag(z)];
zer = eig(E, B);
zer = zer(~isinf(zer));

%% --- Evaluate rational approximant on arbitrary zz (flattened rhandle) ---
r = @(zz) ((1./(zz(:) - z.')) * (w .* f)) ./ ((1./(zz(:) - z.')) * w);

% --- Fix NaN values (Inf/Inf case) ---
zz_test = Z; % example domain for fixing
r_vals = r(zz_test);
ii = find(isnan(r_vals));
for j = 1:length(ii)
    idx = find(zz_test(ii(j)) == z);
    if ~isempty(idx)
        r_vals(ii(j)) = f(idx);
    end
end
r_vals = reshape(r_vals, size(zz_test));

%% --- Done ---
fprintf('AAA completed: type (%d,%d), final error = %.3e\n', ...
    length(z)-1, length(z)-1, err);

