function [Er,Ar,Br,Cr, Da] = production_engineering_2012()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      LOAD DATA           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 100;  % No. of differential variables
l = 10;   % No. of algebraic variables
nin = 1;  % Number of inputs
nout = 1; % Number of outputs
den = 0.01;
% Define matrices
I = speye(n);
M1 = .5 * I + spdiags(-0.2 * ones(n, 1), 2, n, n) + spdiags(-0.2 * ones(n, 1), -2, n, n) + spdiags(0.2 * ones(n, 1), 4, n, n) + spdiags(0.2 * ones(n, 1), -4, n, n);
K11 = spdiags(5 * ones(n, 1), 0, n, n) + spdiags(-1 * ones(n, 1), 2, n, n) + spdiags(-1 * ones(n, 1), -2, n, n) + spdiags(2 * ones(n, 1), 4, n, n) + spdiags(2 * ones(n, 1), -4, n, n);
mu = 0.005;
nu = .1;
D1 = mu * M1 + nu * K11;
K22 = spdiags(5 * ones(l, 1), 0, l, l) + spdiags(-1 * ones(l, 1), 2, l, l) + spdiags(-1 * ones(l, 1), -2, l, l) + spdiags(2 * ones(n, 1), 4, l, l) + spdiags(2 * ones(n, 1), -4, l, l);
K12 = sprand(n, l, den);
K21 = K12';
% Matrix assembly
E1 = [speye(n) sparse(n, n); sparse(n, n) M1];
J1 = [sparse(n, n) speye(n); -K11 -D1];
J2 = [sparse(n, l); -K12];
J3 = [-K21 sparse(l, n)];
J4 = -K22;
B11 = spdiags(ones(n, 1), 0, n, nin);
B1 = [sparse(n, nin); B11];
B2 = spdiags(ones(l, 1), 0, l, nin);
C1 = [B11' sparse(nout, n)];
C2 = B2';
% Form the equation system
E = [E1 spalloc(size(J2, 1), size(J2, 2), 0); spalloc(size(J3, 1), size(J3, 2), 0) spalloc(size(J4, 1), size(J4, 2), 0)];
A = [J1 J2; J3 J4];
B = [B1; B2];
C = B';
Da = sparse(size(B2', 1), size(B2, 2));
%% Shift parameter calculation
l = 40;
n1 = size(A, 1);
B_rand = sprand(n1, 100, .8);
[Uf, S, ~] = svd(full(B_rand), 0);
s_vals = diag(S);
tol = max(size(A)) * eps(max(s_vals));
r = sum(s_vals > tol);
Q = Uf(:, 1:r) * S(1:r, 1:r);
An = Q' * A * Q;
En = Q' * E * Q;
rw0 = eig(full(An), full(En));
rw = rw0(real(rw0) < 0);
p = lp_mnmx(rw, l);
%% Solve ADI
maxiter1 = 100;
restol = 10^(-15);
Z = [];
res1 = zeros(1, maxiter1 + 2);
m1 = size(E1, 1);
W = B1 - J2 * (J4 \ B2);
bnorm = norm(W' * W, 'fro');
i = 1;
while i <= maxiter1
    pc = p(mod(i + l - 1, l) + 1);
    if i == 1
        Xc = [J1 + pc * E1 J2; J3 J4] \ [B1; B2];
        Vc = Xc(1:m1, :);
    else
        Xc = [J1 + pc * E1 J2; J3 J4] \ [W; sparse(size(J4, 1), size(B2, 2))];
        Vc = Xc(1:m1, :);
    end
   
    if isreal(pc)
        Z = [Z sqrt(-2 * real(pc)) * real(Vc)];
        W = W - 2 * real(pc) * E1 * Vc;
    else
        beta = real(pc) / imag(pc);
        gam = 2 * sqrt(-real(pc));
        Z = [Z, gam * (real(Vc) + beta * imag(Vc)), gam * sqrt(beta^2 + 1) * imag(Vc)];
        pc = conj(pc);
        Vc = real(Vc) + beta * imag(Vc);
        W = W - 4 * real(pc) * E1 * Vc;
        i = i + 1;
    end
   
    res1(i) = norm(W' * W, 'fro') / bnorm;
    fprintf(1, 'step: %4d  normalized residual: %d\n', i, res1(i));
   
    if res1(i) < restol
        break;
    end
    i = i + 1;
end
%% Balanced Truncation
Zon=E1*Z;
[Uc,S,~] = svd(full(Zon'*Z),0);
s0=diag(S);
ks=length(s0);
K0=ks;
while (sum(s0(K0-1:ks))<tol/2)&&(K0>2)
  K0=K0-1;
end
r=K0;
  sigma_r=diag(S(1:r,1:r));
Vc=Z*Uc(:,1:r);
T=Vc*(diag(ones(r,1)./sqrt(sigma_r)));
Er=T'*E1*T;J1_til=T'*J1*T;
J2_til=T'*J2;
J3_til=J3*T;
Ar=(J1_til-J2_til*(J4\J3_til));
Br=(T'*B1)-(J2_til*(J4\B2));
Cr=(B1'*T)-(B2'*(J4\J3_til));
Dar=Da-B2'*(J4\B2);
%% Sigma plot (optional)
low = -1; % adjust according to needs
up = 1;
points = 100;
s = logspace(low, up, points);
Ho = zeros(1, points);
Hr = zeros(1, points);
abserr = zeros(1, points);
for k = 1:points
    G1 = C * ((1j * s(k) * E - A) \ B) + Da;
    G2 = Cr * ((1j * s(k) * Er - Ar) \ Br) + Dar;
    Ho(k) = max(svds(G1));
    Hr(k) = max(svds(G2));
    abserr(k) = max(svds(G1 - G2));
end
relaterr = abserr ./ Ho;
% Plot results
figure(1);
loglog(s, Ho, 'r', s, Hr, 'b-.');
xlabel('\omega');
ylabel('\sigma_{max}(G(j\omega))');
title('Transfer function of original and reduced order system');
legend('Full model', 'ROM');
figure(2);
loglog(s, abserr, 'b');
xlabel('\omega');
ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega))');
title('Absolute model reduction error');
figure(3);
loglog(s, relaterr, 'b');
xlabel('\omega');
ylabel('\sigma_{max}(G(j\omega) - G_r(j\omega)) / \sigma_{max}(G(j\omega))');
title('Relative model reduction error');

end

function p = lp_mnmx(rw,l0)
 
% if length(rw)<l0
%   error('length(rw) must be at least l0.');
% end
%
max_rr = +Inf;                       % Choose initial parameter (pair)
for i = 1:length(rw)
  max_r = lp_s(rw(i),rw);
  if max_r < max_rr
    p0 = rw(i);
    max_rr = max_r;
  end
end  
if imag(p0)
  p = [ p0; conj(p0) ];
else
  p = p0;                            
end
[max_r,i] = lp_s(p,rw);         % Choose further parameters.
while size(p,1) < l0
   
  p0 = rw(i);
  if imag(p0)
    p = [ p; p0; conj(p0) ];
  else
    p = [ p; p0];
  end
 
  [max_r,i] = lp_s(p,rw);
   
end
end

function [max_r,ind] = lp_s(p,set)
%
% Computation of the maximal magnitude of the rational ADI function over
% a discrete subset of the left complex half plane.
%
%   Calling sequence:
%
%     [max_r,ind] = lp_s(p,set)
%
%   Input:
%
%     p        vector of ADI parameters;
%     set      vector representing the discrete set.
%
%   Output:
%
%     max_r    maximal magnitude of the rational ADI function over set;
%     ind      index - maximum is attained for set(ind).
%
%  
%   LYAPACK 1.0 (Thilo Penzl, Jan 1999)
%
max_r = -1;
ind = 0;
 
for i = 1:length(set)
 
  x = set(i);
 
  rr = 1;
  for j = 1:length(p)
    rr = rr*abs(p(j)-x)/abs(p(j)+x);
   
  end  
   
  if rr > max_r
   
    max_r = rr;
    ind = i;
   
  end
 
end
end
