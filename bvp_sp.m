clear;

n20 = 20;

XI(n20); 
APHI(n20);  
G(n20);
H(n20,n20);
UGRAPH(101);

PHI(n20);
PHIX(n20);
PHIXX(n20);
% specify parameters
ALPHA = 1.0;    % u[-1]
BETA = -0.5;   % u[1]
N = 22;         % No. of Chebyshev polynomials
NBASIS = N - 2; % No. of basis funcs phi_j, phi(+-1)=0

% funcs
D0 = @(X) 1;
D1 = @(X) 1;
D2 = @(X) 1;
F  = @(X) X;

% compute the interior collocation points and the forcing vector
for I=1:NBASIS
    X(I) = cos(pi * I / (NBASIS+1)); % Cheb points
    X = XI(I);
    B = ALPHA*(1-X)/2. + BETA*(1+X)/2.;
    BX = (-ALPHA + BETA)/2.;
    G(I) = F(X) - D0(X)*B - D1(X)*BX;  
end

% compute the square matrix
for I=1:NBASIS
    X = XI(I);
    
end




















