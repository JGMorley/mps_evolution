%% Find ground state of g0 model, then quench and compare with exact solution

D = 2; % increase this to get a better approximation. Max val about 8
Vmethod = 'nullspace';

g0 = 1.5;
% Run with imag time to find ground state
isingQuench_imagTime

%output = {A};
disp('GS run evolved ok')

%% Assign final state as new initial state

% Clear variables
clear A0 T expSigmaZ rankrl An V_L g rateFnData d h rho En0 dimNullSpaceL k             
clear times J dt minDim y Ln entropy n yData Rn evalues number_samples S g0
A0 = output{end,1};
clear expSigmaX output      

%% Print initial single-site density matrix
[A_init, R_init, L_init] = normalizeMPS(A0);
rho_init = ncon({A_init, R_init, conj(A_init), L_init},...
                {[1 2 -1], [2 4], [3 4 -2], [1 3]});

%% Hamiltonian tensor

J = 1.;
g0 = 1.5;
g = 0.2;
%>> h_{i,i+1} = -J*Z.Z -(Jg/2)(X.I + I.X)]  (v.1) (X = [0 1;1 0] etc)
h(:,1,:,1) = [-J     -J*g/2;-J*g/2     +J];
h(:,1,:,2) = [-J*g/2      0;     0 -J*g/2];
h(:,2,:,1) = [-J*g/2      0;     0 -J*g/2];
h(:,2,:,2) = [+J     -J*g/2;-J*g/2     -J];


%% Run tdvp
T = 3; 
dt = T/600; % Total time, time step
output = tdvpIntegrator(A0,h,T,dt,'INVERSE_FREE',true,'SCHMIDT_TH',1e-8,'BF_TOL',2);
disp('MPS evolved ok')

%% Initialize arrays to be plotted
number_samples = size(output,1);
rateFnData = zeros([number_samples, 1]);
yData = zeros([number_samples, 1]);

yDataExact = zeros([number_samples, 1]);
rightEigenvalues = zeros([number_samples,D]);
leftEigenvalues = zeros([number_samples,D]);
expSigmaX = zeros([number_samples,1]);
expSigmaZ = zeros([number_samples,1]);
rankrl = zeros([number_samples,2]);
entropy = zeros([number_samples,1]);
dimNullSpaceL = zeros([number_samples,1]);

schmidtCoeffs = zeros([number_samples,D]);

%% Exact solution for rate function
epsilon = @(k, J, g) 2.*J.*sqrt((g - cos(k)).^2 + sin(k).^2);
theta = @(k, g) 0.5 .* atan(sin(k) ./ (g - cos(k)));
phi = @(k, g0, g1) theta(k,g0) - theta(k,g1);
integrand = @(z, k, J, g0, g1) log(cos(phi(k,g0,g1)).^2 + ...
                                sin(phi(k,g0,g1)).^2 .* exp(-2.*z.*epsilon(k,J,g1)));
f = @(g0, g1, J, z) -1*integral(@(k)integrand(z,k,J,g0,g1)./(2*pi),0,pi); 
l = @(g0, g1, J, t) f(g0,g1,J,1i*t) + f(g0,g1,J,-1i*t);
ratefnExact = @(t) l(g0,g,J,t);

%% Extract density matrix and derived quantities at each time

progressWindow = waitbar(0, 'postprocessing...');
for n = 1:number_samples
    % Find A, R, L
    An = output{n};
    [~,Rn,Ln] = normalizeMPS(An);
    
    % find rateFnData(n)
    En0 = ncon({An,conj(A0)},{[-2 -4 1],[-1 -3 1]});
    En0 = reshape(En0,[D^2,D^2]);
    y = eigs(En0, 1);
    rateFnData(n) = - 2*log(abs(y));
    yData(n) = abs(y);
    
    t = T * (n - 1) / number_samples;
    yDataExact(n) = ratefnExact(t);
    
    %%%% Plotting extra stuff as a diagnostic
    rho = ncon({An, Rn, conj(An), Ln}, {[1 2 -1], [2 3], [4 3 -2], [1 4]});
    expSigmaX(n) = real(rho(1,2) + rho(2,1));
    expSigmaZ(n) = real(rho(1,1) - rho(2,2));
    rankrl(n,1) = rank(Rn);
    rankrl(n,2) = rank(Ln);
    evalues = eig(rho);
    S = 0;
    for k = 1:size(evalues,1)
        S = S - evalues(k)*log(evalues(k));
    end
    entropy(n) = S;
    V_L = leftGaugeTangentBasis(Ln, An, 2, D, 8);
    dimNullSpaceL(n) = size(V_L, 2);
    
    leftEigenvalues(n,:) = eig(Ln);
    rightEigenvalues(n,:) = eig(Rn);
    
    % Canonicalize and extract Schmidt coefficients
    CFORM_TOL = 6;
    max_abs_el = abs(max(max(max(An))));
    CFORM_TOL = CFORM_TOL - int64(ceil(log10(max_abs_el)));
    
    mpsn = canonicalForm({1,An,0},CFORM_TOL);
    schmidtCoeffs(n,:) = diag(mpsn{3});
    %%%%
    waitbar(n/number_samples)
end
close(progressWindow)

%% Plot 

fig = figure;
fig.OuterPosition = [672 400 680 580];
times = linspace(0,T,number_samples).';

%% Plotting in terms of Schmidt values

subplot(3,2,[1 2])
hold on
scatter(times,rateFnData,'.');
scatter(times,yData,'.');
if g < g0
    plot(times,yDataExact,'k','LineWidth',1.);
    legend('l(t)','|y(t)|','y_{EXACT}')
else
    legend('l(t)','|y(t)|')
end
xlabel('t')

subplot(3,2,3)
hold on
scatter(times,abs(expSigmaX),'.')
xlabel('t')
ylabel('|<\sigma_x>|')

subplot(3,2,4)
hold on
scatter(times,abs(expSigmaZ),'.')
xlabel('t')
ylabel('|<\sigma_z>|')

subplot(3,2,5)
hold on
scatter(times,real(entropy),'.')
ylabel('von Neumann entropy')
xlabel('t')

subplot(3,2,6)
hold on
for k=1:D
    scatter(times,schmidtCoeffs(:,k),'.')
end
xlabel('t')
ylabel('Schmidt Coefficients')