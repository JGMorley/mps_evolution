%% SUBSIDIARY SCRIPT FOR isingQuench_findGSfirst.m

% Clear variables
clear A0 T expSigmaZ rankrl An V_L g rateFnData d h rho En0 dimNullSpaceL k             
clear times J dt minDim y Ln entropy n yData Rn evalues number_samples S  
clear expSigmaX output      

%% Setup initial state

% randomized MPS tensor
A0(:,:,1) = rand(D) + rand(D)*1i; 
A0(:,:,2) = rand(D) + rand(D)*1i; 
A0 = normalizeMPS(A0);

%% Print initial single-site density matrix
[A_init, R_init, L_init] = normalizeMPS(A0);
rho_init = ncon({A_init, R_init, conj(A_init), L_init},...
                {[1 2 -1], [2 4], [3 4 -2], [1 3]})

%% Hamiltonian tensor
J = 1.;
g = g0;

h(:,1,:,1) = [-J     -J*g/2;-J*g/2     +J];
h(:,1,:,2) = [-J*g/2      0;     0 -J*g/2];
h(:,2,:,1) = [-J*g/2      0;     0 -J*g/2];
h(:,2,:,2) = [+J     -J*g/2;-J*g/2     -J];

%% Run tdvp
T = 3; 
dt = T/100; % Total time, time step
output = tdvpIntegrator(A0,h,T,dt,'IMAG_TIME',true,'INVERSE_FREE',true,'SCHMIDT_TH',1e-10,'CFORM_TOL',8);
disp('MPS evolved ok')

%% Initialize arrays to be plotted
number_samples = size(output,1);
rateFnData = zeros([number_samples, 1]);
yData = zeros([number_samples, 1]);

rightEigenvalues = zeros([number_samples,D]);
leftEigenvalues = zeros([number_samples,D]);
expSigmaX = zeros([number_samples,1]);
expSigmaZ = zeros([number_samples,1]);
rankrl = zeros([number_samples,2]);
entropy = zeros([number_samples,1]);
dimNullSpaceL = zeros([number_samples,1]);
energy = zeros([number_samples,1]);

schmidtCoeffs = zeros([number_samples,D]);

%% Extract density matrix and derived quantities at each time

warning('off','ncon:suboptimalsequence') % turn off this warning

progressWindow = waitbar(0, 'postprocessing...');
for n = 1:number_samples
    % Find A, R, L
    An = output{n};
    [~,Rn,Ln] = normalizeMPS(An);
    
    % find rateFnData(n)
    En0 = ncon({An,conj(A0)},{[-2 -4 1],[-1 -3 1]});
    En0 = reshape(En0,[D^2,D^2]);
    y = eigs(En0, 1);
    rateFnData(n) = - 2 * log(abs(y));
    yData(n) = abs(y);
    
    %%%% Plotting extra stuff as a diagnostic
    rho = ncon({An, Rn, conj(An), Ln}, {[1 2 -1], [2 3], [4 3 -2], [1 4]});
    energy(n) = ncon({An,An,Rn,conj(An),conj(An),Ln,h},...
                     {[1 2 5],[2 4 3],[4 6],[8 6 7],[9 8 10],[1 9],[10 7 5 3]});
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
     leftEigenvalues(n,:) = eig(Ln);
     rightEigenvalues(n,:) = eig(Rn);
    waitbar(n/number_samples)
end
close(progressWindow)

%% Plot 

fig = figure;
fig.OuterPosition = [672 400 680 580];
times = linspace(0,T,number_samples).';

%% Plotting in terms of right, left environments

subplot(4,2,[1 2])
hold on
scatter(times,real(rateFnData),'.');
scatter(times,real(yData),'.');
xlabel('t')
legend('l(t)','|y(t)|')

subplot(4,2,3)
hold on
scatter(times,real(expSigmaZ),'.')
xlabel('t')
ylabel('<\sigma_z>')

subplot(4,2,4)
hold on
scatter(times,abs(expSigmaX),'.')
xlabel('t')
ylabel('|<\sigma_x>|')

subplot(4,2,5)
hold on
scatter(times,real(entropy),'.')
ylabel('von Neumann entropy')
xlabel('t')

subplot(4,2,6)
hold on
scatter(times, real(energy),'.')
ylabel('energy density')
xlabel('t')

subplot(4,2,7)
hold on
left_legends = {};
for k=1:D
    maxval = 1.;%max(leftEigenvalues(:,k));
    scatter(times,real(leftEigenvalues(:,k)/maxval),'.')
    %left_legends{k} = ['max e_value = ',num2str(round(maxval,4))];
end
xlabel('t')
ylabel('Norm-zd left e/values')
%legend(left_legends)

subplot(4,2,8)
hold on
right_legends = {};
for k=1:D
    maxval = 1.;%max(rightEigenvalues(:,k));
    %right_legends{k} = ['max e_value = ',num2str(round(maxval,4))];
    %maxval = 1.;
    scatter(times,real(rightEigenvalues(:,k)/maxval),'.')
    
end
xlabel('t')
ylabel('Norm-zd right e/values')
legend(right_legends)