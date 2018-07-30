%% clear variables we'll be using
clear A analyticalSolution ans denominator dt eta expSigmaZ fieldStrength lyapunov1 times
clear h L n number_samples numerator oneSiteDM output R rho rho_exact state T Udt

%% Initialize state
D = 2; % >= 2

% |psi_0> = |+>
%A = ones([D D 2]);

 %|psi_0> = |0>
  A = zeros([D D 2]);
% A(1,1,1) = 1.;
%   A(:,:,1) = [0.5442 + 0.2168i, 0.3642 + 0.8900i; 0.2809 + 0.2573i, 0.5964 + 0.7352i];
%   A(:,:,2) = [0.8324 + 0.5103i, 0.1323 + 0.0289i; 0.7174 + 0.4643i, 0.2236 + 0.4250i];

 Z = [1, 0; 0, -1];

% randomized intial state
A = rand(D,D,2)+ rand(D,D,2)*1i;
% A(:,:,1)=1
% A(:,:,2)=0
% add some random noise - everything breaks!
%A = A + 0.001*(rand(D,D,2) + rand(D,D,2)*1i);
%A = normalizeMPS(A);

%% construct mps cell array and normalize
[~,L,R,eta] = transferMatrix({1,A,[1 0 0;0 1 0;0 0 1]});

% calculate unnormalized single-site density matrix
numerator = ncon({A, R, conj(A), L},...
                 {[1 2 -1], [2 3], [4 3 -2], [1 4]});
denominator = ncon({A, R, conj(A), L},...
                   {[1 2 5], [2 3], [4 3 5], [1 4]});
oneSiteDM = numerator/denominator;

% calculate normalized tensors
A0 = A;
[A, R, L] = normalizeMPS(A);
state = {1,A,R};

% construct single-site density matrix to check we've only change gauge
rho0 = ncon({A, R, conj(A), L}, {[1 2 -1], [2 4], [3 4 -2], [1 3]});

% single-site transverse field Hamiltonian h = ( 1.Z + Z.1 ) / 2
% fieldStrength = 0.1*pi;
 fieldStrength = 1/2;


h(1,1,1,1) = +1;
h(2,2,2,2) = -1;

h = fieldStrength*h;

% Total time, time step
T = 1;
dt = T/100;

% Sample evolution with tdvpIntegrator
%output = tdvpIntegrator(state,h,T,dt);
output = tdvpIntegrator(A,h,T,dt);

% Extract time evolution and plot against analytical solution
number_samples = size(output,1);

expSigmaX = zeros([number_samples,1]);
expSigmaY = zeros([number_samples,1]);
traceRho = zeros([number_samples,1]);
rankrl = zeros([number_samples,2]);
entropy = zeros([number_samples,1]);
analytical_entropy = zeros([number_samples,1]);
w = zeros([number_samples,2*D^2]);
w2 = zeros([number_samples,2*D^2]);

lyapunov1 = zeros([number_samples,1]);
lyapunov2 = zeros([number_samples,1]);
lyapunov3 = zeros([number_samples,1]);
lyapunov4 = zeros([number_samples,1]);

V=eye(2*D^2);

dimNullSpaceL = zeros([number_samples,1]);

% analytical time-evolution operator factorizes into U(T) = product(U(dt))
Udt = cos(fieldStrength*dt)*[1 0;0 1] - 1i*sin(fieldStrength*dt)*[1 0;0 -1];
analyticalSigmaX = zeros([number_samples,1]);
analyticalSigmaY = zeros([number_samples,1]);
analyticalOutput = cell(number_samples,2);

analyticalOutput{1,1} = A;
analyticalOutput{1,2} = 0;

progressWindow = waitbar(0, 'postprocessing...');
for n = 1:number_samples
    A = output{n};
    [A,R,L,Trans] = normalizeMPS(A);
    A;
    rho = ncon({A, R, conj(A), L}, {[1 2 -1], [2 3], [4 3 -2], [1 4]});

    %%%% Plotting extra stuff as a diagnostic
    traceRho(n) = trace(rho);
    rankrl(n,1) = rank(R);
    rankrl(n,2) = rank(L);
    evalues = eig(rho);
    S = 0;
    for k = 1:size(evalues,1)
        S = S - evalues(k)*log(evalues(k));
    end
    entropy(n) = S;
    
    V_L = leftGaugeTangentBasis(L, A, 2, D, 8);
    dimNullSpaceL(n) = size(V_L, 2);
    %%%%
    
    %This is the bit that calculates Lyapunov exponents
    [Exp_H,Exp_M]=CalculateHExpectation(A, Trans, h, R, L, V_L, D,2);
    
    SUPER_H=[imag(Exp_H)+imag(Exp_M), real(Exp_H)-real(Exp_M); -real(Exp_H)-real(Exp_M), imag(Exp_H)-imag(Exp_M)];
    
    U=expm(dt*SUPER_H)*V;
    [q,r]=qr(U);
    dr=diag(r);
    [~,idx] = sort(abs(dr),'descend');
    dr=dr(idx);
    w(n,:)=diag(r);
    w2(n,:)=dr; 
    % The difference between w and w2 is just whether the lyapunov
    % exponents are ordered so they are smoothly varying (w) or ordered by
    % largested to smallest (w2).
    % You can get any lyapunov exponent by just looking at log(abs(w))
    % or log(abs(w2)).
    
    % Update V for the next timestep.
    V=q;
    
    lyapunov1(n)=log(abs(dr(1)));
    lyapunov2(n)=log(abs(dr(2)));

    expSigmaX(n) = real(rho(1,2) + rho(2,1));
    expSigmaY(n) = real(1i*(rho(1,2) - rho(2,1)));
    if imag(rho(1,1) - rho(2,2)) > 1e-5
        'Warning: nonzero imaginary part of rho11 - rho22'
        n
    end
    %break
    % also find analytical solution
    if n == 1
        rho_exact = rho;
    else    
        rho_exact = Udt*rho_exact*Udt';
        
        analyticalOutput{n,1} = rho_exact;
        analyticalOutput{n,2} = (n-1)*dt;
    end
    analyticalSigmaX(n) = real(rho_exact(1,2) + rho_exact(2,1));
    analyticalSigmaY(n) = real(1i*(rho_exact(1,2) - rho_exact(2,1)));
    
    evalues = eig(rho_exact);
    S = 0;
    for k = 1:size(evalues,1)
        S = S - evalues(k)*log(evalues(k));
    end
    analytical_entropy(n) = S;
   
    waitbar(n/number_samples)
end
close(progressWindow)

%% Line plots
% times = linspace(0,T,number_samples).';
% 
% figure
% subplot(2,2,1)
% plot(times,expSigmaZ,times,analyticalSolution);
% legend('TDVP result', 'Analytical result')
% xlabel('t')
% ylabel('<\sigma_x>')
% 
% subplot(2,2,2)
% plot(times,traceRho);
% xlabel('t')
% ylabel('tr\rho')
% 
% subplot(2,2,3)
% plot(times,entropy)
% xlabel('t')
% ylabel('S')
% 
% subplot(2,2,4)
% plot(times,rankrl(:,1),times,rankrl(:,2))
% legend('rank(r)','rank(l)')
% xlabel('t')

%% Scatter plots

%  ncon({L,A0,Z,conj(A0),R},{[1 2],[1 5 4],[3 4],[2 6 3],[5 6]})
%  ncon({L,A,Z,conj(A),R},{[1 2],[1 5 4],[3 4],[2 6 3],[5 6]})

times = linspace(0,T,number_samples).';

figure

subplot(3,2,1)
hold on
scatter(times,expSigmaX,'.')
scatter(times,analyticalSigmaX,'.');
legend('TDVP result', 'Analytical result')
xlabel('t')
ylabel('<\sigma_x>')

subplot(3,2,2)
hold on
scatter(times,expSigmaY,'.')
scatter(times,analyticalSigmaY,'.');
legend('TDVP result', 'Analytical result')
xlabel('t')
ylabel('<\sigma_y>')

subplot(3,2,3)
scatter(times,(expSigmaX - analyticalSigmaX),'.')
xlabel('t')
ylabel('(<\sigma_x>_{TDVP} - <\sigma_x>_{exact})')

% subplot(3,2,3)
% scatter(times,real(traceRho) - 1,'.');
% xlabel('t')
% ylabel('tr\rho - 1')

subplot(3,2,4)
hold on
scatter(times,real(entropy),'.')
scatter(times,real(analytical_entropy),'.')
legend('TDVP entropy', 'Analytical entropy')
xlabel('t')

subplot(3,2,5)
hold on
scatter(times,lyapunov1,'.')
%ylim([-abs(max(lyapunov1))  abs(max(lyapunov1))])
xlabel('t')
ylabel('Lyapunov_1')

subplot(3,2,6)
hold on
scatter(times,lyapunov2,'.')
%ylim([-abs(max(lyapunov1))  abs(max(lyapunov1))])
xlabel('t')
ylabel('Lyapunov_2')
%dim = [.2 .5 .3 .3];
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
% subplot(3,2,6)
% hold on
% scatter(times,dimNullSpaceL,'.')
% xlabel('t')
% ylabel('# null vectors of L')
% %dim = [.2 .5 .3 .3];
% title(['(d = 2, D = ',num2str(D),')']);
% %annotation('textbox',dim,'String',str,'FitBoxToText','on');