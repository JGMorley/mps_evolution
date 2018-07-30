function [H1,H2,H3, GAMMA] = CalculateLyapunov(A, H, right, left, VL, chi, d)
% This function outputs four (chi^2, chi^2) matrices 
% H1 is basically <d_i psi] H [ d_j psi>
% H2 is <d_i d_j psi] H [ psi>
% H3 is the normalization terms from <d_i d_j psi] H [ psi>
% when i and j are on the same site. 
% GAMMA is the term that deals with the parallel transport.
% GAMMA looks like g^km <d_i d_j psi] H [ d_m psi> d(A_k)/dt

% Construct the transfer operator E.
E = ncon({A,conj(A)},{[-2 -4 1],[-1 -3 1]}); 
E=reshape(E,chi^2,chi^2);

% Find the inverse of the transfer operator. Defined such that the left and
% right environment are zero eigenstates of it. 
Q=eye(chi^2)-reshape(right.',[chi^2,1])*reshape(left.',[1,chi^2]);
TransN=reshape(Q*inv(eye(chi^2)-Q*E*Q)*Q,chi,chi,chi,chi);

% Define "B(i)" which spans the tangent space of A. The chi^2 elements of B
% are the chi^2 basis elements of the tangent space.

basis=zeros(chi,chi);

B=zeros(chi,chi,d,chi^2);
for i=1:numel(basis)
    basis=zeros(chi,chi);
    basis(i)=1;
    B(:,:,:,i)=  ncon({sqrtm(inv(left)), VL, basis, sqrtm(inv(right))},{[1 -1],[1 2 -3],[2 3 -4],[3 -2]});
end

% Test they are orthonormal
% test=zeros(numel(basis),numel(basis));
% for i=1:numel(basis)
%      for j=1:numel(basis)
%         test(i,j)=round(ncon({left, B(:,:,:,i), conj(B(:,:,:,j)),right},{[1 2],[1 4 3],[2 5 3],[4 5]}),10);
%      end
%  end
% test;
 

% A2 is used in H3, when you have d_id_j and i=j then instead of have 
% V_L * x  you have something like A * x^dagger x 
basis1=zeros(chi,chi);
basis2=zeros(chi,chi);
A2=zeros(chi,chi,d,chi^2,chi^2);

for i=1:numel(basis1)
    for j=1:numel(basis2)
    basis1=zeros(chi,chi);
    basis1(i)=1;
    basis2=zeros(chi,chi);
    basis2(i)=1;

    A2(:,:,:,i,j)=  (-1/2)*ncon({A,sqrtm(inv(right)),basis1',basis2,sqrtm(inv(right))},{[-1 1 -3],[1 2],[2 3],[3 4],[4 -2]});
    end
end



%Do some contractions which will be used to find all the elements of H and

AA=ncon({A,A},{[-1 1 -3], [1 -2 -4]});
AB=ncon({A,B},{[-1 1 -3], [1 -2 -4 -5]});
BA=ncon({B,A},{[-1 1 -3 -5], [1 -2 -4]});
BB=ncon({B,B},{[-1 1 -3 -5], [1 -2 -4 -6]});

A2A=ncon({A2,A},{[-1 1 -3 -5 -6], [1 -2 -4]});
AA2=ncon({A,A2},{[-1 1 -3], [1 -2 -4  -5 -6]});



TransBB=ncon({TransN,B,conj(B),right},{[-2 -1 5 4],[4 1 3 -3],[5 2 3 -4],[1 2]});

TransA2A=ncon({TransN,A2,conj(A),right},{[-2 -1 5 4],[4 1 3 -3 -4],[5 2 3],[1 2]});


TransAB=ncon({TransN,A,conj(B),right},{[-2 -1 5 4],[4 1 3],[5 2 3 -3],[1 2]});
TransBA=ncon({TransN,B,conj(A),right},{[-2 -1 5 4],[4 1 3 -3],[5 2 3 ],[1 2]});

TransABBA=ncon({TransN,A,conj(B),TransBA},{[-2 -1 5 4],[4 1 3],[5 2 3 -4],[1 2 -3]});
TransABAB=ncon({TransN,A,conj(B),TransAB},{[-2 -1 5 4],[4 1 3],[5 2 3 -3],[1 2 -4]});

% BX is used in the construction of GAMMA, the parallel transport term
% below. I have done this find_F thing in the sloppiest way possible, must
% fix. 

% BX is just VL * x where x is determined using A(t) in the standard TDVP
% method.
x=find_F(A,1,H);
BX= ncon({sqrtm(inv(left)), VL, x, sqrtm(inv(right))},{[1 -1],[1 2 -3],[2 3],[3 -2]});

% Now we do all of the possible contractions to find H1, H2, H3 etc. 
% Some diagrams included in H1 also have their complex conjugate in H1, in
% those cases I have just taken the conjugate transpose. 

% Similarly, there are cases where H2 has diagrams and there there (NOT
% conjugate) transpose, in those cases I have just taken the transpose.
H1=zeros(chi^2,chi^2);
H2=zeros(chi^2,chi^2);
H3=zeros(chi^2,chi^2);
GAMMA=zeros(chi^2,chi^2);

        obv1B=ncon({left,AB,H,conj(AB),right},{[1 2],[1 3 7 8 -1],[5 6 7 8],[2 4 5 6 -2],[3 4]});
        obv1=obv1B;
        
        obv2=ncon({left,AB,H,conj(BA),right},{[1 2],[1 3 7 8 -1],[5 6 7 8],[2 4 5 6 -2],[3 4]});

        obv3=obv2';

        obv4B=ncon({left,BA,H,conj(BA),right},{[1 2],[1 3 7 8 -1],[5 6 7 8],[2 4 5 6 -2],[3 4]});
        obv4=obv4B;
        
        obv6a=ncon({left,AA,H,conj(BB),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6 -1 -2],[3 4]});
        obv6b=obv6a.';
        obv6=(obv6a+obv6b);
        
        obv7B=ncon({left,AA,H,conj(AA),TransBB},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4 -1 -2]});
        obv7=obv7B;
        

        obv9=ncon({left,BA,H,conj(AA),TransAB},{[1 2],[1 3 7 8 -1],[5 6 7 8],[2 4 5 6],[3 4 -2]});

        obv10=ncon({left,AB,H,conj(AA),TransAB},{[1 2],[1 3 7 8 -1],[5 6 7 8],[2 4 5 6],[3 4 -2]});

        obv11=obv9';

        obv12=obv10';

        obv15a=ncon({left,AA,H,conj(BA),TransAB},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6 -2],[3 4 -1]});
        obv15b=obv15a.';
        obv15=(obv15a+obv15b);
        
        obv16a=ncon({left,AA,H,conj(AB),TransAB},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6 -2],[3 4 -1]});
        obv16b=obv16a.';
        obv16=(obv16a+obv16b);
        
        obv17=ncon({left,AA,H,conj(AA),TransABBA},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4 -1 -2]});
        
        obv18=obv17';

        obv20a=ncon({left,AA,H,conj(AA),TransABAB},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4 -1 -2]});
        obv20b=obv20a.';
        obv20=(obv20a+obv20b);

             
        obv21a=ncon({left,AA,H,conj(AA),TransA2A},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4 -2 -1]});
        obv21b=obv21a';
        obv21=obv21a+obv21b;

        obv22a=ncon({left,A2A,H,conj(AA),right},{[1 2],[1 3 7 8 -2 -1],[5 6 7 8],[2 4 5 6],[3 4]});      
        obv22b=obv22a';

        obv22=obv22a+obv22b;

        
        obv23a=ncon({left,AA2,H,conj(AA),right},{[1 2],[1 3 7 8 -2 -1],[5 6 7 8],[2 4 5 6],[3 4]});      
        obv23b=obv23a';
        obv23=obv23a+obv23b;
        
        gam1= ncon({left,BX,conj(B),TransAB},{[1 2],[1 4 3],[2 5 3 -1],[4 5 -2]});
        gam2= gam1.';

        H1=obv1+obv2+obv3+obv4+obv7+obv9+obv10+obv11+obv12+obv17+obv18;%+obv6+obv14+obv15+obv19;

        H2= obv6+obv15+obv16+obv20;
        
        H3=obv21+obv22+obv23;
        
        GAMMA=gam1+gam2;
        


% M2-GAMMA
end


function F = findF( A, t, H, small_change_flag, KWARGS )
%FINDTIMEDERIVATIVE Calculate f:=(d/dt)A(A, t)
    
    [ ~, right, left, E ] = normalizeMPS(A);

    % Calculate V_L, the matrix of null vectors of lambda*A (L)
    if strcmp(KWARGS.Vmethod,'nullspace')
        sqrtleft = squareDecomposition(left, KWARGS.SQRT_TOL, 2);
        [~,D,d] = size(A);
        V_L = leftGaugeTangentBasis(sqrtleft, A, d, D, KWARGS.VL_TOL);            
    else
        if t==0
            warning('Using unitary gauge, requires left canonical form');
            'This method is giving bad results! Need to fix!'
        end
        [~,D,d] = size(A);

        [A, right, U] = tangentBasis_unitaryGauge(A, right);
        U = permute(U, [1 2 4 3]); % U[i j rho sigma]
        V_L = reshape(U(:,:,2:end,:),[D D*(d-1) d]); % U[i (j,rho!=1), sigma]
    end


    % Compute C
    % NB for now considering one two-site operator h
    C = calculateC(A,H);

    % Calculate K
    K = calculateK(C,A,E,left,KWARGS.K_TOL);

    % Find inverses of l^(1/2), r^(1/2)
    lm12 = inv( squareDecomposition(left, KWARGS.SQRT_TOL, 2) );
    rm12 = inv( squareDecomposition(right,KWARGS.SQRT_TOL, 2) );

    % Calculate F
    F = calculateF(V_L,left,C,right,A,K,lm12,rm12,KWARGS.SQRT_TOL);

    
end

function C = calculateC(A,h)
    %% Calculate C, a [DxDxpxp] tensor
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % A: [DxDxp] tensor
    %      MPS tensor
    % h: [pxpxpxp] tensor
    %      Hamiltonian tensor: h(i,j,k,l) = <i|<j|\hat{h}|k>|l>
    %% Diagram
    %          _____
    %  [-1] --|     |-- [-2]     [-1] -- A --- A -- [-2]
    %         |  C  |                    |_____|   
    %         |_____|         =         |       | 
    %          |   |                    |___H___|
    %        [-3] [-4]                   |     |   
    %%  
    
    % check dimensions
    p = size(A,3);
    if ~isequal(size(h), [p p p p])
        err.message = 'Hamiltonian tensor h has wrong dimension';
        err.identifier = 'calculateC:hWrongDimension';
        error(err);
    end
    
    % compute C
    C = ncon({A,A,h},{[-1 1 2], [1 -2 3], [-3 -4 2 3]});
    % (This entire fn could be replaced with one line, however having it as a 
    % separate function allows unit tests to be written for it)
end

function K=calculateK(C,A,E,left,K_TOL)
    %% Calculate K, a [DxD] matrix
    % As shown in e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % C: [DxDxpxp] tensor
    %      Two-site Hamiltonian tensor contracted with two A-tensors
    % A: [DxDxp] tensor
    %      MPS tensor normalized st dominant eigenvalue of transfer matrix = 1
    % E: [D^2xD^2] matrix
    %      Transfer matrix
    % left: [DxD] matrix
    %      Dominant left eigenvector of the transfer matrix E
    %
    %% Diagram:
    %
    %  /  |--      /  |--|     |--|     |--
    % |   |       |   |  |  C  |  |     |
    % | K |     = | l |  |_____|  | E_P |  
    % |   |       |   |   |   |   |     |
    %  \  |--      \  |-- A - A --|_____|--
    %% Another method (using iterative solver): - will be good for checking!
    % K = bicgstab(M,b) solves M*K=b for K
    % M = diag(ones(1,D)) - E + right*left
    % b = (ncon({C,A',A',left},{[1 -1 2 3],[4 5 2],[5 -2 3],[1 4]})) + ...
    % ((ncon({C,A',A',left,right},{[1 2 3 4],[5 6 3],[6 7 4],[1 5],[2 7]}))*left)
    %%
    
    if nargin==4
        K_TOL = 10;
    end
    
    % Find right and left eigenvectors and their associated eigenvalues
    D = sqrt(size(E,1));
    VR = zeros(D^2);
    VL = zeros(D^2);
    idx = 0; % degeneracies sometimes aren't normalized, so try up to 10 times
             % or until VR, VL are full-rank
    while ( rank(VR) < D^2  || rank(VL) < D^2 ) && idx < 10
        [VR,eigenvalues,VL] = eig(E); % cols of VL are e-vectors of E' (not E.')
        idx = idx + 1;
    end
    eigenvalues = diag(eigenvalues);
    maxEValue = max(eigenvalues);
    
    % If E has degeneracies eig() sometimes gives degenerate eigenvectors
    if ( rank(VR) < D^2 ) || ( rank(VL) < D^2 )
        err.identifier = 'TDVP:calculateK:degenerateEigenvectors';
        err.message =['Unable to find ',num2str(D^2),...
                      ' orthogonal eigenvectors after 10 iterations'];
        error(err)
    end
    
    % eig() gives us VL(:,i)'*VR(:,j!=i) = 0, however VL(:,i)'*VR(:,i) != 1,
    % so we need to normalize (but don't need to orthogonalize)
    
    for k=1:D^2
        norm = sqrt(VL(:,k)'*VR(:,k)); % norm is in general complex...
        VR(:,k) = VR(:,k)/norm;
        VL(:,k) = VL(:,k)/conj(norm);   % ...so divide by conjugate here
    end
    
    % check that sum_i |r_i)(l_i| = I
    overlap = @(i,j) VL(:,i)'*VR(:,j);
    iden = overlap(1:D^2, 1:D^2);
    if ~isequal(round(iden,K_TOL), eye(D^2))
        err = abs(max(max(iden - eye(D^2))));
        warning('sum_i |r_i)(l_i| = I + O(%.1d) != I',err);        
        if abs(err) > 1e-3
            [sprintf('iden = '),mat2str(iden),...
             sprintf('\nE = '),mat2str(E),...
             sprintf('\neigenvalues = '),mat2str(eigenvalues),...
             sprintf('\nVL = '),mat2str(VL),...
             sprintf('\nVR = '),mat2str(VR)]
        end
    end
                                             
    % E_P = sum_i [1/(1-lambda_i)]|ri)(li|
    E_P = zeros(D^2);
    for i=1:D^2
        if eigenvalues(i) ~= maxEValue
            coeff = 1 / (1 - eigenvalues(i));
            E_P=E_P + coeff * VR(:,i) * VL(:,i)'; 
            % NB we're using (v| = |v)' != |v).' here
        end
    end

    % reshape into DxDxDxD tensor
    E_P = permute(reshape(E_P,[D,D,D,D]), [2 4 1 3]);
               
    K = ncon({left, C, conj(A), conj(A), E_P},...
             {[1 2], [1 6 3 5], [2 4 3], [4 7 5], [6 -1 7 -2]});
end

function F = calculateF(V_L, left, C, right, A, K, lm12, rm12, TOL)
    %% Calculate F, a [(p-1)DxD] tensor
    % See e.g. 'Summary' in Haegeman et al, PRB 88, 075133 (2013)
    %
    %% Inputs
    % V_L:   [Dx(p-1)D] tensor
    %          columns of V_L form a basis over null space of A'*lambda tensor
    % left:  [DxD] matrix
    %          Dominant left eigenvector of the transfer matrix E
    % C:     [DxDxpxp] tensor
    %          Two-site Hamiltonian tensor contracted with two A-tensors
    % right: [DxD] matrix
    %          Dominant right eigenvector of the transfer matrix E
    % A:     [DxDxp] tensor
    %          MPS tensor
    % K:     [DxD] matrix
    %          tensor arising from sum of terms with B on RHS of Hamiltonian
    % lm12:  [DxD] matrix
    %          left^(-1/2) calculated previously
    % rm12:  [DxD] matrix
    %          right^(-1/2) calculated previously
    %%

    %% Compute x^(-1/2) and x^(1/2) for x = left, right
    lm12c = conj(lm12);
    rm12c = conj(rm12);

    l12 = squareDecomposition(left, TOL, 2);
    r12 = squareDecomposition(right, TOL, 2);
    
    %% Contract the three terms of F and sum
    
    % F1 = (l|H(AA,BA)|r)
    F1 = ncon({conj(V_L),l12,C,right,conj(A),rm12c},...
              {[2 -1 3], [1 2], [1 7 3 5], [7 6], [4 6 5], [-2 4]});
          
    % F2 = (l|H(AA,AB)|r)
    F2 = ncon({conj(V_L),lm12c.',conj(A),left,C,r12},...
              {[1 -1 2], [3 1], [5 3 4], [6 5], [6 7 4 2], [7 -2]});
 
    % F3 = sum_{k=0}^{inf} (l|H(AA,AA) E(A,A)^k E(A,B)|r)
    F3 = ncon({conj(V_L),lm12c.',K,A,r12},...
              {[1 -1 2], [3 1], [4 3], [4 5 2], [5 -2]});

    F = F1 + F2 + F3;
end

function F=find_F(A,t,H, varargin)

    %% Parse input arguments
    
    p = inputParser;
    
    defaultVmethod = 'nullspace';
    validVmethods = {'nullspace','unitarygauge'};
    checkVmethod = @(x) any(validatestring(x,validVmethods));
    
    defaultIMAG_TIME = 0;
    
    defaultSQRT_TOL    = 6;    % # digits tol when taking square root
    defaultRANK_TOL    = 14;   % # digits tol when checking rank of left, right
    defaultK_TOL       = 8;    % # digits tol for checks in calculateK
    defaultVL_TOL      = 6;    % # digits tol for checks in leftGaugeTangentBasis
    defaultBF_TOL      = 6;    % # digits tol for checks in tangentVector
    defaultMAX_EL_TOL  = 1e4;  % Largest el of right, left without error
    defaultNOISE_MAG   = 1e-4; % magnitude of noise added to MPS if chi < D
    defaultCFORM_TOL   = 12;   % # digits tol in checks within canonicalForm
    defaultN_LOOPS     = 10;  % # times to try adding noise to get full rank
     
    addOptional(p,'Vmethod',defaultVmethod,checkVmethod)
    addOptional(p,'IMAG_TIME',defaultIMAG_TIME,@islogical)
    addParameter(p,'SQRT_TOL',defaultSQRT_TOL,@isnumeric)
    addParameter(p,'RANK_TOL',defaultRANK_TOL,@isnumeric)
    addParameter(p,'K_TOL',defaultK_TOL,@isnumeric)
    addParameter(p,'VL_TOL',defaultVL_TOL,@isnumeric)
    addParameter(p,'BF_TOL',defaultBF_TOL,@isnumeric)
    addParameter(p,'MAX_EL_TOL',defaultMAX_EL_TOL,@isnumeric)
    addParameter(p,'NOISE_MAG',defaultNOISE_MAG,@isnumeric)
    addParameter(p,'CFORM_TOL',defaultCFORM_TOL,@isnumeric)
    addParameter(p,'N_LOOPS',defaultN_LOOPS,@isnumeric)

    parse(p,varargin{:})
    small_change_flag = 0;
    
F=findF( A, t, H, small_change_flag, p.Results );
end
