function samples = tdvpIntegrator(A0, H, T, dt, varargin)
    %% Evolve forwards the TDVP equations for MPS. With extra bits
    % Features (in varying combinations) :
    % > RK4
    % > Symmetric timestep
    % > Grow bond dimension
    % > Inverse-free timestep
    % (also imaginary time, but this is trivially implemented with dt->1i*dt)
    %
    % Preferences and tolerances should be defined in varargin - see
    % function 'parse_tdvpinputs' for details
    %%
    
    %% 1. Run setup function: parse varargin and make sampling variables
    [KWARGS, D, d, A, right, z, samples, iSample, time_between_samples] = ...
    tdvpIntegratorSetup(A0, T, dt, varargin{:});

    %% 2. Integration loop
    progressWindow = waitbar(0, 'Time loop 0.00% complete.');
    for t = dt:dt:T
               
        % Find z for first inverse free iteration
        if KWARGS.INVERSE_FREE && t==dt
            [A,r] = normalizeMPS(A);
            [U,S] = svd(r);
            z = U*sqrt(S);
        end
        
        % Find dA
        if KWARGS.INVERSE_FREE && abs(min(min(z*z'))) < KWARGS.SCHMIDT_TH

            [AR,z]=find_AR_z_from_AL(A,z,d,D);
            
            if KWARGS.GROW_BOND_DIM
                if D < KWARGS.DMAX
                    [A,AR,z,D]=increase_bond_dim(A,'none','none',H,D,KWARGS.DMAX,d);
                end
            end   
            
            update_method = 'INV FREE';
            dA = inverseFreeTimestep( A, AR, z, dt, H, d, D, KWARGS );
        
        else
            update_method = 'INV USED';
            dA = inverseUsedTimestep( A, t==0, dt, H, KWARGS );
            if KWARGS.IMAG_TIME, dA = -1i*dA; end
        end
 
        A = A + dA;
        
        if ~KWARGS.INVERSE_FREE
            [~, A] = canonicalForm({1,A,right},KWARGS.CFORM_TOL);
        end
        
        % Output A(t) if next sample time reached
        if t + KWARGS.BUFFER >= (iSample - 1) * time_between_samples
            samples{iSample,1} = A; 
            samples{iSample,2} = t;
            iSample = iSample + 1;
        end
        
        waitbar_msg = ['Time loop ',num2str(100*t/T,'%.2f'),'% complete. ',...
                       'Update method is ',update_method];
        waitbar(t/T, progressWindow, waitbar_msg);
    end
    close(progressWindow)
end


%% setup function - keeping the messy bits separate

function [KWARGS, D, d, A, R, z, samples, iSample, t_b_s] = ...
         tdvpIntegratorSetup(A0, T, dt, varargin)
    % 1. Parse optional inputs and set defaults
    KWARGS = parse_tdvpinputs(T, dt, varargin);
    
    if KWARGS.GROW_BOND_DIM && strcmpi(KWARGS.DMAX,'none')
        % no bond dimension supplied but growth asked for
        warning('No DMAX specified, so running without bond dimension growth')
        KWARGS.GROW_BOND_DIM = false;
    end
    % 2. Gauge transform input state and check full rank within INV_TOL
    [~,D,d] = size(A0);
    if KWARGS.INVERSE_FREE
       [A, R] = normalizeMPS(A0);
    else
        mps = canonicalForm({1,A0,0},KWARGS.CFORM_TOL);
        [~,A,R] = mps{:};
    end
    [U,S] = svd(R);
    z = U*sqrt(S); 
    
    smallest_Schmidt = min(diag(R));
    if (smallest_Schmidt < 10^-(KWARGS.RANK_TOL)) && ~(KWARGS.INVERSE_FREE)
        error('Environment not full rank. Try using inverse-free integrator')
    end
    
    % 3. Intialize time-loop and sampling variables    
    if KWARGS.nSAMPLES > (T / dt)
        warning(['Number of samples larger than number of timesteps - ',...
                 'sampling at every timestep of ',num2str(dt)]);
        KWARGS.nSAMPLES = int64(T / dt);
    end
    t_b_s = T / double(KWARGS.nSAMPLES);
    
    samples = cell(KWARGS.nSAMPLES + 1, 2); 
    samples{1,1} = A; % initial state as first element - not counted as a sample
    samples{1,2} = 0; % initial time = 
    iSample = 2; % index for first sample from the integrator
end

%% parsing function

function KWARGS = parse_tdvpinputs( T, dt, varargin )
    %% PARSE_TDVPINPUTS Parse the keyword arguments for tdvpIntegrator
    %
    %% Compatability table 26/09/17
    %
    %      Feature   | Valid combinations (not finished yet!)
    %  --------------|-----------------------------------
    %  SYM_STEP      | N N N N Y N N
    %  RK4           | N N N Y N N Y
    %  INVERSE_FREE  | N N Y N N Y Y
    %  GROW_BOND_DIM | N Y N N N Y N
    
    %% Parse
    p = inputParser;
    
    defaultVmethod = 'nullspace';
    validVmethods = {'nullspace','unitarygauge'};
    checkVmethod = @(x) any(validatestring(x,validVmethods));
    
    checkDMAX = @(x) any([strcmpi(x,'none'),isnumeric(x)]);
    
    defaultIMAG_TIME = 0;
    
    defaultSCHMIDT_TH    = 1e-10;          % min schmidt val for inverse used  
    defaultBUFFER        = 100*eps;         % time buffer for sampling
    defaultnSamples      = int64(T/dt); % number of samples
    defaultSQRT_TOL      = 6;     % # digits tol when taking square root
    defaultRANK_TOL      = 14;    % # digits tol when checking rank of left, right
    defaultK_TOL         = 8;     % # digits tol for checks in calculateK
    defaultVL_TOL        = 6;     % # digits tol for checks in leftGaugeTangentBasis
    defaultBF_TOL        = 6;     % # digits tol for checks in tangentVector
    defaultMAX_EL_TOL    = 1e4;   % Largest el of right, left without error
    defaultNOISE_MAG     = 1e-4;  % magnitude of noise added to MPS if chi < D
    defaultCFORM_TOL     = 12;    % # digits tol in checks within canonicalForm
    defaultN_LOOPS       = 10;    % # times to try adding noise to get full rank
    defaultDMAX          = 8; % max bond dimension to grow to 
    defaultSYM_STEP      = false; % Take symmetric time step
    defaultRK4           = false; % Use RK4 integrator where appropriate
    defaultINVERSE_FREE  = false; % Use inverse-free integrator where appropriate
    defaultGROW_BOND_DIM = false;  % Grow the bond dimension up to a user-defined limit
     
    addOptional(p,'Vmethod',defaultVmethod,checkVmethod)
    addOptional(p,'IMAG_TIME',defaultIMAG_TIME,@islogical)
    addOptional(p,'GROW_BOND_DIM',defaultGROW_BOND_DIM,@islogical)
    
    addOptional(p,'SCHMIDT_TH',defaultSCHMIDT_TH,@isnumeric)
    addOptional(p,'DMAX',defaultDMAX,checkDMAX)
    addOptional(p,'RK4',defaultRK4,@islogical)
    addOptional(p,'INVERSE_FREE',defaultINVERSE_FREE,@islogical)
    addOptional(p,'SYM_STEP',defaultSYM_STEP,@islogical)
    addParameter(p,'nSAMPLES',defaultnSamples,@isnumeric)
    addParameter(p,'BUFFER',defaultBUFFER,@isnumeric)
    addParameter(p,'SQRT_TOL',defaultSQRT_TOL,@isnumeric)
    addParameter(p,'RANK_TOL',defaultRANK_TOL,@isnumeric)
    addParameter(p,'K_TOL',defaultK_TOL,@isnumeric)
    addParameter(p,'VL_TOL',defaultVL_TOL,@isnumeric)
    addParameter(p,'BF_TOL',defaultBF_TOL,@isnumeric)
    addParameter(p,'MAX_EL_TOL',defaultMAX_EL_TOL,@isnumeric)
    addParameter(p,'NOISE_MAG',defaultNOISE_MAG,@isnumeric)
    addParameter(p,'CFORM_TOL',defaultCFORM_TOL,@isnumeric)
    addParameter(p,'N_LOOPS',defaultN_LOOPS,@isnumeric)

    parse(p,varargin{1}{:})
    
    %% Assign values to KWARGS strut object
    KWARGS = p.Results;
    KWARGS.nSAMPLES = int64(KWARGS.nSAMPLES); % cast nSamples to integer
    
    %% Raise warnings and correct for incompatible inputs
    combo = [KWARGS.SYM_STEP,KWARGS.RK4,KWARGS.INVERSE_FREE,KWARGS.GROW_BOND_DIM];
    valid_feature_combos = logical([0 0 0 0;1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1;...
                                    0 0 1 1;0 1 1 0]);
    combo_is_valid = false;
    for k=1:length(valid_feature_combos)
        if isequal(combo, valid_feature_combos(k,:))
            combo_is_valid = true;
        end
    end
    if ~combo_is_valid
        error('Invalid combination of feature inputs')
    end
    
    if isequal(combo,[0 1 1 0])
        warning('Inverse-free RK4 not working properly yet, so using Euler step')
        KWARGS.RK4 = false;
    end
end