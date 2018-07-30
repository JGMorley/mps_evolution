%% Main function to generate tests
% run by executing 'run(unitTesting);' in the Command Window
function tests = tdvpTesting
    tests = functiontests(localfunctions);
end

%% Test functions

function test_normalizeMPS_randomInput(testCase)
    % given random input, check outputs have correct properties
    A = rand([3,3,2]) + rand([3,3,2])*1i;
    [A,R,L,actE] = normalizeMPS(A);
    
    actOverlap = trace(L.'*R); % should equal 1
    
    verifyEqual(testCase, actOverlap, 1, 'AbsTol', 1e-10);
    
    E = ncon({A,conj(A)},{[-2 -4 1],[-1 -3 1]});
    E = reshape(E,[9,9]);
    
    verifyEqual(testCase, actE, E, 'AbsTol', 1e-10);
    
    actDomEValue = eigs(E,1); % should equal 1
    
    verifyEqual(testCase, actDomEValue, 1., 'AbsTol', 1e-10);
    
    E = permute(reshape(E, [3 3 3 3]), [2 1 4 3]);
    LEminusL = ncon({L, E}, {[1 2], [1 2 -1 -2]}) - L;
    ERminusR = ncon({E, R}, {[-1 -2 1 2], [1 2]}) - R;
   
    verifyEqual(testCase, LEminusL, zeros(3), 'AbsTol', 1e-10);
    verifyEqual(testCase, ERminusR, zeros(3), 'AbsTol', 1e-10);
end

function test_leftGaugeTangentBasis_fixture(testCase)
    % check we get null vectors for the fixture case
    A = testCase.TestData.state1{2};
    left = [sqrt(2) 0; 0 sqrt(2)];
    
    VL = leftGaugeTangentBasis(left, A);
    
    verifyEqual(testCase, size(VL), [2 4 3]);
    
    p = 3;
    D = 2;
    
    L = ncon({conj(A),left}, {[1 -3 -2], [-1 1]});  %L_{a, i, b}  
    L = permute(L, [3 2 1]); % L_{b, i, a}
    L = reshape(L, [D,p*D]); % L_{b,(a,i)}
    
    for b = 1:D        % indexes elements
        for k = 1:D*(p-1) % indexes null vectors
            % check b'th component of L*[kth null vector] is zero
            el = L(b,:)*reshape(permute(VL(:,k,:), [3 2 1]),[p*D,1]);
            verifyEqual(testCase, el, 0., 'AbsTol', 1e-13);
        end
    end
end

function test_leftGaugeTangentBasis_errors(testCase)
    % check errors messages are called ok
    verifyError(testCase, @()leftGaugeTangentBasis(0),...
                'leftGaugeTangentBasis:WrongNoInputs'); 
end

function test_transferMatrix_fixtureTest(testCase)
    % test outputs against previously verified fixture outputs
    [E, VL, VR, eValue] = transferMatrix(testCase.TestData.state1);

    verifyEqual(testCase, E, ...
     [1.0000 - 0.0000i  0.0001 - 0.0000i  0.0001 + 0.0000i  0.9988 + 0.0000i;...
     -0.0000 - 0.0000i -0.0212 + 0.0180i  0.0012 - 0.0046i  0.0189 + 0.0029i;...
     -0.0000 + 0.0000i  0.0012 + 0.0046i -0.0212 - 0.0180i  0.0189 - 0.0029i;...
      0.0001 + 0.0000i -0.0001 + 0.0001i -0.0001 - 0.0001i  0.0013 + 0.0000i],...
     'AbsTol', 1e-4);
 
    verifyEqual(testCase, VL, ...
     [ 0.7071 + 0.0000i  0.0000 - 0.0000i;...
       0.0000 + 0.0000i  0.7071 - 0.0000i],...
     'AbsTol', 1e-4);
  
    verifyEqual(testCase, VR, ...
     [ 1.0000 + 0.0000i -0.0000 - 0.0000i
      -0.0000 + 0.0000i  0.0001 + 0.0000i],...
     'AbsTol', 1e-4);
  
    verifyEqual(testCase, eValue, 1., 'RelTol', 1e-3);
end

function test_transferMatrix_errorMessage(testCase)
    % test a Vidal form mps throws the right error
    mps = {3,[1 0;0 1], [1 0;0 1]};
    
    verifyError(testCase, @()transferMatrix(mps),...
                'transferMatrix:VidalFormIncompatible');
end

%% setup function

function setup(testCase)
    % set up an input fixture for one loop of tdvpIntegrator

    % 1. p=3, D=2, FORM=1. Already in canonical form (to 3dp) and normalized
    testCase.TestData.state1 = ...
    {1,...
     cat(3, [-0.0096 + 0.3572i  0.3527 + 0.2393i; ...
             -0.0048 + 0.0048i  0.0199 - 0.0131i],...
            [-0.0552 + 0.6803i  0.3741 - 0.4573i; ...
              0.0031 - 0.0014i  0.0201 - 0.0176i],...
            [-0.2675 + 0.5787i -0.6735 + 0.1203i; ...
             -0.0001 - 0.0012i  -0.0002 - 0.0057i]),...
     diag([0.9999, 0.0001])};

    % ... and a two-site Hamiltonian to go with it
    loop = [1 0 -1 1];
    H1 = [];
    for n = 1:20
    H1 = cat(2, H1, loop);
    end
    H1 = cat(2, H1, 0);
    H1 = reshape(H1, [3 3 3 3]);
    testCase.TestData.H1 = H1;

    % ... and the C-tensor
    Clist = [-1.0679 - 0.4260i, -0.0021 + 0.0034i,  0.7262 - 0.2681i,...
         -0.0060 + 0.0072i, -0.4030 - 0.2119i, -0.0053 - 0.0068i,...
         -0.1285 - 0.9885i,  0.0004 + 0.0010i, -0.1756 - 0.0492i,...
          0.0004 - 0.0036i, -0.6301 + 0.5037i,  0.0055 - 0.0035i,...
         -0.8551 - 0.3867i,  0.0041 + 0.0035i,  0.1712 + 0.8657i,...
          0.0005 - 0.0073i, -1.0679 - 0.4260i, -0.0021 + 0.0034i,...
          0.7262 - 0.2681i, -0.0060 + 0.0072i, -0.4030 - 0.2119i,...
         -0.0053 - 0.0068i, -0.1285 - 0.9885i,  0.0004 + 0.0010i,...
         -0.1756 - 0.0492i,  0.0004 - 0.0036i, -0.6301 + 0.5037i,...
          0.0055 - 0.0035i, -0.8551 - 0.3867i,  0.0041 + 0.0035i,...
          0.1712 + 0.8657i,  0.0005 - 0.0073i, -0.8048 - 0.1172i,...
         -0.0028 + 0.0031i,  0.6148 + 0.1500i, -0.0062 + 0.0064i];

    testCase.TestData.C1 = reshape(Clist,[2 2 3 3]);

    % ... and the transfer matrix E
    E = [1.0000-0.0000i  0.0001-0.0000i  0.0001+0.0000i  0.9988+0.0000i;...
     -0.0000-0.0000i -0.0212+0.0180i  0.0012-0.0046i  0.0189+0.0029i;...
     -0.0000+0.0000i  0.0012+0.0046i -0.0212-0.0180i  0.0189-0.0029i;...
      0.0001+0.0000i -0.0001+0.0001i -0.0001-0.0001i  0.0013+0.0000i];

    testCase.TestData.E1 = E; 

    % ... and the highest-weight left/right eigenvectors as DxD matrices
    testCase.TestData.left1  = diag([0.7071,0.7071]);
    testCase.TestData.right1 = diag([1., 0.0001]);

    % ... and the K-tensor
    testCase.TestData.K1 = [0.0002 + 0.000i, -0.3284 - 0.088i; ...
                        -0.3285 + 0.088i  -1.8060 + 0.000i];

    % ... and the VL tensor (null vectors)
    Vlist =  [-0.6820 - 0.4747i,  0.0001 - 0.0001i,  -0.0001 + 0.0000i,...
            1.0000 + 0.0000i,  0.0000 + 0.0001i,  -0.0000 + 0.0000i,...
           -0.0000 + 0.0000i, -0.0000 - 0.0000i,   0.4264 - 0.0521i,...
            0.0002 - 0.0001i, -0.0002 + 0.0000i,  -0.0000 + 0.0000i,...
           -0.0002 - 0.0000i,  1.0000 + 0.0000i,  -0.0000 - 0.0000i,...
           -0.0000 - 0.0000i, -0.2017 + 0.2904i,   0.0000 - 0.0001i,...
            0.0002 + 0.0001i, -0.0000 + 0.0000i,   0.0002 + 0.0001i,...
           -0.0000 + 0.0000i,  0.0000 + 0.0001i,   1.0000 + 0.0000i];
    testCase.TestData.VL1 = reshape(Vlist, [2 4 3]);

    % ... and the F-tensor
    testCase.TestData.F1 = [-0.2817 - 0.6497i  -0.3104 - 0.2386i;...
                         -0.3643 + 0.0151i   0.0135 - 0.0066i;...
                         -0.3157 - 0.1393i   0.0047 + 0.0144i;...
                          0.1778 - 0.1214i  -0.0114 + 0.0035i];

    % 2. p=2, D=3, FORM=1. Example worke through by hand (not including F),
    %    starting with the intial state in canonical form
    elp = (sqrt(3)/6) * (2 + sqrt(2)); % Schmidt values
    elm = (sqrt(3)/6) * (2 - sqrt(2)); 
    
    A2(:,:,1) = [0  -elp;-elm 0  ];
    A2(:,:,2) = [elp 0  ; 0  -elm];

    left2  = eye(2);
    right2 = diag([elp^2 elm^2]);
    
    testCase.TestData.state2 = {1, A2, left2};
    
    testCase.TestData.left2  = left2;
    testCase.TestData.right2 = right2;

    testCase.TestData.H2 = 0.05 * cat(4, cat(3, [0 1;1 0], [1 0;0 1]),...
                                      cat(3, [1 0;0 1], [0 1;1 0]));                            

    C2(:,:,1,1) = 0.05 * [0 -(elp^2 - elp*elm);(elm^2 - elp*elm) 0];
    C2(:,:,2,2) = 0.05 * [0 -(elp^2 - elp*elm);(elm^2 - elp*elm) 0];
    C2(:,:,1,2) = 0.05 * [(elp^2 + elp*elm)  0;0 (elm^2 + elp*elm)];
    C2(:,:,2,1) = 0.05 * [(elp^2 + elp*elm)  0;0 (elm^2 + elp*elm)];                                 
     
    testCase.TestData.C2 = C2;
    
    testCase.TestData.E2 = [elp^2  0        0       elp^2;...
                            0     -elp*elm  elp*elm 0    ;...
                            0      elp*elm -elp*elm 0    ;...
                            elm^2 0         0       elm^2];
                        
    k0 = -0.05 * (elp^4 - elm^4);
    testCase.TestData.K2 = k0 * [0 1;1 0];
    
    % 3. p=2, D=4. Un-normalized, un-canonicalized MPS. Recorded 06/02/17 11:16
    testCase.TestData.A3 =  cat(3,...
            [0.1799+0.2185i 0.3715+0.0448i 0.7636+0.2992i 0.6493+0.9729i; ...
             0.2752+0.3588i 0.4661+0.4263i 0.0803+0.5472i 0.9591+0.686i;  ...
             0.3382+0.8888i 0.6336+0.0004i 0.2568+0.9193i 0.6316+0.3286i; ...
             0.1576+0.0635i 0.3647+0.8033i 0.6449+0.5648i 0.1577+0.7812i],...
            [0.7932+0.4179i 0.4789+0.435i 0.5904+0.8771i 0.2179+0.4615i;  ...
             0.35+0.9058i 0.6971+0.1925i 0.6878+0.1377i 0.173+0.956i;     ...
             0.5825+0.0388i 0.4548+0.8934i 0.3607+0.2559i 0.4049+0.0935i; ...
             0.4239+0.2424i 0.9348+0.926i 0.1629+0.9357i 0.6007+0.7328i]);

    % ... and a two-site Hamiltonian to go with it
    testCase.TestData.H3 = 0.15 * cat(4, cat(3, [0 1;1 0], [1 0;0 1]),...
                                      cat(3, [1 0;0 1], [0 1;1 0])); 

    % ... and the C-tensor
    %C1 =  

    testCase.TestData.C1 = reshape(Clist,[2 2 3 3]);

    % ... and the transfer matrix E
    E = [1.0000-0.0000i  0.0001-0.0000i  0.0001+0.0000i  0.9988+0.0000i;...
     -0.0000-0.0000i -0.0212+0.0180i  0.0012-0.0046i  0.0189+0.0029i;...
     -0.0000+0.0000i  0.0012+0.0046i -0.0212-0.0180i  0.0189-0.0029i;...
      0.0001+0.0000i -0.0001+0.0001i -0.0001-0.0001i  0.0013+0.0000i];

    testCase.TestData.E1 = E; 

    % ... and the highest-weight left/right eigenvectors as DxD matrices
    testCase.TestData.left1  = diag([0.7071,0.7071]);
    testCase.TestData.right1 = diag([1., 0.0001]);

    % ... and the K-tensor
    testCase.TestData.K1 = [0.0002 + 0.000i, -0.3284 - 0.088i; ...
                        -0.3285 + 0.088i  -1.8060 + 0.000i];

    % ... and the VL tensor (null vectors)
    Vlist =  [-0.6820 - 0.4747i,  0.0001 - 0.0001i,  -0.0001 + 0.0000i,...
            1.0000 + 0.0000i,  0.0000 + 0.0001i,  -0.0000 + 0.0000i,...
           -0.0000 + 0.0000i, -0.0000 - 0.0000i,   0.4264 - 0.0521i,...
            0.0002 - 0.0001i, -0.0002 + 0.0000i,  -0.0000 + 0.0000i,...
           -0.0002 - 0.0000i,  1.0000 + 0.0000i,  -0.0000 - 0.0000i,...
           -0.0000 - 0.0000i, -0.2017 + 0.2904i,   0.0000 - 0.0001i,...
            0.0002 + 0.0001i, -0.0000 + 0.0000i,   0.0002 + 0.0001i,...
           -0.0000 + 0.0000i,  0.0000 + 0.0001i,   1.0000 + 0.0000i];
    testCase.TestData.VL1 = reshape(Vlist, [2 4 3]);

    % ... and the F-tensor
    testCase.TestData.F1 = [-0.2817 - 0.6497i  -0.3104 - 0.2386i;...
                         -0.3643 + 0.0151i   0.0135 - 0.0066i;...
                         -0.3157 - 0.1393i   0.0047 + 0.0144i;...
                          0.1778 - 0.1214i  -0.0114 + 0.0035i];
end