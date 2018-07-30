function [ OUTPUT,energy ] = tdvpIntegrator_variabletimestep(A0, H, T, dt, bins, M, Energy_threshold)
% tdvpIntegrator_variabletimestep does integration using a variable
% timestep if a certain energy threshold is exceeded in a certain time. 
%
% A0 is the initial state.
%
% H is the Hamiltonian.
%
% T is the total time. 
%
% dt is the "default timestep". 
% The evolution is divided in to "bins". After each bin the energy will be
% checked. 
%
% If the energy change exceeded a certain amount it will be repeated with a
% timestep that is dt/M. 
%
% Energy_threshold determines the allowed energy threshold.

N=round(T/(dt*bins));
DT=dt;
energy=zeros([N+1,1]);

OUTPUT={A0,0};

[A,R,L,~] = normalizeMPS(A0);


% Find the initial energy. 
E=ncon({L,A,A,H,conj(A),conj(A),R},{[1 2],[1 5 3],[5 9 6],[4 8 3 6],[2 7 4],[7 10 8],[9 10]});
energy(1)=E;



progressWindow = waitbar(0, 'Variable timestep progress');
for n=1:N
    Ainit=OUTPUT{end,1};
%     [Ainit, R, L] = normalizeMPS(Ainit);
    Ainit;
    % Do the standard time evolution. 
    output_n = tdvpIntegrator(Ainit,H,bins*dt,dt,'INVERSE_FREE',true,'SYM_STEP',true,'DMAX',3);
    A_n=output_n{end,1};
    [A,R,L,~] = normalizeMPS(A_n);
    % Find the energy after this evolution. 
    E=ncon({L,A,A,H,conj(A),conj(A),R},{[1 2],[1 5 3],[5 9 6],[4 8 3 6],[2 7 4],[7 10 8],[9 10]});
    while abs(E-energy(n)) > Energy_threshold 
        % If the energy change is too large then we will repeat this
        % evolution.
        X=['Change in energy ', num2str(abs(E-energy(n))), 'at time ', num2str(n*bins*DT), char(10),'changing timestep to ', num2str(dt/M), char(10)];
        disp(X)
        dt=dt/M;
        % The function calls itself recursively, this may lead to four or
        % five layers of recursion. 
        [output_n,~]=tdvpIntegrator_variabletimestep(Ainit, H, bins*DT, dt, bins,M,Energy_threshold);

        A_n=output_n{end,1};
        [A,R,L,~] = normalizeMPS(A_n);
        % Check the new energy change. 
        E=ncon({L,A,A,H,conj(A),conj(A),R},{[1 2],[1 5 3],[5 9 6],[4 8 3 6],[2 7 4],[7 10 8],[9 10]});
    end
    dt=DT;
    energy(n+1)=E;
    % Form a new output with the correct times. The times from output_n
    % need to be offset by some amount. 
    output_temp1=output_n(2:end,1);
    output_temp2=cellfun(@(x) x+dt*bins*(n-1),output_n(2:end,2));
    output_temp2=mat2cell(output_temp2,ones(size(output_temp2,1),1));
    
    output_temp=horzcat(output_temp1,output_temp2);
    
    % The output from this section is then combine with the rest of the
    % time evolution. 
    OUTPUT=vertcat(OUTPUT,output_temp);
    waitbar_msg = ['Progress ',num2str(100*n/N,'%.2f'),'%'];
    waitbar(n/N,progressWindow,waitbar_msg)

end
    close(progressWindow)

end

