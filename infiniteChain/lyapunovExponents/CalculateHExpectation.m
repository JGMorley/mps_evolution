function [Exp_H,Exp_M] = CalculateHExpectation (A, E, H, right, left, VL, chi, d)
   
% This function outputs two (chi^2, chi^2) matrices 
% Exp_H is ~<d_A psi]H[d_A psi>
% Exp_M is ~<d_A d_A psi]H[psi>


% Find the inverse of the transfer operator. Defined such that the left and
% right environment are zero eigenstates of it. 

 Q=eye(chi^2)-reshape(right.',[chi^2,1])*reshape(left.',[1,chi^2]);
TransN=reshape(Q*inv(eye(chi^2)-Q*E*Q)*Q,chi,chi,chi,chi);

% TransN=pinv(eye(chi^2,chi^2)-E);
% TransN=reshape(TransN,chi,chi,chi,chi);



% Define "B(i)" which spans the tangent space of A. The chi^2 elements of B
% are the chi^2 basis elements of the tangent space.

basis=zeros(chi,chi);

B=zeros(chi,chi,d,chi^2);
for i=1:numel(basis)
    basis=zeros(chi,chi);
    basis(i)=1;
    B(:,:,:,i)=  ncon({sqrtm(inv(left)), VL, basis, sqrtm(inv(right))},{[1 -1],[1 2 -3],[2 3 -4],[3 -2]})/chi;
end

% Test they are orthonormal
% test=zeros(numel(basis),numel(basis));
% for i=1:numel(basis)
%      for j=1:numel(basis)
%         test(i,j)=round(ncon({left, B(:,:,:,i), conj(B(:,:,:,j)),right},{[1 2],[1 4 3],[2 5 3],[4 5]}),10);
%      end
%  end
% test;
 

%Do some contractions which will be used to find all the elements of H and
%M.

AA=ncon({A,A},{[-1 1 -3], [1 -2 -4]});
AB=ncon({A,B},{[-1 1 -3], [1 -2 -4 -5]});
BA=ncon({B,A},{[-1 1 -3 -5], [1 -2 -4]});
BB=ncon({B,B},{[-1 1 -3 -5], [1 -2 -4 -6]});

LeftTransBB=ncon({left,B,conj(B),TransN},{[1 2],[1 4 3 -3],[2 5 3 -4],[5 4 -2 -1]});

LeftTransAA=ncon({left,A,conj(A),TransN},{[1 2],[1 4 3],[2 5 3],[5 4 -2 -1]});


TransBB=ncon({TransN,B,conj(B),right},{[-1 -2 5 4],[4 1 3 -3],[5 2 3 -4],[1 2]});

TransAA=ncon({TransN,A,conj(A),right},{[-1 -2 5 4],[4 1 3],[5 2 3],[1 2]});

TransAB=ncon({TransN,A,conj(B),right},{[-1 -2 5 4],[4 1 3],[5 2 3 -3],[1 2]});
TransBA=ncon({TransN,B,conj(A),right},{[-1 -2 5 4],[4 1 3 -3],[5 2 3],[1 2]});

TransABBA=ncon({TransN,A,conj(B),TransBA},{[-1 -2 5 4],[4 1 3],[5 2 3 -4],[1 2 -3]});
TransABAB=ncon({TransN,A,conj(B),TransAB},{[-1 -2 5 4],[4 1 3],[5 2 3 -3],[1 2 -4]});
TransBABA=ncon({TransN,B,conj(A),TransBA},{[-1 -2 5 4],[4 1 3 -3],[5 2 3],[1 2 -4]});
TransBAAB=ncon({TransN,B,conj(A),TransAB},{[-1 -2 5 4],[4 1 3 -3],[5 2 3],[1 2 -4]});



% Exp_H is ~<d_A psi]H[d_A psi>
Exp_H=zeros(chi^2,chi^2);

% Exp_M is ~<d_A d_A psi]H[psi>

Exp_M=zeros(chi^2,chi^2);


for i=1:numel(basis)

    for j=1:numel(basis)
        obv1B=ncon({left,AB(:,:,:,:,i),H,conj(AB(:,:,:,:,j)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv1A=ncon({left,AA,H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv1=-obv1A+obv1B;
        
        obv2=ncon({left,AB(:,:,:,:,i),H,conj(BA(:,:,:,:,j)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv3=ncon({left,BA(:,:,:,:,i),H,conj(AB(:,:,:,:,j)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv4B=ncon({left,BA(:,:,:,:,i),H,conj(BA(:,:,:,:,j)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv4A=ncon({left,AA,H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv4=-obv4A+obv4B;
        
        obv5a=ncon({left,BB(:,:,:,:,i,j),H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv5b=ncon({left,BB(:,:,:,:,j,i),H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv5=(obv5a+obv5b);
        
        obv6a=ncon({left,AA,H,conj(BB(:,:,:,:,i,j)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv6b=ncon({left,AA,H,conj(BB(:,:,:,:,j,i)),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv6=(obv6a+obv6b);
        
        obv7B=ncon({left,AA,H,conj(AA),TransBB(:,:,i,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv7A=ncon({left,AA,H,conj(AA),TransAA},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv7=-obv7A+obv7B;
        
        obv8B=ncon({LeftTransBB(:,:,i,j),AA,H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv8A=ncon({LeftTransAA,AA,H,conj(AA),right},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv8=-obv8A+obv8B;

        obv9=ncon({left,BA(:,:,:,:,i),H,conj(AA),TransAB(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv10=ncon({left,AB(:,:,:,:,i),H,conj(AA),TransAB(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv11=ncon({left,AA,H,conj(BA(:,:,:,:,j)),TransBA(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv12=ncon({left,AA,H,conj(AB(:,:,:,:,j)),TransBA(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        
        obv13a=ncon({left,BA(:,:,:,:,i),H,conj(AA),TransBA(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv13b=ncon({left,BA(:,:,:,:,j),H,conj(AA),TransBA(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv13=(obv13a+obv13b);
        
        obv14a=ncon({left,AB(:,:,:,:,i),H,conj(AA),TransBA(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv14b=ncon({left,AB(:,:,:,:,j),H,conj(AA),TransBA(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv14=(obv14a+obv14b);
        
        obv15a=ncon({left,AA,H,conj(BA(:,:,:,:,j)),TransAB(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv15b=ncon({left,AA,H,conj(BA(:,:,:,:,i)),TransAB(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv15=(obv14a+obv14b);
        
        obv16a=ncon({left,AA,H,conj(AB(:,:,:,:,j)),TransAB(:,:,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv16b=ncon({left,AA,H,conj(AB(:,:,:,:,i)),TransAB(:,:,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv16=(obv16a+obv16b);
        
        obv17=ncon({left,AA,H,conj(AA),TransABBA(:,:,i,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        
        obv18=ncon({left,AA,H,conj(AA),TransBAAB(:,:,i,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});

        obv19a=ncon({left,AA,H,conj(AA),TransBABA(:,:,i,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv19b=ncon({left,AA,H,conj(AA),TransBABA(:,:,j,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv19=(obv19a+obv19b);
        
        obv20a=ncon({left,AA,H,conj(AA),TransABAB(:,:,i,j)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv20b=ncon({left,AA,H,conj(AA),TransABAB(:,:,j,i)},{[1 2],[1 3 7 8],[5 6 7 8],[2 4 5 6],[3 4]});
        obv20=(obv20a+obv20b);
        
        Exp_H(i,j)=obv1+obv2+obv3+obv4+obv7+obv8+obv9+obv10+obv11+obv12+obv17+obv18;%+obv6+obv14+obv15+obv19;

        Exp_M(i,j)= obv6+obv15+obv16+obv20;
    end
end


end

