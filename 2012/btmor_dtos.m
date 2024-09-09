function [Er,Ar,Br,Cr,Dar]=btmor_dtos(E1,J1,J2,J3,J4,B1,B2,Da,Zc,Zo,tol)
% 
Zon=E1*Zo;
[Uc,S,Ub] = svd(full(Zon'*Zc),0);
s0=diag(S);
ks=length(s0);
K0=ks;
while (sum(s0(K0-1:ks))<tol/2)&&(K0>2)
  K0=K0-1; 
end
r=K0;
%r=min(K0,fix_order);
  sigma_r=diag(S(1:r,1:r));
%S1=[speye(size(J1));-(J3*J2)\(J3*J1)];
Vc=Zc*Ub(:,1:r);
Vo=Zo*Uc(:,1:r);
TR=Vc*(diag(ones(r,1)./sqrt(sigma_r)));
TL=Vo*(diag(ones(r,1)./sqrt(sigma_r)));
% if(isempty(Da))
Er=TL'*E1*TR;J1_til=TL'*J1*TR;
J2_til=TL'*J2;
J3_til=J3*TR;
B1=[spalloc(size(B1,1),size(B1,2),0);B1];
%C1=B1';
B1_til=TL'*B1;
C1_til=B1'*TR;
Ar=(J1_til-J2_til*(J4\J3_til));
Br=(B1_til-J2_til*(J4\B2));
Cr=C1_til-B2'*(J4\J3_til);
Dar=-B2'*(J4\B2);
end
 