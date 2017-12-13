function [Ncm_mx,Ncp_mx] = tNP_lhs( r,D1,D2,M,par_dt,par_epsi,par_gamma,int_ext_opt )
%MX4CHEBOPS Summary of this function goes here
%   Detailed explanation goes here
N=length(r);
Ncm_mx=zeros(N,N,M);
Ncp_mx=zeros(N,N,M);
I=eye(N);
par_NP=par_epsi*par_gamma;
for j=1:M
    grad2=D2+diag(2./r)*D1-j*(j-1)*diag(1./r.^2);
    Ncm_mx(:,:,j)=1/par_dt*I-par_NP*grad2;
    Ncp_mx(:,:,j)=1/par_dt*I-par_NP*grad2;
    Ncm_mx(N,:,j)=D1(N,:);
    Ncp_mx(N,:,j)=D1(N,:);
    Ncm_mx(1,:,j)=I(1,:);
    Ncp_mx(1,:,j)=I(1,:);
    if int_ext_opt==1 && j==2
            Ncm_mx(1,:,2)=D1(1,:)+2/r(1)*I(1,:);
            Ncp_mx(1,:,2)=D1(1,:)+2/r(1)*I(1,:);
    end
end
end

