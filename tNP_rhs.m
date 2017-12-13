function [Lcm,Lcp] = tNP_rhs( r,cm_L,cp_L,gphi_r_R,gphi_th_R,d2phi_R,u_R,v_R,par_dt,par_epsi,par_gamma,D1 )
%NP_RHS_CP Summary of this function goes here
%   Detailed explanation goes here
global P_mx;

[N,M]=size(cm_L);
par_NP=par_epsi*par_gamma;

[cp_r_R,cp_th_R]=Calc_grad_R(r,cp_L,D1);
[cm_r_R,cm_th_R]=Calc_grad_R(r,cm_L,D1);
gcpgphi=cp_r_R.*gphi_r_R+cp_th_R.*gphi_th_R;
gcmgphi=cm_r_R.*gphi_r_R+cm_th_R.*gphi_th_R;
Lcp_R=-u_R.*cp_r_R-v_R.*cp_th_R+par_NP*...
    (gcpgphi+(cp_L*P_mx(1:M,:)).*d2phi_R);
Lcm_R=-u_R.*cm_r_R-v_R.*cm_th_R-par_NP*...
    (gcmgphi+(cm_L*P_mx(1:M,:)).*d2phi_R);
Lcp=R2L(Lcp_R,M);
Lcm=R2L(Lcm_R,M);
Lcp=Lcp+cp_L/par_dt;
Lcm=Lcm+cm_L/par_dt;

%outpart=R2L(-u_R.*cm_r_R-v_R.*cm_th_R,M);
%outpart1=R2L(cm_r_R.*gphi_r_R+cm_th_R.*gphi_th_R,M);
%outpart2=R2L(cp_r_R.*gphi_r_R+cp_th_R.*gphi_th_R,M);

    Lcp(1,:)=zeros(1,M);
    Lcm(1,:)=zeros(1,M);
    Lcp(1,1)=1;
    Lcm(1,1)=1;
Lcp(N,:)=-R2L((cp_L(N,:)*P_mx(1:M,:)).*gphi_r_R(N,:),M);
Lcm(N,:)=R2L((cm_L(N,:)*P_mx(1:M,:)).*gphi_r_R(N,:),M);
end

