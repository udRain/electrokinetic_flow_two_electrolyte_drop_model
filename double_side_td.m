%% Example of double side problem
% Summary of example objective
% An example with fluids on both sides and rigid spherical memberane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up all useful variables
clear all;
format long;
N_int=256;
N_ext=256;
r_min=0.3;
r_mid=1;
r_max=1.7;
ra_int=(r_min-r_mid)/2;
ra_ext=(r_max-r_mid)/2;
[x_int,DM_int]=chebdif(N_int,3);
[x_ext,DM_ext]=chebdif(N_ext,3);
D1_int=DM_int(:,:,1)/ra_int;
D2_int=DM_int(:,:,2)/ra_int^2;
D3_int=DM_int(:,:,3)/ra_int^3;
D1_ext=DM_ext(:,:,1)/ra_ext;
D2_ext=DM_ext(:,:,2)/ra_ext^2;
D3_ext=DM_ext(:,:,3)/ra_ext^3;
r_int=r_min-ra_int+ra_int*x_int;
r_ext=r_max-ra_ext+ra_ext*x_ext;

par_dt=1e-3;
N_dt=1.5e2;

par_epsi=1/500;
par_gamma=1;
par_beta=1/4;
par_alpha=1/2;
par_lambda=2;
par_nu=1;
par_psi=10;

I_int=eye(N_int);
I_ext=eye(N_ext);
IN_int=D1_int;
IN_int(1,:)=I_int(1,:);
IN_ext=D1_ext;
IN_ext(1,:)=I_ext(1,:);

global P_mx;        %matrix of Legendre polynomial P_j
global dP_mx;       %matrix of 1st derivative of Legendre polynomial dP_j/dtheta
global absc;        %Abscissae of Gauss Legendre quadrature
global wts;         %Weights of Gauss Legendre quadrature

setup_modes;
n_modes=6;

tt=0:par_dt:par_dt*(N_dt-1);
par_A=par_beta*2;
par_B=4;
par_a1=(par_A+par_B)/par_beta;
par_b1=(4*sqrt(par_alpha*par_beta)+4)/sqrt(par_alpha*par_beta);
c2_coef=3*par_B/2/(par_A+par_B);
c1_coef=3*par_A/(par_A+par_B);
c2_int_v1=par_psi*(1-c2_coef*exp(-par_b1/par_a1*tt));
c1_ext_v1=par_psi*(1-c1_coef*exp(-par_b1/par_a1*tt))/2;


c20_coef=-3/40*par_psi^2*par_nu*sqrt(par_alpha*par_beta)...
    /(1+sqrt(par_alpha*par_beta))*(1-exp(-par_b1/par_a1*tt)).^2/(1+par_lambda);
c20_lt=c20_coef(end);
psi_ext_lt=c20_lt*(1-1./r_ext.^2);
psi_int_lt=-c20_lt*(r_int.^5-r_int.^3);

delta_diff=sqrt(6*pi/5)*par_nu*par_epsi/(1+sqrt(par_alpha*par_beta))...
    *3/2*(1-exp(-par_b1/par_a1*tt))*(sqrt(par_alpha*par_beta)...
    *(1-par_A)+par_alpha*par_beta/par_lambda*par_A)/(par_A+par_B)...
    .*exp(-par_b1/par_a1*tt);

cm_int=zeros(N_int,n_modes);
cp_int=zeros(N_int,n_modes);
cm_ext=zeros(N_ext,n_modes);
cp_ext=zeros(N_ext,n_modes);
cm_int(:,1)=ones(N_int,1);
cp_int(:,1)=ones(N_int,1);
cm_ext(:,1)=ones(N_ext,1);
cp_ext(:,1)=ones(N_ext,1);
phi_int=zeros(N_int,n_modes);
%phi_int(:,2)=-par_psi*c2_coef*r_int;
phi_ext=zeros(N_ext,n_modes);
%phi_ext(:,2)=-par_psi*(r_ext+(1-c1_coef)/2./r_ext.^2);
psi_int=zeros(N_int,n_modes);
psi_ext=zeros(N_ext,n_modes);
dphi_int=zeros(N_int,n_modes);
%dphi_int(:,2)=-par_psi*ones(size(r_int));
dphi_ext=zeros(N_ext,n_modes);
%dphi_ext(:,2)=-par_psi*ones(size(r_ext));
gphi_r_R_int=dphi_int*P_mx(1:n_modes,:);
gphi_th_R_int=diag(1./r_int)*phi_int*dP_mx(1:n_modes,:);
gphi_r_R_ext=dphi_ext*P_mx(1:n_modes,:);
gphi_th_R_ext=diag(1./r_ext)*phi_ext*dP_mx(1:n_modes,:);
d2phi_R_int=zeros(N_int,16);
d2phi_R_ext=zeros(N_ext,16);
u_R_int=zeros(N_int,16);
u_R_ext=zeros(N_int,16);
v_R_int=zeros(N_int,16);
v_R_ext=zeros(N_int,16);
c1_ext_v=zeros(N_dt,1);
c2_int_v=zeros(N_dt,1);
q1_ext_v=zeros(N_dt,1);
q2_int_v=zeros(N_dt,1);
q4_ext_v=zeros(N_dt,1);
q3_int_v=zeros(N_dt,1);
c1_f=zeros(N_dt,n_modes);
c2_f=zeros(N_dt,n_modes);
q1_f=zeros(N_dt,n_modes);
q2_f=zeros(N_dt,n_modes);
q3_f=zeros(N_dt,n_modes);
q4_f=zeros(N_dt,n_modes);
u_diff1=zeros(N_dt,1);
u_diff2=zeros(N_dt,1);
fj1m=zeros(N_dt,1);
fj1p=zeros(N_dt,1);

%% System simulation
% We begin to iterate the whole system. Within each iteration, we first update
% the ion concentrations using the Nernst-Planck equations, then update the
% electric potential using the Poisson equation. Finally we update the
% fluid velocity field using the Stokes equations.

for i=1:N_dt
    %%%ion concentrations
    [NPm_mx_int,NPp_mx_int]=tNP_lhs(r_int,D1_int,D2_int,n_modes,par_dt,...
        par_epsi,par_gamma,0);
    [NPm_mx_ext,NPp_mx_ext]=tNP_lhs(r_ext,D1_ext,D2_ext,n_modes,par_dt,...
        par_epsi,1,1);
    [rhscm_int,rhscp_int]=tNP_rhs(r_int,cm_int,cp_int,gphi_r_R_int,...
        gphi_th_R_int,d2phi_R_int,u_R_int,v_R_int,par_dt,par_epsi,par_gamma,D1_int);
    [rhscm_ext,rhscp_ext]=tNP_rhs(r_ext,cm_ext,cp_ext,gphi_r_R_ext,...
        gphi_th_R_ext,d2phi_R_ext,u_R_ext,v_R_ext,par_dt,par_epsi,1,D1_ext);
    for j=1:n_modes
        cm_int(:,j)=NPm_mx_int(:,:,j)\rhscm_int(:,j);
        cp_int(:,j)=NPp_mx_int(:,:,j)\rhscp_int(:,j);
        cm_ext(:,j)=NPm_mx_ext(:,:,j)\rhscm_ext(:,j);
        cp_ext(:,j)=NPp_mx_ext(:,:,j)\rhscp_ext(:,j);
    end
    %plot(r_int,cm_int(:,1));
    %hold on
    %plot(r_int,cp_int(:,1));
    %%electric potential
    d2phi_L_int=(cm_int-cp_int)/par_epsi^2*par_beta/par_alpha/2;
    d2phi_L_ext=(cm_ext-cp_ext)/par_epsi^2/2;
    d2phi_R_int=d2phi_L_int*P_mx(1:n_modes,:);
    d2phi_R_ext=d2phi_L_ext*P_mx(1:n_modes,:);
    for j=0:n_modes-1
        in1_int=IN_int\(r_int.^(j+2).*[0;d2phi_L_int(2:N_int,j+1)]);
        in2p_int=r_int.^(-2-2*j).*in1_int;
        in2_int=IN_int\([0;in2p_int(2:N_int)]);
        in1_ext=IN_ext\(r_ext.^(j+2).*[0;d2phi_L_ext(2:N_ext,j+1)]);
        in2p_ext=r_ext.^(-2-2*j).*in1_ext;
        in2_ext=IN_ext\([0;in2p_ext(2:N_ext)]);
        c1_ext=(par_alpha*in1_int(N_int)-in1_ext(N_ext)+(par_alpha-1)*...
            (j*in2_ext(N_ext)-par_psi*(j==1)))/(j+1+par_alpha*j);
        c2_int=in2_ext(N_ext)-in2_int(N_int)-c1_ext;
        phi_int(:,j+1)=r_int.^j.*in2_int+c2_int*r_int.^j;
        phi_ext(:,j+1)=r_ext.^j.*in2_ext-c1_ext./r_ext.^(j+1);
        if j==1
            c1_ext_v(i)=c1_ext;
            c2_int_v(i)=c2_int;
        end
        c1_f(i,j+1)=c1_ext;
        c2_f(i,j+1)=c2_int;
        %dphi_int(:,j+1)=in1_int.*r_int.^(-j-2)+j*phi_int(:,j+1)./r_int;
        %dphi_ext(:,j+1)=(in1_ext+c1_ext).*r_int.^(-j-2)+j*phi_ext(:,j+1)./r_ext;
        %if j==1
        %    plot(r_int,dphi_int(:,2));
        %    hold on
            %plot(r_ext,cm_ext(:,1));
        %    plot(r_ext,dphi_ext(:,2));
            %plot(r_ext,cp_ext(:,1));
        %    hold off
        %    title(['t=' num2str(i)])
        %    drawnow
        %end
    end
    
    phi_int(:,2)=phi_int(:,2)-par_psi*r_int;
    phi_ext(:,2)=phi_ext(:,2)-par_psi*r_ext;
    dphi_int=D1_int*phi_int;
    dphi_ext=D1_ext*phi_ext;
    %plot(r_int,phi_int(:,2));
            %hold on
            %plot(r_ext,cm_ext(:,1));
            %plot(r_ext,phi_ext(:,2));
            %plot(r_ext,cp_ext(:,1));
            %hold off
            %title(['t=' num2str(i)])
            drawnow
    %phi_int(:,2)=phi_int(:,2)-par_psi*(r_int);
    %phi_ext(:,2)=phi_ext(:,2)-par_psi*(r_ext);
    %dphi_int(:,2)=dphi_int(:,2)-par_psi*ones(size(r_int));
    %dphi_ext(:,2)=dphi_ext(:,2)-par_psi*ones(size(r_ext));
    %plot(r_ext,phi_ext(:,2))
    %hold on
    %plot(r_int,phi_int(:,2))
    %hold off
    %drawnow
    
    c2_int_lt=c2_int_v1(i);
    c1_ext_lt=c1_ext_v1(i);
    phi_int_lt=(-par_psi+c2_int_lt)*r_int;
    phi_ext_lt=-(par_psi*r_ext+c1_ext_lt./r_ext.^2); 
   
    gphi_r_R_int=dphi_int*P_mx(1:n_modes,:);
    gphi_r_R_ext=dphi_ext*P_mx(1:n_modes,:);
    gphi_th_R_int=diag(1./r_int)*phi_int*dP_mx(1:n_modes,:);
    gphi_th_R_ext=diag(1./r_ext)*phi_ext*dP_mx(1:n_modes,:);
    
    %%fluid velocity field
    force_r_R_int=-gphi_r_R_int.*d2phi_R_int;
    force_th_R_int=-gphi_th_R_int.*d2phi_R_int;
    force_r_L_int=R2L(force_r_R_int,n_modes);
    force_th_L_int=R2dL(force_th_R_int,n_modes);
    r_force_th_L_int=diag(r_int)*force_th_L_int;
    dr_force_th_L_int=D1_int*r_force_th_L_int;
    rhs_int=force_r_L_int-dr_force_th_L_int;
    rhs_int=rhs_int*par_alpha*par_epsi*par_nu/par_lambda;
    
    force_r_R_ext=-gphi_r_R_ext.*d2phi_R_ext;
    force_th_R_ext=-gphi_th_R_ext.*d2phi_R_ext;
    force_r_L_ext=R2L(force_r_R_ext,n_modes);
    force_th_L_ext=R2dL(force_th_R_ext,n_modes);
    r_force_th_L_ext=diag(r_ext)*force_th_L_ext;
    dr_force_th_L_ext=D1_ext*r_force_th_L_ext;
    rhs_ext=force_r_L_ext-dr_force_th_L_ext;
    rhs_ext=rhs_ext*par_nu*par_epsi;
    
    for j=1:n_modes-1
        in1s_int=IN_int\([0;rhs_int(2:N_int,j+1)].*r_int.^(j+3));
        in2s_int=IN_int\([0;in1s_int(2:N_int)]./r_int.^(2*j+4));
        in3s_int=IN_int\([0;in2s_int(2:N_int)].*r_int);
        in4p_int=in3s_int.*r_int.^(2*j-2);
        in4s_int=IN_int\([0;in4p_int(2:N_int)]);
        ps_int=in4s_int.*r_int.^(2-j);
        in1s_ext=IN_ext\([0;rhs_ext(2:N_ext,j+1)].*r_ext.^(j+3));
        in2s_ext=IN_ext\([0;in1s_ext(2:N_ext)]./r_ext.^(2*j+4));
        in3s_ext=IN_ext\([0;in2s_ext(2:N_ext)].*r_ext);
        in4p_ext=in3s_ext.*r_ext.^(2*j-2);
        in4s_ext=IN_ext\([0;in4p_ext(2:N_ext)]);
        ps_ext=in4s_ext.*r_ext.^(2-j);
        
        q_mat=[0,1/2/(2*j+1),1/(2*j-1),0;1/2/(2*j+1)/(2*j+3),0,0,-1;...
            1/(2*j+3)/(2*j+1),-1/2,-1,0;1/(2*j+3),par_lambda,0,0];
        q_rhs=[-in4s_int(end);in4s_ext(end);in3s_int(end)-in3s_ext(end);...
            in2s_ext(end)-par_lambda*in2s_int(end)];
        q_vec=q_mat\q_rhs;
        q1_ext=q_vec(1);
        q2_int=q_vec(2);
        q3_int=q_vec(3);
        q4_ext=q_vec(4);
        
        %q2_int=1\(1+par_lambda)*(in2s_ext(N_ext)-par_lambda*in2s_int(N_int)-...
        %    (2*j+1)*(in3s_int(N_int)-in3s_ext(N_ext)-(2*j-1)*in4s_int(N_int)));
        %q3_int=-(2*j-1)*(in4s_int(N_int)+q2_int/2/(2*j+1));
        %q1_ext=(2*j+3)*(in2s_ext(N_ext)-par_lambda*in2s_int(N_int)-...
        %    par_lambda*q2_int);
        %q4_ext=q1_ext/2/(2*j+1)/(2*j+3)-in4s_ext(N_ext);
        
        if j==2
            %in2s_ext(end)-in2s_int(end)
            q1_ext_v(i)=q1_ext;
            q2_int_v(i)=q2_int;
            q3_int_v(i)=q3_int;
            q4_ext_v(i)=q4_ext;
            fj1m(i)=in3s_int(end)+q2_int/2+q3_int;
            fj1p(i)=in3s_ext(end)+q1_ext/35;
        end
        
            q1_f(j+1)=q1_ext;
            q2_f(j+1)=q2_int;
            q3_f(j+1)=q3_int;
            q4_f(j+1)=q4_ext;
        
        psi_int(:,j+1)=ps_int+q2_int/2/(2*j+1)*r_int.^(j+3)+q3_int/(2*j-1)*...
            r_int.^(j+1);
        psi_ext(:,j+1)=ps_ext-q1_ext/2/(2*j+1)/(2*j+3)*r_ext.^(-j)+q4_ext*...
            r_ext.^(2-j);
    end
    
    dpsi_int=D1_int*psi_int;
    dpsi_ext=D1_ext*psi_ext;
    u_diff1(i)=dpsi_int(end,3);
    u_diff2(i)=dpsi_ext(end,3);
    [u_R_int,v_R_int]=Calc_u_v(r_int,psi_int,dpsi_int);
    [u_R_ext,v_R_ext]=Calc_u_v(r_ext,psi_ext,dpsi_ext);
    
    trunc=1;
    if mod(i,N_dt)==0
        tool_phi_plot(r_int(trunc:end),r_ext(trunc:end),...
            phi_int(trunc:end,:),phi_ext(trunc:end,:),...
            phi_int_lt,phi_ext_lt,...
            64,1,i,par_dt);
        %tool_cpm_plot_l(r_int(trunc:end),r_ext(trunc:end),...
        %    (cm_int(trunc:end,:)+cp_int(trunc:end,:))/2,...
        %    (cm_ext(trunc:end,:)+cp_ext(trunc:end,:))/2,64,1,i,par_dt);
        %tool_q_plot_l(r_int(trunc:end),r_ext(trunc:end),...
        %    (cp_int(trunc:end,:)-cm_int(trunc:end,:))/2,...
        %    (cp_ext(trunc:end,:)-cm_ext(trunc:end,:))/2,64,1,i,par_dt);
        % tool_psi_plot_l(r_int(trunc:end),r_ext(trunc:end),...
        %     psi_int(trunc:end,:),psi_ext(trunc:end,:),psi_int_lt,-psi_ext_lt,...
        %     64,1,i,par_dt);
    end
end

