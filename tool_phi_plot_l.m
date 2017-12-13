function tool_phi_plot_l( r_int,r_ext,phi_int,phi_ext,phi_int_lt,phi_ext_lt,n_th,label_opt,cur_time_i,par_dt )
%TOOL_CPM_PLOT Summary of this function goes here
%   Detailed explanation goes here

[N_int,M]=size(phi_int);
[N_ext,~]=size(phi_ext);
tP_mx=zeros(M,n_th);
dth=2*pi/(n_th-1);
theta=0:dth:2*pi;

 x_int=zeros(n_th,N_int);
 y_int=zeros(n_th,N_int);
 x_ext=zeros(n_th,N_ext);
 y_ext=zeros(n_th,N_ext);
 for j=1:n_th
     [tP_mx(:,j),~]...
         =Calc_LP(M-1,cos(theta(j)));
     for i=1:N_int
         x_int(j,i)=theta(j);
         y_int(j,i)=r_int(i);
     end
     for i=1:N_ext
         x_ext(j,i)=theta(j);
         y_ext(j,i)=r_ext(i);
     end
 end
data_int=phi_int*(tP_mx);
data_ext=phi_ext*(tP_mx);


figure;
hold on

[C1,h1]=contour(x_int,y_int,data_int',[5,-5,10,-10],'b');
[C2,h2]=contour(x_ext,y_ext,data_ext',[5,-5,10,-10],'b');
plot(x_int(:,N_int),y_int(:,N_int),'b','LineWidth',1.0);
h1.LineWidth=1;
h2.LineWidth=1;

xlabel('\theta');
ylabel('r');
title(['t=' num2str(cur_time_i*par_dt)]);
if label_opt==1
    clabel(C1,h1,'manual','FontSize',16);
    clabel(C2,h2,'manual','FontSize',16);
end
jpg_name=strcat('l_fig_t_phi_',num2str(cur_time_i),'.fig');
saveas(gcf,jpg_name);
hold off
drawnow

end

