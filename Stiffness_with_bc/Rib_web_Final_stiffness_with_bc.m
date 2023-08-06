clc
final_result=zeros(1,5);
gauss_points_in_x=[-1/sqrt(3) 1/sqrt(3)];
gauss_points_in_z=[-1/sqrt(3) 1/sqrt(3)];
 gauss_points_trans_x=0;
 gauss_points_in_y=[-1/sqrt(3) 1/sqrt(3)];
gauss_points_trans_y=0;
weights_in_x=[1 1];
 weights_in_y=[1 1];
weights_trans_x=2;
 weights_trans_y=2;
sweep=30;%degree
 b=13.56*10^3;%scaled to meter
 cr=6*10^3;%scaled to meter
 sweep_rad=30*(pi/180);
 b=b/2;
 tip=cr-tan(sweep_rad)*(b/2);
 Aspect_ratio=(b^2)/(0.5*b*0.5*(cr+tip));
 Area=0.5*b*0.5*(cr+tip);
 weight=38000*10;
 theta=[0 90 90 0];
 thickness_of_each_ply=0.001*ones(1,length(theta));
l_1=1000;
 c=1000;
 h=1000;
  h_1=1000;
  t_1=1000;
  h_2=1000;
  t_2=1000;
 force_local=[0;0;weight/Area;0;0];
 Final_result=0;
 Final_converge=zeros(1,5);
 iter=0;
%  force_vector_final=zeros(rep*rep+rep*rep+rep*rep+rep*rep+rep*S,1);
po=1;
tic
for boom=2:8
    for po=1:length(thickness_of_each_ply)
 I=boom;J=boom;K=boom;L=boom; M=boom;N=boom;P=boom;Q=boom;R=boom;S=boom;%skin
 I_spar=boom;J_spar=boom;K_spar=boom;L_spar=boom; M_spar=boom;N_spar=boom;P_spar=boom;Q_spar=boom;R_spar=boom;S_spar=boom;%spar
 I_rib=boom;J_rib=boom;K_rib=boom;L_rib=boom; M_rib=boom;N_rib=boom;P_rib=boom;Q_rib=boom;R_rib=boom;S_rib=boom;%ribs
 %Stiffness_of_skin=FEM_Main_file_skin(gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,sweep,b,cr,I,J,K,L,M,N,P,Q,R,S,thickness_of_each_ply(po),theta(po),weights_trans_x,weights_trans_y,t_1,h_1,h);
 % Stiffness_of_spars_web=FEM_Main_file_spars_web(gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,weights_trans_x,weights_trans_y,sweep,b,cr,I_spar,J_spar,K_spar,L_spar,M_spar,N_spar,P_spar,Q_spar,R_spar,S_spar,h,h_1,t_1,theta,thickness_of_each_ply);
%  Stiffness_of_spars_cap=FEM_Main_file_spars_cap(gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,weights_trans_x,weights_trans_y,sweep,b,cr,I_spar,J_spar,K_spar,L_spar,M_spar,N_spar,P_spar,Q_spar,R_spar,S_spar,l_1,h,h_1);
  Stiffness_of_rib_web=FEM_Main_file_rib_web(gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,weights_trans_x,weights_trans_y,sweep,b,cr,I,J,K,L,M,N,P,Q,R,S,h_2,t_2,h,theta(po),thickness_of_each_ply(po));
%  Stiffness_of_rib_cap=FEM_Main_file_ribs_cap(gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,weights_trans_x,weights_trans_y,sweep,b,cr,I_rib,J_rib,K_rib,L_rib,M_rib,N_rib,P_rib,Q_rib,R_rib,S_rib,l_1,h_2,h,theta,thickness_of_each_ply);
  Stiffness_of_bc=Artificial_Stiffness(I,J,K,L,M,N,P,Q,R,S,sweep,b,cr,gauss_points_in_x,gauss_points_in_y,gauss_points_trans_x,gauss_points_trans_y,weights_in_x,weights_in_y,weights_trans_x,weights_trans_y);
 Final_Stiffness=Stiffness_of_bc+Stiffness_of_rib_web;%+Stiffness_of_rib_cap+Stiffness_of_rib_web+Stiffness_of_spars_cap+Stiffness_of_spars_web;Stiffness_of_skin
 force_local=[0;0;(5.2079);0;0];
 force_vector_final=force_vector(sweep,cr,b,I,J,K,L,M,N,P,Q,R,S,gauss_points_in_x,gauss_points_in_y,gauss_points_in_z,weights_in_x,weights_in_y,force_local);
 bc_u=linspace(1,I*J,I*J);%,I*J+K*L+M*N+P*Q+R*S];
 bc_v=linspace(I*J+1,I*J+K*L,I*J);
 bc_w=I*J+K*L+2:2:I*J+K*L+M*N-1;
 bc_theta_x=linspace(I*J+K*L+M*N+1,I*J+K*L+M*N+P*Q,I*J);
 bc_theta_y=linspace(I*J+K*L+M*N+P*Q+1,I*J+K*L+M*N+P*Q+R*S,I*J);
ibc=0;
for i=1:length(bc_u)
    
force_vector_final=force_vector_final-(Final_Stiffness(:,bc_u(i))*ibc);
  Final_Stiffness(:,bc_u(i))=0;
  Final_Stiffness(bc_u(i),:)=0;
  Final_Stiffness(bc_u(i),bc_u(i))=1;
  force_vector_final(bc_u(i),1)=ibc;
end
for i=1:length(bc_v)
force_vector_final=force_vector_final-(Final_Stiffness(:,bc_v(i))*ibc);
  Final_Stiffness(:,bc_v(i))=0;
  Final_Stiffness(bc_v(i),:)=0;
  Final_Stiffness(bc_v(i),bc_v(i))=1;
  force_vector_final(bc_v(i),1)=ibc;
end
for i=1:length(bc_w)
force_vector_final=force_vector_final-(Final_Stiffness(:,bc_w(i))*ibc);
  Final_Stiffness(:,bc_w(i))=0;
  Final_Stiffness(bc_w(i),:)=0;
  Final_Stiffness(bc_w(i),bc_w(i))=1;
  force_vector_final(bc_w(i),1)=ibc;
end
for i=1:length(bc_theta_x)
force_vector_final=force_vector_final-(Final_Stiffness(:,bc_theta_x(i))*ibc);
  Final_Stiffness(:,bc_theta_x(i))=0;
  Final_Stiffness(bc_theta_x(i),:)=0;
  Final_Stiffness(bc_theta_x(i),bc_theta_x(i))=1;
  force_vector_final(bc_theta_x(i),1)=ibc;
end
for i=1:length(bc_theta_y)
force_vector_final=force_vector_final-(Final_Stiffness(:,bc_theta_y(i))*ibc);
  Final_Stiffness(:,bc_theta_y(i))=0;
  Final_Stiffness(bc_theta_y(i),:)=0;
  Final_Stiffness(bc_theta_y(i),bc_theta_y(i))=1;
  force_vector_final(bc_theta_y(i),1)=ibc;
end
 goal=det(Final_Stiffness);
 if goal==0
     disp("Fail")
     break
 else
     iter=iter+1;
     disp("Success")
     
%%bc
 final_disp=pinv(Final_Stiffness)*force_vector_final;
 values=find(final_disp>=0.001);
 convergence=final_disp(values);
 values_transverse=find(force_vector_final>0);
 end
 final_result(po)=max(final_disp(values_transverse));
 for i=1:length(gauss_points_in_x)
     Final_result=Final_result+final_result(po)*(sum(thickness_of_each_ply)/2)*weights_in_x(i);
 end
    end
     if goal==0
     disp("Fail")
     break
     end
Final_converge(boom)=Final_result;
end
toc
plot(2:boom,Final_converge(2:boom),"-o");
xlabel("Number of nodes for each of the displacement")
ylabel("Transverse displacement Rib web in mm")
title("Convergence test")