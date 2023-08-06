function K_33_alpha_tip_trans=Artificial_stiffness_tip_K_11(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_x,gauss_points_y,weights_x,weights_y)
%  I=4;
% J=4;
% K=4;
%  L=4;
%  M=4;
%  N=4;
%  P=4;
%  Q=4;
%  R=4;
%  S=4;
%  sweep_angle=30;%degree
%  span_of_the_entire_wing=4;
 sweep_angle_in_radians=sweep_angle*pi/180;
%   root_of_chord=5;
 tip_chord=root_of_chord-((span_of_the_entire_wing/2)*tan(sweep_angle_in_radians));
%  gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
%  weights_x=[1 1];
% weights_y=[1 1];
%  u_o=zeros(I,1);
 P_non_derive_I_J_x=zeros(1,I*J);
P_non_derive_K_L_x=zeros(1,K*L);
P_non_derive_M_N_x=zeros(1,M*N);
P_non_derive_P_Q_x=zeros(1,P*Q);
P_non_derive_R_S_x=zeros(1,R*S);
P_non_derive_I_J_y=zeros(1,I*J);
P_non_derive_K_L_y=zeros(1,K*L);
P_non_derive_M_N_y=zeros(1,M*N);
P_non_derive_P_Q_y=zeros(1,P*Q);
P_non_derive_R_S_y=zeros(1,R*S);
P_non_derive_I_J_x_root=zeros(1,I*J);
P_non_derive_K_L_x_root=zeros(1,K*L);
P_non_derive_M_N_x_root=zeros(1,M*N);
P_non_derive_P_Q_x_root=zeros(1,P*Q);
P_non_derive_R_S_x_root=zeros(1,R*S);
P_non_derive_I_J_y_root=zeros(1,I*J);
P_non_derive_K_L_y_root=zeros(1,K*L);
P_non_derive_M_N_y_root=zeros(1,M*N);
P_non_derive_P_Q_y_root=zeros(1,P*Q);
 P_non_derive_R_S_y_root=zeros(1,R*S);
number_of_nodes_I_J_le_x=I*J;
number_of_nodes_K_L_le_x=K*L;
number_of_nodes_M_N_le_x=M*N;
number_of_nodes_P_Q_le_x=P*Q;
number_of_nodes_R_S_le_x=R*S;
number_of_nodes_I_J_le_y=I*J;
number_of_nodes_K_L_le_y=K*L;
number_of_nodes_M_N_le_y=M*N;
number_of_nodes_R_S_le_y=R*S;
 number_of_nodes_I_J=((I*J));
number_of_nodes_K_L=((K*L));
number_of_nodes_M_N=((M*N));
number_of_nodes_P_Q=((P*Q));
number_of_nodes_R_S=((R*S));
 number_of_nodes_I_J_le=((I*J));
number_of_nodes_K_L_le=((K*L));
number_of_nodes_M_N_le=((M*N));
number_of_nodes_P_Q_le=((P*Q));
number_of_nodes_R_S_le=((R*S));
coordinates_x=[0 root_of_chord root_of_chord root_of_chord-tip_chord]';
coordinates_y=[0 0 span_of_the_entire_wing/2 span_of_the_entire_wing/2]';
x_3=coordinates_x(3);
y_3=coordinates_y(3);
x_4=coordinates_x(4);
y_4=coordinates_y(4);
K_33_alpha_tip_trans=zeros(1,I*J);
for op=1:length(gauss_points_y)
    for j=1:length(gauss_points_x)
        P_non_derive_I_J_x(1,1)=1;
        P_non_derive_I_J_x(1,2)=(gauss_points_x(j));
        P_non_derive_K_L_x(1,1)=1;
P_non_derive_K_L_x(1,2)=(gauss_points_x(j));
P_non_derive_M_N_x(1,1)=1;
P_non_derive_M_N_x(1,2)=gauss_points_x(j);
P_non_derive_P_Q_x(1,1)=1;
P_non_derive_P_Q_x(1,2)=(gauss_points_x(j));
P_non_derive_R_S_x(1,1)=1;
P_non_derive_R_S_x(1,2)=(gauss_points_x(j));
P_non_derive_I_J_y(1,1)=1;
P_non_derive_I_J_y(1,2)=(gauss_points_y(op));
P_non_derive_K_L_y(1,1)=1;
P_non_derive_K_L_y(1,2)=(gauss_points_y(op));
P_non_derive_M_N_y(1,1)=1;
P_non_derive_M_N_y(1,2)=gauss_points_y(op);
P_non_derive_P_Q_y(1,1)=1;
P_non_derive_P_Q_y(1,2)=gauss_points_y(op);
P_non_derive_R_S_y(1,1)=1;
P_non_derive_R_S_y(1,2)=(gauss_points_y(op));
P_non_derive_I_J_x_root(1,1)=1;
P_non_derive_I_J_x_root(1,2)=(gauss_points_x(j));
P_non_derive_I_J_y_root(1,1)=1;
P_non_derive_I_J_y_root(1,2)=(gauss_points_y(op));
 P_non_derive_K_L_x_root(1,1)=1;
P_non_derive_K_L_x_root(1,2)=(gauss_points_x(j));
P_non_derive_K_L_y_root(1,1)=1;
P_non_derive_K_L_y_root(1,2)=(gauss_points_x(j));
P_non_derive_M_N_x_root(1,1)=1;
P_non_derive_M_N_x_root(1,2)=gauss_points_x(j);
P_non_derive_M_N_y_root(1,1)=1;
P_non_derive_M_N_y_root(1,2)=gauss_points_y(op);
P_non_derive_P_Q_x_root(1,1)=1;
P_non_derive_P_Q_x_root(1,2)=(gauss_points_x(j));
P_non_derive_P_Q_y_root(1,1)=1;
P_non_derive_P_Q_y_root(1,2)=gauss_points_y(op);
P_non_derive_R_S_x_root(1,1)=1;
P_non_derive_R_S_x_root(1,2)=(gauss_points_y(op));
P_non_derive_R_S_y_root(1,1)=1;
P_non_derive_R_S_y_root(1,2)=(gauss_points_y(op));
for k=3:number_of_nodes_I_J
    P_non_derive_I_J_x(1,k)=((((2*k+1)/(k+1))*(gauss_points_x(j))*P_non_derive_I_J_x(1,k-1)-(k/(k+1)*P_non_derive_I_J_x(1,k-2)))); 
end
for k=3:number_of_nodes_K_L
    P_non_derive_K_L_x(1,k)=((((2*k+1)/(k+1))*(gauss_points_x(j)))*P_non_derive_K_L_x(1,k-1)-(k/(k+1)*P_non_derive_K_L_x(1,k-2))); 
end
for k=3:number_of_nodes_M_N
     P_non_derive_M_N_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(op))*P_non_derive_M_N_x(1,k-1)-(k/(k+1)*P_non_derive_M_N_x(1,k-2))); 
end
 for k=3:number_of_nodes_P_Q
      P_non_derive_P_Q_x(1,k)=((((2*k+1)/(k+1))*(gauss_points_x(j)))*P_non_derive_P_Q_x(1,k-1)-(k/(k+1)*P_non_derive_P_Q_x(1,k-2))); 
 end
  for k=3:number_of_nodes_R_S
       P_non_derive_R_S_x(1,k)=((((2*k+1)/(k+1))*(gauss_points_x(j)))*P_non_derive_R_S_x(1,k-1)-(k/(k+1)*P_non_derive_R_S_x(1,k-2))); 
  end
   for k=3:number_of_nodes_I_J
       P_non_derive_I_J_y(1,k)=((((2*k+1)/(k+1))*(gauss_points_y(op)))*P_non_derive_I_J_y(1,k-1)-(k/(k+1)*P_non_derive_I_J_y(1,k-2))); 
   end
   for k=3:number_of_nodes_K_L
       P_non_derive_K_L_y(1,k)=((((2*k+1)/(k+1))*(gauss_points_y(op)))*P_non_derive_K_L_y(1,k-1)-(k/(k+1)*P_non_derive_K_L_y(1,k-2))); 
   end
    for k=3:number_of_nodes_M_N
        P_non_derive_M_N_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_M_N_y(1,k-1)-(k/(k+1)*P_non_derive_M_N_y(1,k-2))); 
    end
    for k=3:number_of_nodes_P_Q
        P_non_derive_P_Q_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_P_Q_y(1,k-1)-(k/(k+1)*P_non_derive_P_Q_y(1,k-2))); 
    end
     for k=3:number_of_nodes_R_S
         P_non_derive_R_S_y(1,k)=((((2*k+1)/(k+1))*(gauss_points_y(op)))*P_non_derive_R_S_y(1,k-1)-(k/(k+1)*P_non_derive_R_S_y(1,k-2))); 
     end
     for k=3:number_of_nodes_I_J_le_x
    P_non_derive_I_J_x_root(1,k)=((((2*k+1)/(k+1))*(1)*P_non_derive_I_J_x_root(1,k-1)-(k/(k+1)*P_non_derive_I_J_x_root(1,k-2)))); 
end
for k=3:number_of_nodes_K_L_le_x
    P_non_derive_K_L_x_root(1,k)=((((2*k+1)/(k+1))*(1))*P_non_derive_K_L_x_root(1,k-1)-(k/(k+1)*P_non_derive_K_L_x_root(1,k-2))); 
end
for k=3:number_of_nodes_M_N_le_x
     P_non_derive_M_N_x_root(1,k)=((((2*k+1)/(k+1))*1)*P_non_derive_M_N_x_root(1,k-1)-(k/(k+1)*P_non_derive_M_N_x_root(1,k-2))); 
end
 for k=3:number_of_nodes_P_Q_le_x
      P_non_derive_P_Q_x_root(1,k)=((((2*k+1)/(k+1))*1)*P_non_derive_P_Q_x_root(1,k-1)-(k/(k+1)*P_non_derive_P_Q_x_root(1,k-2))); 
 end
  for k=3:number_of_nodes_R_S_le_x
       P_non_derive_R_S_x_root(1,k)=((((2*k+1)/(k+1))*1)*P_non_derive_R_S_x_root(1,k-1)-(k/(k+1)*P_non_derive_R_S_x_root(1,k-2))); 
  end
   for k=3:number_of_nodes_I_J_le_y
       P_non_derive_I_J_y_root(1,k)=((((2*k+1)/(k+1))*(1))*P_non_derive_I_J_y_root(1,k-1)-(k/(k+1)*P_non_derive_I_J_y_root(1,k-2))); 
   end
   for k=3:number_of_nodes_K_L_le_y
       P_non_derive_K_L_y_root(1,k)=((((2*k+1)/(k+1))*(1))*P_non_derive_K_L_y_root(1,k-1)-(k/(k+1)*P_non_derive_K_L_y_root(1,k-2))); 
   end
    for k=3:number_of_nodes_M_N_le_y
        P_non_derive_M_N_y_root(1,k)=((((2*k+1)/(k+1))*1)*P_non_derive_M_N_y_root(1,k-1)-(k/(k+1)*P_non_derive_M_N_y_root(1,k-2))); 
    end
     for k=3:number_of_nodes_P_Q
        P_non_derive_P_Q_y_root(1,k)=((((2*k+1)/(k+1))*1)*P_non_derive_P_Q_y_root(1,k-1)-(k/(k+1)*P_non_derive_P_Q_y_root(1,k-2))); 
    end
     for k=3:number_of_nodes_R_S_le_y
         P_non_derive_R_S_y_root(1,k)=((((2*k+1)/(k+1))*(1))*P_non_derive_R_S_y_root(1,k-1)-(k/(k+1)*P_non_derive_R_S_y_root(1,k-2))); 
     end
     final_IJ=P_non_derive_I_J_x.*P_non_derive_I_J_y;
        final_KL=P_non_derive_K_L_x.*P_non_derive_K_L_y;
        final_MN=P_non_derive_M_N_x.*P_non_derive_M_N_y;
        final_PQ=P_non_derive_P_Q_x.*P_non_derive_P_Q_y;
        final_RS=P_non_derive_R_S_x.*P_non_derive_R_S_y;
        final_IJ_root=P_non_derive_I_J_x_root.*P_non_derive_I_J_y_root;
        final_KL_root=P_non_derive_K_L_x_root.*P_non_derive_K_L_y_root;
        final_MN_root=P_non_derive_M_N_x_root.*P_non_derive_M_N_y_root;
        final_PQ_root=P_non_derive_P_Q_x_root.*P_non_derive_P_Q_y_root;
        final_RS_root= P_non_derive_R_S_x_root.*P_non_derive_R_S_y_root;
        K_33_alpha_tip_trans=K_33_alpha_tip_trans+((x_3-x_4)/4)*(final_IJ_root.*final_KL_root.*final_PQ_root.*final_MN_root.*final_RS_root)'*(final_PQ)*weights_x(j)*weights_y(op);
    end
end
end