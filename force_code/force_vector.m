function fg=force_vector(sweep_angle,root_of_chord,span_of_the_entire_wing,I,J,K,L,M,N,P,Q,R,S,gauss_points_x,gauss_points_y,gauss_points_z,weights_x,weights_y,force_local)
% sweep_angle=30;%degree
% span_of_the_entire_wing=4;
%sweep_angle_in_radians=sweep_angle*pi/180;
% root_of_chord=5;
%tip_chord=root_of_chord-((span_of_the_entire_wing/2)*tan(sweep_angle_in_radians));
% area=0.5*(span_of_the_entire_wing/2)*(root_of_chord+tip_chord);
% aspect=span_of_the_entire_wing/area;
%  coordinates_left_top_corner=[0 root_of_chord];
%  coordinates_left_bottom_corner=[0 0];
%  coordinates_right_bottom_corner=[span_of_the_entire_wing/2 0];
%  coordinates_right_top_corner=[span_of_the_entire_wing/2 tip_chord];
%  I=5;
%  J=5;
%  K=5;
%  L=5;
%  M=5;
%  N=5;
%  P=5;
%  Q=5;
%  R=5;
%  S=5;
 number_of_nodes_I_J=((I*J));
number_of_nodes_K_L=((K*L));
number_of_nodes_M_N=((M*N));
number_of_nodes_P_Q=((P*Q));
number_of_nodes_R_S=((R*S));
% weights_x=[1 1];
% weights_y=[1 1];
%% Transformed coordinates
%lagrange shape functions
% gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_z=[-1/sqrt(3) 1/sqrt(3)];
% coordinate_trans_x=zeros(1,length(coordinates_x));
% coordinate_trans_y=zeros(1,length(coordinates_y));
% for i=1:length(gauss_points_x)
%     for j=1:length(gauss_points_y)
%         B_1=0.25*(1-gauss_points_y(j))*(1-gauss_points_x(i));
%         B_2=0.25*(1+gauss_points_y(j))*(1-gauss_points_x(i));
%         B_3=0.25*(1+gauss_points_y(j))*(1+gauss_points_x(i));
%         B_4=0.25*(1-gauss_points_y(j))*(1+gauss_points_x(i));
%         B_matrix=[B_1 B_2 B_3 B_4];
%         for op=1:length(x)
%         coordinate_trans_x(1,op)=(B_1+B_1+B_1+B_1)*coordinates_x(1,op);
%         coordinate_trans_y(1,op)=(B_1+B_1+B_1+B_1)*coordinates_y(1,op);
%         end
%     end
% end
%hope=zeros((I*J+K*L+M*N+P*Q+R*S),1);
P_non_derive_I_J_x=zeros(1,I*J);
P_derive_I_J_x=zeros(1,I*J);
P_derive_K_L_x=zeros(1,K*L);
P_non_derive_K_L_x=zeros(1,K*L);
P_derive_M_N_x=zeros(1,M*N);
P_non_derive_M_N_x=zeros(1,M*N);
P_non_derive_P_Q_x=zeros(1,P*Q);
P_derive_P_Q_x=zeros(1,P*Q);
P_non_derive_R_S_x=zeros(1,R*S);
P_derive_R_S_x=zeros(1,R*S);
P_non_derive_I_J_y=zeros(1,I*J);
P_derive_I_J_y=zeros(1,I*J);
P_derive_K_L_y=zeros(1,K*L);
P_non_derive_K_L_y=zeros(1,K*L);
P_derive_M_N_y=zeros(1,M*N);
P_non_derive_M_N_y=zeros(1,M*N);
P_non_derive_P_Q_y=zeros(1,P*Q);
P_derive_P_Q_y=zeros(1,P*Q);
P_non_derive_R_S_y=zeros(1,R*S);
P_derive_R_S_y=zeros(1,R*S);
hope_1=zeros(5,(I*J+K*L+M*N+P*Q+R*S));
hope_local=zeros((I*J+K*L+M*N+P*Q+R*S),1);
%mn=1:1:I*J;
for op=1:length(gauss_points_y)
    for j=1:length(gauss_points_x)
        for ki=1:length(gauss_points_z)
P_non_derive_I_J_x(1,1)=1;
P_non_derive_I_J_x(1,2)=gauss_points_x(j);
P_derive_I_J_x(1,1)=0;
P_derive_I_J_x(1,2)=1;
P_non_derive_K_L_x(1,1)=1;
P_non_derive_K_L_x(1,2)=gauss_points_x(j);
P_derive_M_N_x(1,1)=0;
P_derive_M_N_x(1,2)=1;
P_non_derive_M_N_x(1,1)=1;
P_non_derive_M_N_x(1,2)=gauss_points_x(j);
P_derive_M_N_x(1,1)=0;
P_derive_M_N_x(1,2)=1;
P_non_derive_P_Q_x(1,1)=1;
P_non_derive_P_Q_x(1,2)=gauss_points_x(j);
P_derive_P_Q_x(1,1)=0;
P_derive_P_Q_x(1,2)=1;
P_non_derive_R_S_x(1,1)=1;
P_non_derive_R_S_x(1,2)=gauss_points_x(j);
P_derive_R_S_x(1,1)=0;
P_derive_R_S_x(1,2)=1;
P_non_derive_I_J_y(1,1)=1;
P_non_derive_I_J_y(1,2)=gauss_points_y(op);
P_derive_I_J_y(1,1)=0;
P_derive_I_J_y(1,2)=1;
P_non_derive_K_L_y(1,1)=1;
P_non_derive_K_L_y(1,2)=gauss_points_y(op);
P_derive_M_N_y(1,1)=0;
P_derive_M_N_y(1,2)=1;
P_non_derive_M_N_y(1,1)=1;
P_non_derive_M_N_y(1,2)=gauss_points_y(op);
P_derive_M_N_y(1,1)=0;
P_derive_M_N_y(1,2)=1;
P_non_derive_P_Q_y(1,1)=1;
P_non_derive_P_Q_y(1,2)=gauss_points_y(op);
P_derive_P_Q_y(1,1)=0;
P_derive_P_Q_y(1,2)=1;
P_non_derive_R_S_y(1,1)=1;
P_non_derive_R_S_y(1,2)=gauss_points_y(op);
P_derive_R_S_y(1,1)=0;
P_derive_R_S_y(1,2)=1;
[detJacob,inerse]=coordinate_transformation_Jacobian(sweep_angle,root_of_chord,span_of_the_entire_wing,gauss_points_x,gauss_points_y);
        for k=3:number_of_nodes_I_J
           P_derive_I_J_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_derive_I_J_x(1,k-1)-(k/(k+1)*P_derive_I_J_x(1,k-2))); 
           P_non_derive_I_J_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(op))*P_non_derive_I_J_x(1,k-1)-(k/(k+1)*P_non_derive_I_J_x(1,k-2))); 
        end
        for k=3:number_of_nodes_K_L
           P_derive_K_L_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_derive_K_L_x(1,k-1)-(k/(k+1)*P_derive_K_L_x(1,k-2))); 
           P_non_derive_K_L_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(op))*P_non_derive_K_L_x(1,k-1)-(k/(k+1)*P_non_derive_K_L_x(1,k-2))); 
        end
        for k=3:number_of_nodes_M_N
            P_derive_M_N_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_derive_M_N_x(1,k-1)-(k/(k+1)*P_derive_M_N_x(1,k-2))); 
            P_non_derive_M_N_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_non_derive_M_N_x(1,k-1)-(k/(k+1)*P_non_derive_M_N_x(1,k-2))); 
         end
        for k=3:number_of_nodes_P_Q
           P_derive_P_Q_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_derive_P_Q_x(1,k-1)-(k/(k+1)*P_derive_P_Q_x(1,k-2))); 
           P_non_derive_P_Q_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(op))*P_non_derive_P_Q_x(1,k-1)-(k/(k+1)*P_non_derive_P_Q_x(1,k-2))); 
        end
        for k=3:number_of_nodes_R_S
           P_derive_R_S_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(j))*P_derive_R_S_x(1,k-1)-(k/(k+1)*P_derive_R_S_x(1,k-2))); 
           P_non_derive_R_S_x(1,k)=((((2*k+1)/(k+1))*gauss_points_x(op))*P_non_derive_R_S_x(1,k-1)-(k/(k+1)*P_non_derive_R_S_x(1,k-2))); 
        end
        for k=3:number_of_nodes_I_J
           P_derive_I_J_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(j))*P_derive_I_J_y(1,k-1)-(k/(k+1)*P_derive_I_J_y(1,k-2))); 
           P_non_derive_I_J_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_I_J_y(1,k-1)-(k/(k+1)*P_non_derive_I_J_y(1,k-2))); 
        end
        for k=3:number_of_nodes_K_L
           P_derive_K_L_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(j))*P_derive_K_L_x(1,k-1)-(k/(k+1)*P_derive_K_L_x(1,k-2))); 
           P_non_derive_K_L_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_K_L_y(1,k-1)-(k/(k+1)*P_non_derive_K_L_y(1,k-2))); 
        end
         for k=3:number_of_nodes_M_N
            P_derive_M_N_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_derive_M_N_y(1,k-1)-(k/(k+1)*P_derive_M_N_y(1,k-2))); 
            P_non_derive_M_N_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_M_N_y(1,k-1)-(k/(k+1)*P_non_derive_M_N_y(1,k-2))); 
         end
        for k=3:number_of_nodes_P_Q
           P_derive_P_Q_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(j))*P_derive_P_Q_y(1,k-1)-(k/(k+1)*P_derive_P_Q_y(1,k-2))); 
           P_non_derive_P_Q_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_P_Q_y(1,k-1)-(k/(k+1)*P_non_derive_P_Q_y(1,k-2))); 
        end
        for k=3:number_of_nodes_R_S
           P_derive_R_S_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(j))*P_derive_R_S_y(1,k-1)-(k/(k+1)*P_derive_R_S_y(1,k-2))); 
           P_non_derive_R_S_y(1,k)=((((2*k+1)/(k+1))*gauss_points_y(op))*P_non_derive_R_S_y(1,k-1)-(k/(k+1)*P_non_derive_R_S_y(1,k-2))); 
        end
        final_IJ=P_non_derive_I_J_x.*P_non_derive_I_J_y;
        final_KL=P_non_derive_K_L_x.*P_non_derive_K_L_y;
        final_MN=P_non_derive_M_N_x.*P_non_derive_M_N_y;
        final_PQ=P_non_derive_P_Q_x.*P_non_derive_P_Q_y;
        final_RS=P_non_derive_R_S_x.*P_non_derive_R_S_y;
        opq=1;
        pew=1;
        omega=1;
        yt=1;
        Ltd=10;
        while opq<=(I*J+K*L+M*N+P*Q+R*S)
            final=[final_IJ(1,yt) final_KL(1,yt) final_MN(1,yt) final_PQ(1,yt) final_RS(1,yt)];
            hope_1(pew,opq)=final(omega);
            yt=yt+1;
            if rem(opq,I*J)==0
                pew=pew+1;
                yt=1;
                omega=omega+1;
            end
                opq=opq+1;
        end
        %hope(u,I*J)=[final_IJ(u) zeros(1,length(mn));zeros(1,m_n) final_IJ(u) zeros(1,length)]
        %force_local=[0;0;1;gauss_points_z(ki)*0;gauss_points_z(ki)*0];
        hope_local=hope_local+hope_1'*force_local*detJacob*weights_x(j)*weights_y(op);
    end
    
end
end
fg=hope_local;
end
