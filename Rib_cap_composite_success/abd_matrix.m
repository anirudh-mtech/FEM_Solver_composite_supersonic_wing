function stiffness=abd_matrix(theta_of_each_ply,~)
%theta_of_each_ply=[0 45 45 0];
%thickness_of_each_ply=ones(1,length(theta_of_each_ply));
% height_of_lamina=sum(thickness_of_each_ply);
% ords1=zeros(1,(length(thickness_of_each_ply)+1));
% ords1(1,1)=-height_of_lamina/2;
% for i=2:(length(thickness_of_each_ply))
%     ords1(1,i)=ords1(1,i-1)+thickness_of_each_ply(1,i);
% end
% ords1(1,end)=abs(ords1(1,1));
Q_transformed_1=zeros(3,3);
for i=1:length(theta_of_each_ply)
    Q_transformed=inv(ques_3_with_arguments_transformed_in_plane(theta_of_each_ply(i)))*ques_2_Q_matrix_in_plane*(inv(ques_3_with_arguments_transformed_in_plane(theta_of_each_ply(i)))');
%     abd(1:3,1:3)=abd(1:3,1:3)+Q_transformed*((ords1(1,i+1))-ords1(1,i));
%     abd(4:end,1:3)=abd(4:end,1:3)+0.5*Q_transformed*((ords1(1,i+1))^2-ords1(1,i)^2);
%     abd(1:3,4:6)=abd(1:3,4:6)+0.5*Q_transformed*((ords1(1,i+1))^2-ords1(1,i)^2);
%     abd(4:6,4:6)=abd(4:6,4:6)+(1/3)*Q_transformed*((ords1(1,i+1))^3-ords1(1,i)^3);
Q_transformed_1=Q_transformed_1+Q_transformed;
end
% z_k=0;
% z_k_2=0;
% for i=1:length(theta_of_each_ply)
% z_k=z_k+(ords1(1,i+1)-ords1(1,i));
% z_k_2=z_k_2+0.5*((ords1(1,i+1))^2-(ords1(1,i))^2);
% end
% complaince=[1 0 0 1 0 0;0 1 0 0 1 0;0 0 1 0 0 1]*inv(abd)*[z_k 0 0;0 z_k 0;0 0 z_k;(((z_k_2)^2)/2) 0 0;0 (((z_k_2)^2)/2) 0;0 0 (((z_k_2)^2)/2)];
stiffness=inv(Q_transformed_1);
end