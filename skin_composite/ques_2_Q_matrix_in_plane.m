function Q=ques_2_Q_matrix_in_plane()
Stiffness=ques_1_compliance_in_plane();
Q=zeros(3,3);
Q(1,1)=Stiffness(1,1);%-(Stiffness(1,3)^2/Stiffness(3,3));
Q(1,2)=Stiffness(1,2);%-((Stiffness(1,3)*Stiffness(2,3))/Stiffness(3,3));
Q(2,2)=Stiffness(2,2);%-(Stiffness(2,3)^2/Stiffness(3,3));
Q(3,3)=Stiffness(3,3);
Q(2,1)=Q(1,2);
end
