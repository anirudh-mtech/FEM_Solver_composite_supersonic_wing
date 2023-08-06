function C=ques_3_with_arguments_transformed(theta)
%prompt="Enter the value of theta in degree";
%theta=input(prompt);
theta=(pi*theta)/180;
C=zeros(3,3);
C(1,1)=cos(theta)^2;
C(1,2)=sin(theta)^2;
C(1,3)=2*cos(theta)*sin(theta);
C(2,1)=sin(theta)^2;
C(2,2)=cos(theta)^2;
C(2,3)=-2*cos(theta)*sin(theta);
C(3,1)=-cos(theta)*sin(theta);
C(3,2)=cos(theta)*sin(theta);
C(3,3)=cos(theta)^2-sin(theta)^2;
end