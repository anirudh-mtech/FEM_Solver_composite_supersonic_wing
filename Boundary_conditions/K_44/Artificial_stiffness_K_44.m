function K_44=Artificial_stiffness_K_11(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_x,gauss_points_y,gauss_points_trans_x,gauss_points_trans_y,weights_x,weights_y,weights_trans_x,weights_trans_y)
%  I=4;
%  J=4;
%  K=4;
%   L=4;
%  M=4;
%   N=4;
%  P=4;
%  Q=4;
%   R=4;
%   S=4;
%   sweep_angle=30;%degree
%  span_of_the_entire_wing=4;
%  root_of_chord=5;
%  sweep_angle_in_radians=sweep_angle*pi/180;
%  tip_chord=root_of_chord-((span_of_the_entire_wing/2)*tan(sweep_angle_in_radians));
%   gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
%  gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
%  gauss_points_trans_x=0;
%  gauss_points_trans_y=0;
%  weights_trans_x=2;
%  weights_trans_y=2;
%   weights_x=[1 1];
%  weights_y=[1 1];
K_11_1=Artificial_stiffness_leading_edge_K_44(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_x,gauss_points_y,weights_x,weights_y);
K_11_2=Artificial_stiffness_root_K_44(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_x,gauss_points_y,weights_x,weights_y);
K_11_3=Artificial_stiffness_tip_K_44(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_trans_x,gauss_points_trans_y,weights_trans_x,weights_trans_y);
K_11_4=Artificial_stiffness_trailing_edge_K_44(I,J,K,L,M,N,P,Q,R,S,sweep_angle,span_of_the_entire_wing,root_of_chord,gauss_points_x,gauss_points_y,weights_x,weights_y);
K_44=K_11_1+K_11_2+K_11_3+K_11_4;
end