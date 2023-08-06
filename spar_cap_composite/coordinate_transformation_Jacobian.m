function [det_jacobian,Jacobian_inverse]=coordinate_transformation_Jacobian(sweep_angle,root_of_chord,span_of_the_entire_wing,gauss_points_x,gauss_points_y)
%sweep_angle=30;%degree
%span_of_the_entire_wing=10;
sweep_angle_in_radians=sweep_angle*pi/180;
%root_of_chord=5;
tip_chord=root_of_chord-((span_of_the_entire_wing/2)*tan(sweep_angle_in_radians));
% coordinates_left_top_corner=[0 root_of_chord];
% coordinates_left_bottom_corner=[0 0];
% coordinates_right_bottom_corner=[span_of_the_entire_wing/2 0];
% coordinates_right_top_corner=[span_of_the_entire_wing/2 tip_chord];
coordinates_x=[0 root_of_chord root_of_chord root_of_chord-tip_chord]';
coordinates_y=[0 0 span_of_the_entire_wing/2 span_of_the_entire_wing/2]';
% gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
% coordinates_transformation_along_x=0;
% coordinates_transformation_along_y=0;
coordinates=[coordinates_x coordinates_y];
for i=1:length(gauss_points_x)
    for j=1:length(gauss_points_y)
%N_1=0.25*(1-gauss_points_x(i))*(1-gauss_points_y(j));
%N_2=0.25*(1+gauss_points_x(i))*(1-gauss_points_y(j));
%N_3=0.25*(1+gauss_points_x(i))*(1+gauss_points_y(j));
%N_4=0.25*(1-gauss_points_x(i))*(1+gauss_points_y(j));
B_1_x=-0.25*(1-gauss_points_y(j));
B_1_y=-0.25*(1-gauss_points_x(i));
B_2_x=0.25*(1-gauss_points_y(j));
B_2_y=-0.25*(1+gauss_points_x(i));
B_3_x=0.25*(1+gauss_points_x(i));
B_3_y=0.25*(1+gauss_points_y(j));
B_4_x=-0.25*(1+gauss_points_y(j));
B_4_y=0.25*(1+gauss_points_x(i));
Jacobian=[B_1_x B_2_x B_3_x B_4_x;B_1_y B_2_y B_3_y B_4_y]*coordinates;
    end
end
det_jacobian=det(Jacobian);
Jacobian_inverse=inv(Jacobian);
end
