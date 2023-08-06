function D_matrix_trans_1=D_matrix_lower_airfoil_trans_plane(sweep,cr,b,theta,thickness_of_each_ply,gauss_points_x,gauss_points_y,coeff_of_z_lower,h)
% coeff_of_z_lower=Lower_fitting_a_curve();
% coeff_of_z_upper=Upper_fitting_a_curve();
% sweep=30;
% h=1;
% cr=10*10^-3;
% b=20*10^-3;
thickness_lower=0.01;
Z=zeros(5,12);
% gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
z_2=0;
z_derive_1=0;
z_derive_2=0;
z_derive_3=0;
n=10;
z_lower=[0.000000  0.000000;
  0.002000 -0.009300;
  0.005000 -0.016000;
  0.010000 -0.022100;
  0.020000 -0.029500;
  0.030000 -0.034400;
  0.040000 -0.038100;
  0.050000 -0.041200;
  0.070000 -0.046200;
  0.100000 -0.051700;
  0.120000 -0.054700;
  0.150000 -0.058500;
  0.170000 -0.060600;
  0.200000 -0.063300;
  0.220000 -0.064700;
  0.250000 -0.066600;
  0.280000 -0.068000;
  0.300000 -0.068700;
  0.320000 -0.069200;
  0.350000 -0.069600;
  0.370000 -0.069600;
  0.400000 -0.069200;
  0.420000 -0.068800;
  0.450000 -0.067600;
  0.480000 -0.065700;
  0.500000 -0.064400;
  0.530000 -0.061400;
  0.550000 -0.058800;
  0.580000 -0.054300;
  0.600000 -0.050900;
  0.630000 -0.045100;
  0.650000 -0.041000;
  0.680000 -0.034600;
  0.700000 -0.030200;
  0.730000 -0.023500;
  0.750000 -0.019200;
  0.770000 -0.015000;
  0.800000 -0.009300;
  0.830000 -0.004800;
  0.850000 -0.002400;
  0.870000 -0.001300;
  0.890000 -0.000800;
  0.920000 -0.001600;
  0.940000 -0.003500;
  0.950000 -0.004900;
  0.960000 -0.006600;
  0.970000 -0.008500;
  0.980000 -0.010900;
  0.990000 -0.013700;
  1.000000 -0.016300];
x=z_lower(1:end,1);
z_lower_1=z_lower(1:end,1);
for po=1:length(coeff_of_z_lower)
    z_derive_1=z_derive_1+(coeff_of_z_lower(po)*(((factorial(length(coeff_of_z_lower))/(factorial(length(coeff_of_z_lower)-po)*factorial(po))))*po*(length(coeff_of_z_lower)-po)*((z_lower_1(po))^(po-1)))*((z_lower_1(po))^0.75)*((1-z_lower(po))^0.75));
    z_derive_2=z_derive_2+(coeff_of_z_lower(po))*(factorial(length(coeff_of_z_lower))/((factorial(length(coeff_of_z_lower)-po)*factorial(po))))*(0.75)*(z_lower_1(po))^(0.75-1)*(1-z_lower_1(po))^0.75;
     z_derive_3=z_derive_3+(coeff_of_z_lower(po)*(((factorial(length(coeff_of_z_lower))/(factorial(length(coeff_of_z_lower))*(0.75*(1-z_lower_1(po))^(0.75-1)*(z_lower_1(po))^(0.75)))))));
end
z_derive=z_derive_1+z_derive_2+z_derive_3;
y=atan(z_derive);
y=(y*180/pi);
thickness_lower_1=thickness_lower*sqrt(1+(tan(y))^2);
%step_size=thickness_lower_1/n; 
%for k=1:length(coeff_of_z_lower)
 theta=[0 90 90 0];
 thickness_of_each_ply=0.1*(ones(1,4));
z_1=(((x).^0.75).*(1-x).^(0.75));
for k=1:length(x)
z_2=z_2+(factorial(length(coeff_of_z_lower))/(factorial(length(coeff_of_z_lower)-k)*factorial(k)))*coeff_of_z_lower*(((x).^k).*(ones(1,length(x))-x).^(length(coeff_of_z_lower)-k))';
end
z_1=z_1';
z=z_1+z_2;
z_upper=z+ones(1,length(x))*(thickness_lower_1/2)+ones(1,length(x))*h;
z_l_lower=z+(ones(1,length(x))*(thickness_lower_1/2));
h=(z_upper-z_l_lower)./n;
integral_final=letsintegrate(z_upper,z_l_lower,sweep,cr,b,gauss_points_x,gauss_points_y);
for i=1:length(n)
    k=z_l_lower+i*h;
    for lo=1:length(z_l_lower)
    if rem(k(lo),3)==0
        integral_final=integral_final+2*(function_1(k(lo),sweep,cr,b,gauss_points_x,gauss_points_y));
    else
        integral_final=integral_final+3*(function_1(k(lo),sweep,cr,b,gauss_points_x,gauss_points_y));
    end
    end
    integral_final=integral_final*3*h(1)/8;
end
[stiffness]=abd_matrix_trans(theta,thickness_of_each_ply);
D_matrix_trans_1=integral_final'*stiffness*integral_final;
function func=function_1(k,sweep,cr,b,gauss_points_x,gauss_points_y)
% sweep=30;
% cr=10;
% b=20;
% gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
inverse_trans=inverse';
func=[0 0 0 0 0 (inverse_trans(1,2)+inverse_trans(2,2)) 0 0 0 0 0 1;0 0 0 0 (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 1 0];
end
function integral=letsintegrate(z_pp,z_l_lower,sweep,cr,b,gauss_points_x,gauss_points_y)
% sweep=30;
% cr=10*10^-3;
% b=20*10^-3;
% gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
% gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
inverse_trans=inverse';
integration_1=zeros(2,12);
integration_2=zeros(2,12);
for i=1:length(z_pp)
integration_1=integration_1+[0 0 0 0 0 (inverse_trans(1,2)+inverse_trans(2,2)) 0 0 0 0 0 1;0 0 0 0 (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 1 0];
integration_2=integration_2+[0 0 0 0 0 (inverse_trans(1,2)+inverse_trans(2,2)) 0 0 0 0 0 1;0 0 0 0 (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 1 0];
end
integral=integration_1+integration_2;
end
end