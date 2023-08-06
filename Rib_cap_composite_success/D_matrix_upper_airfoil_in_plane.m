function D_matrix_2=D_matrix_upper_airfoil_in_plane(sweep,cr,b,theta,thickness_of_each_ply,gauss_points_x,gauss_points_y,coeff_of_z_upper,h_1)
% coeff_of_z_lower=Lower_fitting_a_curve();
% coeff_of_z_upper=Upper_fitting_a_curve();
% sweep=30;
% cr=10*10^-3;
% b=20*10^-3;
% h_1=1;
gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
thickness_lower=0.01;
Z=zeros(5,12);
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
z_2=0;
z_derive_1=0;
z_derive_2=0;
z_derive_3=0;
n=10;
z_uppe_1=[0.000000  0.000000;
  0.002000  0.009500;
  0.005000  0.015800;
  0.010000  0.021900;
  0.020000  0.029300;
  0.030000  0.034300;
  0.040000  0.038100;
  0.050000  0.041100;
  0.070000  0.046200;
  0.100000  0.051800;
  0.120000  0.054800;
  0.150000  0.058500;
  0.170000  0.060600;
  0.200000  0.063200;
  0.220000  0.064600;
  0.250000  0.066400;
  0.270000  0.067300;
  0.300000  0.068500;
  0.330000  0.069200;
  0.350000  0.069600;
  0.380000  0.069800;
  0.400000  0.069700;
  0.430000  0.069500;
  0.450000  0.069200;
  0.480000  0.068400;
  0.500000  0.067800;
  0.530000  0.066600;
  0.550000  0.065600;
  0.570000  0.064500;
  0.600000  0.062500;
  0.620000  0.061000;
  0.650000  0.058500;
  0.680000  0.055500;
  0.700000  0.053300;
  0.720000  0.050900;
  0.750000  0.046900;
  0.770000  0.043900;
  0.800000  0.038900;
  0.820000  0.035300;
  0.850000  0.029400;
  0.870000  0.025100;
  0.900000  0.018100;
  0.920000  0.013100;
  0.950000  0.004900;
  0.970000 -0.000900;
  0.980000 -0.003900;
  0.990000 -0.007100;
  1.000000 -0.010400;
  ];
x=z_uppe_1(1:end,1);
z_lower_actual=z_uppe_1(1:end,1);
for po=1:length(coeff_of_z_upper)
    z_derive_1=z_derive_1+(coeff_of_z_upper(po)*(((factorial(length(coeff_of_z_upper))/(factorial(length(coeff_of_z_upper)-po)*factorial(po))))*po*(length(coeff_of_z_upper)-po)*((z_lower_actual(po))^(po-1)))*((z_lower_actual(po))^0.75)*((1-z_lower_actual(po))^0.75));
    z_derive_2=z_derive_2+(coeff_of_z_upper(po))*(factorial(length(coeff_of_z_upper))/((factorial(length(coeff_of_z_upper)-po)*factorial(po))))*(0.75)*(z_lower_actual(po))^(0.75-1)*(1-z_lower_actual(po))^0.75;
     z_derive_3=z_derive_3+(coeff_of_z_upper(po)*(((factorial(length(coeff_of_z_upper))/(factorial(length(coeff_of_z_upper))*(0.75*(1-z_lower_actual(po))^(0.75-1)*(z_lower_actual(po))^(0.75)))))));
end
z_derive=z_derive_1+z_derive_2+z_derive_3;
y=atan(z_derive);
y=(y*180/pi);
thickness_lower_1=thickness_lower*sqrt(1+(tan(y))^2);
%step_size=thickness_lower_1/n; 
%for k=1:length(coeff_of_z_lower)
% theta=[0 90 90 0];
% thickness_of_each_ply=0.1*(ones(1,4));
z_1=(((x).^0.75).*(1-x).^(0.75));
for k=1:length(x)
z_2=z_2+(factorial(length(coeff_of_z_upper))/(factorial(length(coeff_of_z_upper)-k)*factorial(k)))*coeff_of_z_upper*(((x).^k).*(ones(1,length(x))-x).^(length(coeff_of_z_upper)-k))';
end
z_1=z_1';
z=z_1+z_2;
z_upper=z-ones(1,length(x))*(thickness_lower_1/2);
z_l_lower=z-(ones(1,length(x))*(thickness_lower_1/2))-(ones(1,length(x))*h_1);
h=(z_upper-z_l_lower)./n;
integral_final=letsintegrate(z_upper,z_l_lower,sweep,cr,b,gauss_points_x,gauss_points_y);
for i=1:length(n)
    k=z_l_lower+i*h;
    for lo=1:length(z_l_lower)
    if rem(k(lo),3)==0
        integral_final=integral_final+2*(function_1(k(lo),sweep,cr,b));
    else
        integral_final=integral_final+3*(function_1(k(lo),sweep,cr,b));
    end
    end
    integral_final=integral_final*3*h(1)/8;
end
[stiffness]=abd_matrix(theta,thickness_of_each_ply);
D_matrix_2=integral_final'*stiffness*integral_final;
function func=function_1(k,sweep,cr,b)
% sweep=30;
% cr=10;
% b=20;
gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
inverse_trans=inverse';
func=[(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 k*(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0;0 0 0 (inverse_trans(1,2)+inverse_trans(2,1)) 0 0 0 0 0 (k)*(inverse_trans(1,2)+inverse_trans(2,2)) 0 0;0 (inverse_trans(1,2)+inverse_trans(2,1)) (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 k*(inverse_trans(2,2)+inverse_trans(2,1)) k*(inverse_trans(1,1)+inverse_trans(1,2)) 0 0 0];
end
function integral=letsintegrate(z_pp,z_l_lower,sweep,cr,b,gauss_points_x,gauss_points_y)
% sweep=30;
% cr=10*10^-3;
% b=20*10^-3;
gauss_points_x=[-1/sqrt(3) 1/sqrt(3)];
gauss_points_y=[-1/sqrt(3) 1/sqrt(3)];
[detJacobian,inverse]=coordinate_transformation_Jacobian(sweep,cr,b,gauss_points_x,gauss_points_y);
inverse_trans=inverse';
integration_1=zeros(3,12);
integration_2=zeros(3,12);
for i=1:length(z_pp)
integration_1=integration_1+[(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 z_pp(i)*(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0;0 0 0 (inverse_trans(1,2)+inverse_trans(2,1)) 0 0 0 0 0 (z_pp(i))*(inverse_trans(1,2)+inverse_trans(2,2)) 0 0;0 (inverse_trans(1,2)+inverse_trans(2,1)) (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 z_pp(i)*(inverse_trans(2,2)+inverse_trans(2,1)) z_pp(i)*(inverse_trans(1,1)+inverse_trans(1,2)) 0 0 0];
integration_2=integration_2+[(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0 z_l_lower(i)*(inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 0;0 0 0 (inverse_trans(1,2)+inverse_trans(2,1)) 0 0 0 0 0 (z_l_lower(i))*(inverse_trans(1,2)+inverse_trans(2,2)) 0 0;0 (inverse_trans(1,2)+inverse_trans(2,1)) (inverse_trans(1,1)+inverse_trans(2,1)) 0 0 0 0 z_l_lower(i)*(inverse_trans(2,2)+inverse_trans(2,1)) z_l_lower(i)*(inverse_trans(1,1)+inverse_trans(1,2)) 0 0 0];
end
integral=integration_1+integration_2;
end
end