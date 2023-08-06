function curve_fit_Lower=Lower_fitting_a_curve()
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
x_guess=ones(1,length(z_lower));
airfoil_equation_lower_curve_2=@airfoil_equation_lower_curve;
curve_fit_Lower=lsqcurvefit(airfoil_equation_lower_curve_2,x_guess,z_lower(1:end,1),z_lower(1:end,2)');
z_predict=z_lower(1:end,1)';
thickness_ratio_lower_surface=0.01;
S_predict=zeros(length(z_lower),length(z_lower));
for i=1:length(z_predict)
   S_predict(i,1:length(z_predict))=(factorial(length(z_lower))/(factorial(length(z_lower)-i)*factorial(i)))*(((z_predict).^i).*((ones(1,length(z_predict))-z_predict).^(length(z_predict)-i))); 
end
% predicited_graph=((z_predict).^0.75).*((ones(1,length(z_predict))-z_predict).^0.75).*(curve_fit_Lower*S_predict)+ones(1,length(z_predict))*thickness_ratio_lower_surface;
% plot(z_lower(1:end,1)',predicited_graph,"b")
% hold on
% plot(z_lower(1:end,1),z_lower(1:end,2),"bo")
% hold off
function z=airfoil_equation_lower_curve(x,x_data)
x_data=x_data';
thickness_ratio_lower_surface=0.01;
S=zeros(length(x_data),length(x_data));
for i=1:length(x_data)
   S(i,1:length(x_data))=(factorial(length(x_data))/(factorial(length(x_data)-i)*factorial(i)))*(((x_data).^i).*((ones(1,length(x_data))-x_data).^(length(x_data)-i))); 
end
matrix=zeros(1,length(x_data));
for i=1:length(x_data)
    matrix(1,i)=x(i);
end
z=((x_data).^0.75).*((ones(1,length(x_data))-x_data).^0.75).*(matrix*S)+ones(1,length(x_data))*thickness_ratio_lower_surface;
end
end
