function curve_fit_upper=Upper_fitting_a_curve()
z_upper=[0.000000  0.000000;
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
x_guess=ones(1,length(z_upper));
airfoil_equation_lower_curve_2=@airfoil_equation_lower_curve;
curve_fit_upper=lsqcurvefit(airfoil_equation_lower_curve_2,x_guess,z_upper(1:end,1),z_upper(1:end,2)');
z_predict=z_upper(1:end,1)';
thickness_ratio_lower_surface=0.01;
S_predict=zeros(length(z_upper),length(z_upper));
for i=1:length(z_predict)
   S_predict(i,1:length(z_predict))=(factorial(length(z_upper))/(factorial(length(z_upper)-i)*factorial(i)))*(((z_predict).^i).*((ones(1,length(z_predict))-z_predict).^(length(z_predict)-i))); 
end
% predicited_graph=((z_predict).^0.75).*((ones(1,length(z_predict))-z_predict).^0.75).*(curve_fit_upper*S_predict)+ones(1,length(z_predict))*thickness_ratio_lower_surface;
% plot(z_upper(1:end,1)',predicited_graph,"b")
% hold on
% plot(z_upper(1:end,1),z_upper(1:end,2),"bo")
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
