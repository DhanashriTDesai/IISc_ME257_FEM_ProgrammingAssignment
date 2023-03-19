clear; clc; 

load coordinates5.dat 
load connectivity5.dat 

a = coordinates5(:,:); 
b = connectivity5(:,:); 

coordinates = transpose(a); 
connectivity = transpose(b); 

nNodes = length(coordinates); 
nElements = length(connectivity); 

xh = coordinates(1,:); 
yh = coordinates(2,:);   

 

%% Exact Solution 'u=x' 
u = yh'; 

%% FE Solution 'uh' 
uh = [ 


];  

err = u-uh

err_metric  = norm(err)^2 