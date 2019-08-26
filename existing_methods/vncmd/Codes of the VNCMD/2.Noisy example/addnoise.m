function noise = addnoise(N,Mean,STD)
% Generate Gaussian white noise with mean values Mean and standard deviation STD
% N: the number of samples
y=randn(1,N); 
y=y/std(y); 
y=y-mean(y); 
a=Mean; 
b=STD; 
y=a+b*y; 

noise = y;
end