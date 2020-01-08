w1= (0:0.01:pi);
f = (0:100:5000);
fs = 10000;
w = (2*pi*f/fs);
b0 = 1/4;
b1 = (1*cos(-w)+i*sin(-w))/4;
b2 = (1*cos(-2*w)+i*sin(-2*w))/4;
b3 = (1*cos(-3*w)+i*sin(-3*w))/4;

ab = (b0+b1+b2+b3);
%plot(w, abs(ab));

%clear; 
num=[1 1 1 1]; 
den=[1 0 0 0]; 
[z p k] = tf2zp(num,den);
zplane(z,p)