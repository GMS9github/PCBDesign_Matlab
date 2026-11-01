%Experiment 1
clear all
clc

Centerfreq=2e9;
freq=sort([linspace(1e9,3e9,20),2e9]);
mline=microstripLine(Width=0.002,Length=0.05);
mline=design(microstripLine,Centerfreq);
figure
show(mline)

%% 

spar = sparameters(mline,freq);
s=spar.Parameters;
f=spar.Frequencies;
abcd=s2abcd(s,50);
A=squeeze(abcd(1,1,:));
B=squeeze(abcd(1,2,:));
C=squeeze(abcd(2,1,:));
D=squeeze(abcd(2,2,:));
Zc=sqrt(B./C);
theta=angle(A+B./Zc);

%% 

figure
plot(f,real(Zc),'LineWidth',2);
hold on 
plot(f,imag(Zc), 'b--','LineWidth',2);
xlabel('Frequency (GHz)');
ylabel('Characteristic Impedance');
legend('Real Part (Resistance)', 'Imaginary Part (Reactance)');
title('Characteristic Impedance of Microstrip Line');
grid on;
%% 

figure
rfplot(spar);

%% 

figure
plot(f,theta*180/pi,'LineWidth',2)
xlabel('Frequency (GHz)');
ylabel('Electrical Length');
title('Electrical Length of Microstrip Line');

l=mline.Length
w=mline.Width
sub=mline.Substrate
