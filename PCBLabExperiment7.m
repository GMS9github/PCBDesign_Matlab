%Experiment 7
clear all
clc
%% 

function [w,l] = w_l_from_z0_h_epsr(z0,h,epsr,theta,f)
    A=(z0/60)*sqrt((epsr+1)/2)+((epsr-1)/(epsr+1))*(0.23+(0.11/epsr));
    wh1=8*exp(A)/(exp(2*A)-2);
    B=(377*pi)/(2*z0*sqrt(epsr));
    wh2=(2/pi)*(B-1-log(2*B-1)+((epsr-1)/(2*epsr))*(log(B-1)+0.39-0.61/epsr));
    if wh1<2
        wh=wh1;
    else
        wh=wh2;
    end
    w = wh*h ;
    eps_eff = ((epsr+1)/2)+((epsr-1)/2)*(1/sqrt(1+12/wh));
    %epsr is relative permittivity of dielectric 
    c=3e8;
    k0=(2*pi*f)/c;
    beta = k0*sqrt(eps_eff);
    %beta is propagation constant and eps_eff is effective permittivity
    l=theta/beta;
end

%% 

function [trace] = tracepoint2rectangle(p1,p2,w)
%It will create rectangular trace from point p1 to p2 with width=w
% Author: Rakesh Sinha, NIT Rourkela
% Detailed explanation goes here
d=(p2-p1)/norm(p2-p1);
n=[-d(2),d(1)];
c1=p1+n*w/2;
c2=p1-n*w/2;
c3=p2-n*w/2;
c4=p2+n*w/2;
trace = antenna.Polygon(Vertices=[c1(1) c1(2) 0; c2(1) c2(2) 0; c3(1) c3(2) 0;c4(1),c4(2),0]);
end

%% 

function [trace] = diagonalpoint2rectangle(p1,p2)
%It will create rectangular trace from point p1 to p2 with width=w
% Author: Rakesh Sinha, NIT Rourkela
% Detailed explanation goes here
d=(p2-p1)/norm(p2-p1);
n=[-d(2),d(1)];
c1=p1;
c2=[p2(1),p1(2)];
c3=p2;
c4 = [p1(1), p2(2)];
trace = antenna.Polygon(Vertices=[c1(1) c1(2) 0; c2(1) c2(2) 0; c3(1) c3(2) 0;c4(1),c4(2),0]);
end

%% 

Z0=50;
Z1=Z0*sqrt(2);

%% 

[w0,l0] = w_l_from_z0_h_epsr(Z0,1.58e-3,4.4,pi/2,2.4e9);
[w1,l1] = w_l_from_z0_h_epsr(Z1,1.58e-3,4.4,pi/2,2.4e9);
[w2,l2] = w_l_from_z0_h_epsr(Z0,1.58e-3,4.4,pi/6,2.4e9);
%2*la+lb=l1+2;
[la,lb]=deal((l1+2e-3/2)/3);

%% 

ZA_Width =  deal(w0);
ZA_Length = deal(l0);
[ZB_Width,ZC_Width,ZD_Width,ZG_Width,ZH_Width,ZI_Width] =  deal(w1);
[ZB_Length,ZG_Length] = deal(la);
[ZC_Length,ZH_Length] = deal(lb);
[ZD_Length,ZI_Length] = deal(la-2e-3);
[ZE_Width,ZF_Width,ZJ_Width,ZK_Width] =  deal(w2);
[ZE_Length,ZF_Length,ZJ_Length,ZK_Length] = deal(l2);

%% 

P1=[0,0];
P2=[ZA_Length,0];
P3=[ZA_Length,ZB_Length];
P4=[ZA_Length+ZC_Length,ZB_Length];
P5=[ZA_Length+ZC_Length,2e-3];
P6=[ZA_Length+ZC_Length+(ZE_Length/sqrt(2)),(ZE_Length/sqrt(2))+2e-3];
P7=[ZA_Length+ZC_Length+(ZE_Length/sqrt(2)),(ZE_Length/sqrt(2))+2e-3+ZF_Length];
P8=[ZA_Length,-ZG_Length];
P9=[ZA_Length+ZH_Length,-ZG_Length];
P10=[ZA_Length+ZH_Length,-2e-3];
P11=[ZA_Length+ZH_Length+(ZJ_Length/sqrt(2)),-(ZJ_Length/sqrt(2))-2e-3];
P12=[ZA_Length+ZH_Length+(ZJ_Length/sqrt(2)),-(ZJ_Length/sqrt(2))-2e-3-ZK_Length];
P13=[0,-(ZJ_Length/sqrt(2))-2e-3-ZK_Length];
P14=[ZA_Length+ZH_Length+(ZJ_Length/sqrt(2))+ZF_Width,(ZE_Length/sqrt(2))+2e-3+ZF_Length];
ZA=tracepoint2rectangle(P1,P2,ZA_Width);
ZB=tracepoint2rectangle(P2,P3,ZB_Width);
ZC=tracepoint2rectangle(P3,P4,ZC_Width);
ZD=tracepoint2rectangle(P4,P5,ZD_Width);
ZE=tracepoint2rectangle(P5,P6,ZE_Width);
ZF=tracepoint2rectangle(P6,P7,ZF_Width);
ZG=tracepoint2rectangle(P2,P8,ZG_Width);
ZH=tracepoint2rectangle(P8,P9,ZH_Width);
ZI=tracepoint2rectangle(P9,P10,ZI_Width);
ZJ=tracepoint2rectangle(P10,P11,ZJ_Width);
ZK=tracepoint2rectangle(P11,P12,ZK_Width);
strip=ZA+ZB+ZC+ZD+ZE+ZF+ZG+ZH+ZI+ZJ+ZK;
ground = diagonalpoint2rectangle(P13,P14);
substrate = dielectric(EpsilonR = 4.4,LossTangent = 0.0027,Name = "custom",Thickness = 1.27e-3);
pcb = pcbComponent;
% Define nodes where the lumped element connects
%pcb.Nodes = [
%     node("P5",  [P5  0])     % [x y z]
%     node("P10", [P10 0])
% ];
% le = lumpedElement(Impedance = 50,TerminalNodes = ["P5","P10"] );
pcb.BoardShape = ground;
pcb.BoardThickness = 1.27e-3;
pcb.Layers ={strip,substrate,ground};
% add(pcb,le);
pcb.FeedDiameter = w0*0.7;
pcb.FeedLocations = [0 0 1 3; 
    ZA_Length+ZC_Length+(ZE_Length/sqrt(2)) (ZE_Length/sqrt(2))+2e-3+ZF_Length 1 3;
    ZA_Length+ZH_Length+(ZJ_Length/sqrt(2)) -(ZJ_Length/sqrt(2))-2e-3-ZK_Length 1 3];
figure();
show(pcb)

%% ----------------------------------------------
%  Simulate S-parameters
% -----------------------------------------------
freq = linspace(0.1e9, 6e9, 100);  % change as required
spar = sparameters(pcb, freq);      % your PCB model
freqGHz = freq/1e9;
figure;
rfplot(spar);

% Extract S-parameters
S11 = squeeze(spar.Parameters(1,1,:));
S21 = squeeze(spar.Parameters(2,1,:));
S31 = squeeze(spar.Parameters(3,1,:));
S22 = squeeze(spar.Parameters(2,2,:));
S33 = squeeze(spar.Parameters(3,3,:));
S23 = squeeze(spar.Parameters(2,3,:));   % isolation between 2 & 3
S32 = squeeze(spar.Parameters(3,2,:));

%% -----------------------------------------------
%  Return Loss (S11, S22, S33)
% -----------------------------------------------
figure;
plot(freqGHz, 20*log10(abs(S11)), 'LineWidth', 2); hold on;
plot(freqGHz, 20*log10(abs(S22)),'+', 'LineWidth', 2);
plot(freqGHz, 20*log10(abs(S33)),'--', 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Return Loss (dB)');
legend('S11', 'S22', 'S33');
title('Return Loss');

%% -----------------------------------------------
%  Isolation Plots (S21, S31, S23)
% -----------------------------------------------
figure;
plot(freqGHz, 20*log10(abs(S21)), 'LineWidth', 2); hold on;
plot(freqGHz, 20*log10(abs(S31)),'--', 'LineWidth', 2);
plot(freqGHz, 20*log10(abs(S23)),'+', 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Isolation (dB)');
legend('S21', 'S31', 'S23');
title('Isolation');

%% -----------------------------------------------
%  Power Split (|S21|², |S31|²)
% -----------------------------------------------
P21 = 10*log10(abs(S21).^2);
P31 = 10*log10(abs(S31).^2);

figure;
plot(freqGHz, P21, 'LineWidth', 2); hold on;
plot(freqGHz, P31,'--', 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Power (dB)');
legend('|S21|^2', '|S31|^2');
title('Power Split');

%% -----------------------------------------------
%  Phase Plots
% -----------------------------------------------
phaseS21 = unwrap(angle(S21))*180/pi;
phaseS31 = unwrap(angle(S31))*180/pi;

figure;
plot(freqGHz, phaseS21, 'LineWidth', 2); hold on;
plot(freqGHz, phaseS31,'--', 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Phase (deg)');
legend('∠S21', '∠S31');
title('Phase Response');

%% -----------------------------------------------
%  Phase Imbalance (Δφ = ∠S21 − ∠S31)
% -----------------------------------------------
phaseImbalance = phaseS21 - phaseS31;

figure;
plot(freqGHz, phaseImbalance, 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Phase Imbalance (deg)');
title('Phase Imbalance: ∠S21 - ∠S31');

%% -----------------------------------------------
%  Amplitude Imbalance (optional)
% -----------------------------------------------
ampImbalance = 20*log10(abs(S21)) - 20*log10(abs(S31));

figure;
plot(freqGHz, ampImbalance, 'LineWidth', 2);
grid on; xlabel('Frequency (GHz)'); ylabel('Amplitude Imbalance (dB)');
title('Amplitude Imbalance');

%% -----------------------------------------------
%  Smith Plot (S11 and S22 for example)
% -----------------------------------------------
figure;
smithplot(S11); hold on;
smithplot(S22);  
title('Smith Chart: S11 and S22');
