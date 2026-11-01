%Experiment 2
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
[w0,l0] = w_l_from_z0_h_epsr(70.7,1.58e-3,4.4,pi/2,2.4e9)
w1 = 0.003; l1 = 0.01; l2 = 0.02;
[ZB_Width,ZC_Width] =  deal(w0);
[ZB_Length,ZC_Length] = deal(l0);
[ZA_Width,ZD_Width] =  deal(w1);
[ZA_Length,ZD_Length] = deal(l1);
ZE_Length = l2 , ZE_Width = w1;
gndL=ZA_Length+ZB_Length+ZC_Length+ZD_Length;
gndW=2*ZE_Length;

%% 
ZA = traceRectangular(Length = ZA_Length,Width = ZA_Width, ...
    Center = [-((ZA_Length/2)+ZB_Length); 0]);
ZB = traceRectangular(Length = ZB_Length,Width = ZB_Width, ...
    Center = [-(ZB_Length/2); 0]);
ZC = traceRectangular(Length = ZC_Length,Width = ZC_Width, ...
    Center = [(ZC_Length/2); 0]);
ZD = traceRectangular(Length = ZD_Length,Width = ZD_Width, ...
    Center = [(ZC_Length+(ZD_Length/2)); 0]);
ZE = traceRectangular(Length = ZE_Width,Width = ZE_Length, ...
    Center = [0; -ZE_Length/2]);
filter = ZA+ZB+ZC+ZD+ZE;
show(filter)
%% 
substrate = dielectric(EpsilonR = 4.4,LossTangent = 0.0027,Thickness = 1.58e-3);
ground = traceRectangular(Length = gndL,Width = gndW,Center = [0 0]);
figure; 
show(ground)
%% 
pcb = pcbComponent;
pcb.BoardShape = ground;
pcb.BoardThickness = 1.58e-3;
pcb.Layers ={filter,substrate,ground};
pcb.FeedDiameter = w0/2;
pcb.FeedLocations = [0 -ZE_Length 1 3;-(ZA_Length+ZB_Length) 0 1 3;
    (ZC_Length+ZD_Length) 0 1 3];
figure;
show(pcb)

%% 

freqs=linspace(0.1e9,6e9,200);
s_par= sparameters(pcb, freqs);
spar=s_par.Parameters;
freq=s_par.Frequencies;

%% 

figure;
rfplot(s_par)

%% 

% Extract S-parameters
S11 = squeeze(spar(1,1,:));
S21 = squeeze(spar(2,1,:));
S31 = squeeze(spar(3,1,:));
S23 = squeeze(spar(2,3,:));

%% Convert to dB
S11_dB = 20*log10(abs(S11));
S21_dB = 20*log10(abs(S21));
S31_dB = 20*log10(abs(S31));
S23_dB = 20*log10(abs(S23));

%% -------------------------------
% 1. S-PARAMETER MAGNITUDE PLOT
% -------------------------------
figure;
plot(freq/1e9, S11_dB, 'b--'); hold on;
plot(freq/1e9, S21_dB, 'r+');
plot(freq/1e9, S31_dB, 'y');
grid on;
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('S-Parameter Magnitude of T-Junction');
legend('S11', 'S21', 'S31', 'Location','best');
% set(gca,'FontSize',14);

%% -------------------------------
% 2. PHASE PLOTS OF S21, S31
% -------------------------------
figure;
plot(freq/1e9, unwrap(angle(S21))*180/pi,'LineWidth',2); hold on;
plot(freq/1e9, unwrap(angle(S31))*180/pi, 'g--','LineWidth',2);
grid on;
xlabel('Frequency (GHz)');
ylabel('Phase (Degrees)');
title('Phase of S21 and S31');
legend('Phase(S21)', 'Phase(S31)');
% set(gca, 'FontSize', 14);
% Set the font size of the current axes (gca) to 14 points

%% 
% Phase Imbalance

Phase_Imbalance = unwrap(angle(S21))*180/pi - unwrap(angle(S31))*180/pi;
plot(freq/1e9, Phase_Imbalance,'LineWidth',2);
xlabel('Frequency (GHz)');
ylabel('Phase Imbalance (Degrees)');
title('Phase Imbalance Between S21 and S31');

%% -------------------------------
% 3. POWER SPLIT VERIFICATION
% -------------------------------
P21 = 10*log10(abs(S21).^2);
P31 = 10*log10(abs(S31).^2);

figure;
plot(freq/1e9, P21,'b--','LineWidth',2); hold on;
plot(freq/1e9, P31,'LineWidth',2);
grid on;
xlabel('Frequency (GHz)');
ylabel('Power (dB)');
title('Power Split: |S21|^2 and |S31|^2');
legend('Power to Port 2','Power to Port 3');

%% -------------------------------
% 4. RETURN LOSS ONLY (S11)
% -------------------------------
figure;
% plot(freq/1e9, S11_dB, 'Color',[0.8500 0.3250 0.0980]);
plot(freq/1e9, S11_dB,'LineWidth',2);

grid on;
xlabel('Frequency (GHz)');
ylabel('Return Loss (dB)');
title('Return Loss S11 of T-Junction');
% set(gca,'FontSize',14);

%% -------------------------------
% 5. ISOLATION PLOT: S23
% -------------------------------
figure;
plot(freq/1e9, S23_dB, 'g--','LineWidth',2);
grid on;
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
title('Isolation Between Output Ports: S23');




