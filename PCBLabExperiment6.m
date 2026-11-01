%Experiment 6
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

Z0=50;
Z1=Z0*sqrt(2);

%% 

[w0,l0] = w_l_from_z0_h_epsr(Z0,1.58e-3,4.4,pi/2,2.4e9);
[w1,l1] = w_l_from_z0_h_epsr(Z1,1.58e-3,4.4,pi/2,2.4e9);
[ZA_Width,ZC_Width,ZF_Width,ZG_Width] =  deal(w0);
[ZA_Length,ZC_Length,ZF_Length,ZG_Length] = deal(l0);
[ZB_Width,ZD_Width,ZE_Width,ZH_Width,ZI_Width,ZJ_Width] =  deal(w1);
[ZB_Length,ZD_Length,ZE_Length,ZH_Length,ZI_Length,ZJ_Length] = deal(l1);
gndL = ZA_Length+ZB_Length+ZC_Length;
gndW = 1.5*(ZD_Length+ZH_Length+ZB_Width);

%% 
ZA = traceRectangular(Length = ZA_Length,Width = ZA_Width, ...
    Center = [-(ZA_Length/2); 0]);
ZB = traceRectangular(Length = ZB_Length,Width = ZB_Width, ...
    Center = [(ZB_Length/2); 0]);
ZC = traceRectangular(Length = ZC_Length,Width = ZC_Width, ...
    Center = [ZB_Length+(ZC_Length/2);0]);
ZD = traceRectangular(Length = ZD_Width,Width = ZD_Length, ...
    Center = [0;-(ZD_Length/2)]);
ZE = traceRectangular(Length = ZE_Width,Width = ZD_Length, ...
    Center = [ZB_Length;-(ZE_Length/2)]);
ZF = traceRectangular(Length = ZF_Length,Width = ZF_Width, ...
    Center = [-(ZF_Length/2);-ZD_Length]);
ZG = traceRectangular(Length = ZG_Length,Width = ZG_Width, ...
    Center = [ZB_Length+(ZG_Length/2);-ZD_Length]);
ZH = traceRectangular(Length = ZH_Width,Width = ZH_Length, ...
    Center = [0; -(ZD_Length+ZH_Length/2)]);
ZI = traceRectangular(Length = ZI_Length,Width = ZI_Width, ...
    Center = [ZB_Length/2; -(ZD_Length+ZH_Length)]);
ZJ = traceRectangular(Length = ZJ_Width,Width = ZJ_Length, ...
    Center = [ZB_Length; -(ZD_Length+ZH_Length/2)]);
filter = ZA+ZB+ZC+ZD+ZE+ZF+ZG+ZH+ZI+ZJ;
show(filter)

%% 
substrate = dielectric(EpsilonR = 4.4,LossTangent = 0.0027,Name = "custom", ...
    Thickness = 1.58e-3);
ground = traceRectangular(Length = gndL,Width = gndW,Center = [(ZB_Length/2) -ZD_Length]);
figure; 
show(ground)

%% 
pcb = pcbComponent;
pcb.BoardShape = ground;
pcb.BoardThickness = 1.58e-3;
pcb.Layers ={filter,substrate,ground};
pcb.FeedDiameter = w0/2;
pcb.FeedLocations = [-ZA_Length 0 1 3;(ZB_Length+ZC_Length) 0 1 3;
    (ZB_Length+ZC_Length) -ZD_Length 1 3;-ZA_Length -ZD_Length 1 3];
figure;
show(pcb)

%% 

freqs=linspace(0.1e9,6e9,100);
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
S32 = squeeze(spar(3,2,:));
S41 = squeeze(spar(4,1,:));

%%

% Convert to dB
S11_dB = 20*log10(abs(S11));
S21_dB = 20*log10(abs(S21));
S31_dB = 20*log10(abs(S31));
S32_dB = 20*log10(abs(S32));
S41_dB = 20*log10(abs(S41));

%% -------------------------------
% 1. S-PARAMETER MAGNITUDE PLOT
% -------------------------------

figure;
plot(freq/1e9, S11_dB, 'b--'); hold on;
plot(freq/1e9, S21_dB, 'r+');
plot(freq/1e9, S31_dB, 'y*');
plot(freq/1e9, S41_dB, 'g');
grid on;
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
legend('S11','S21','S31','S41');
title('S-Parameter Magnitude of RRC');
legend('S11', 'S21', 'S31', 'S41' , 'Location','best');
% set(gca,'FontSize',14);

%% -------------------------------
% 2. PHASE Plots OF S21 and S31
% -------------------------------

figure;
plot(freq/1e9, (angle(S21))*180/pi,'LineWidth',2); hold on;
plot(freq/1e9, (angle(S31))*180/pi,'--','LineWidth',2);
grid on; 
xlabel('GHz'); 
ylabel('Phase (deg)');
title(['RRC Phase of S21 and S31']);
legend('Phase S21','Phase S31');

%% -------------------------------
% 3. PHASE IMBALANCE BETWEEN S21 and S31
% ------------------------------- 

Phase_Imbalance = unwrap(angle(S21))*180/pi - unwrap(angle(S31))*180/pi;
plot(freq/1e9, Phase_Imbalance,'LineWidth',2);
xlabel('Frequency (GHz)');
ylabel('Phase Imbalance (Degrees)');
title('Phase Imbalance Between S21 and S31');

%% -------------------------------
% 4. POWER SPLIT VERIFICATION
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
% 5. RETURN LOSS ONLY (S11)
% -------------------------------
figure;
% plot(freq/1e9, S11_dB, 'Color',[0.8500 0.3250 0.0980]);
plot(freq/1e9, S11_dB,'LineWidth',2);

grid on;
xlabel('Frequency (GHz)');
ylabel('Return Loss (dB)');
title('Return Loss S11 of Rat Race Coupler');

%% -------------------------------
% 6. ISOLATION PLOT: S32
% -------------------------------
figure;
plot(freq/1e9, S32_dB, 'g--','LineWidth',2);
grid on;
xlabel('Frequency (GHz)');
ylabel('Isolation (dB)');
title('RRC Isolation S32');


%% -------------------------------
% 7. SMITH PLOT : S11
% -------------------------------%

figure;
smithplot(S11);
title('Smith Chart - BLC S11');

