%%
%Greg Kojadinovich
%Peng Robinson EOS
%Z-Factor for Methane

clear;
close all;

%% Variable Library
omegabeta = 0.07780; %Co-Volume Parameter Constant
omegaalpha = 0.45724; %Attraction Parameter Coeffificent
R = 10.73; %gas coefficient psia*cuft / lbmol*R
Tc = 343.01; %critical temperature Rankine
Pc = 667.0; %psia
omega = 0.0115; %accentric factor
func_omega = 0.374640+1.54226*omega-0.26992*omega^2; %f(omega)
Tr = [1.05 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3.0]; #pseudo-reduced temperatures
Pr = [0:0.05:20]'; #pseudo-reduced pressures


%% Looping for each temperature and pressure
for i = 1:16 # for each fixed isotherm
    for j = 1:401 # for each pressure

        % Calculcate Temperature and Pressure
        T = Tr(i)*Tc;
        P = Pr(j)*Pc;

        %Calculate coefficent b
        b = (omegabeta*(R*Tc))/Pc;

        %Calculate coefficent a
        a_alpha = omegaalpha*((Tc^2)*(R^2)/Pc)*(1+func_omega*(1-Tr(i)^0.5))^2; %alpha

        %Pure Substane Dimensionless Attraction Parameters A
        A = (a_alpha*P)/((R^2)*(T^2));

        %Covolume Parameter B
        B = (b*P)/(R*T);

        %% Caculate Peng-Robinson Cubic Coefficients a1,b1,c1
        a1 = -(1-B);
        b1 = A-3*B^2-2*B;
        c1 = -(A*B-B^2-B^3);

        %% Solve Compressibility
        pp = [1 a1 b1 c1]; %polynomials
        z = roots(pp); %root of poly for z

        #Find real roots
        index = find(imag(z) ==0);

        #Saving real and largest root
        if length(index)>1
          Zreal = real(z(index));
          Zfactors(i,j) = [max(Zreal)];
        else
          Zfactors(i,j) = real(z(index));
        end
  end
end

%% Plot Compressiblity Vs Pr
figure(1)
plot(Pr,Zfactors);
xlabel('Pseudo Reduced Pressure, Pr')
ylabel('Compressibiliy  Factor, Z')
title('Compressibility Factor For Natural Gas')
legend('1.05', '1.1', '1.2','1.3', '1.4' ,'1.5','1.6','1.7','1.8','1.9','2.0','2.2','2.4','2.6', '2.8','3.0');

#Calculate temperature and pressure for graph
T = Tr*Tc; # temperature in Rankine
P = Pr*Pc; # pressure in PSIA


figure(2)
plot(P,Zfactors);
xlabel('Pressure, PSIA')
ylabel('Compressibiliy  Factor, Z')
title('Compressibility Factor For Natural Gas')
legend('360.16',	'377.31',	'411.61',	'445.91',	'480.21',	'514.51',	'548.82',	'583.12',	'617.42',	'651.72',	'686.02',	'754.62',	'823.22',	'891.83',	'960.43',	'1029');
