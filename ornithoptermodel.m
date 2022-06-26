clear
clc
clf


 %% 
    discreteBins = 1000; %We will use this number of bins for plotting and calculating all functions such as Torque, speed etc.
% Input part of the main function
 
    StallTorque = 110;
    StallCurrent = 400;
    RatedVoltage = 6;
    NoLoadCurrent = 35;
    NoLoadSpeed = 2500;
    
 %%
    %Some basic input error checking is here.
    if or(not(isfloat(StallTorque)), isempty(StallTorque))
        StallTorque = 17;
        fprintf('\nUsing default value for StallTorque');
    end
    if or(not(isfloat(StallCurrent)), isempty(StallCurrent))
        StallCurrent = 700;
        fprintf('\nUsing default value for StallCurrent');
    end
    if or(not(isfloat(RatedVoltage)), isempty(RatedVoltage))
        RatedVoltage = 6;
        fprintf('\nUsing default value for RatedVoltage');
    end
    if or(not(isfloat(NoLoadCurrent)), isempty(NoLoadCurrent))
        NoLoadCurrent = 40;
        fprintf('\nUsing default value for NoLoadCurrent');
    end
    if or(not(isfloat(NoLoadSpeed)), isempty(NoLoadSpeed))
        NoLoadSpeed = 290;
        fprintf('\nUsing default value for NoLoadSpeed');
    end
    %
   

 %% 
    %Here we calculate basic stuff to get all the variables and outputs.
    Resistance = RatedVoltage / (StallCurrent/1000);
    %Torque line
    TorqueLine = 0:(StallTorque/discreteBins):StallTorque;
    %Current Line
    CurrentLine = NoLoadCurrent:(StallCurrent-NoLoadCurrent)/discreteBins:StallCurrent;
    %Speed Line
    SpeedLine = NoLoadSpeed: (0-NoLoadSpeed)/discreteBins : 0;
    % Torque Constant in Torque per current is
    SlopeOfTorqueVsCurrent = (StallCurrent - NoLoadCurrent) / (StallTorque);

    %Output Mechanical Power in watts is Torque * Speed * 0.00074 watts
    OutputPower = ((TorqueLine* 0.0000980665) .* SpeedLine) ./ 9.5488;
    %Input Electrical Power to the motor is Voltage * Current
    InputPower = CurrentLine * RatedVoltage / 1000; %We are dividing by 1000 as the input was in mA and we need power in Watts.
 
    TorqueLine_Nm = TorqueLine *0.0000980665;


%% 

%Plot part of the functions 
    figure(1)
    subplot(3,1,1)
    [hAx, hLine1, hLine2] = plotyy([0 StallTorque], [NoLoadSpeed 0], [0 StallTorque], [NoLoadCurrent StallCurrent]); %This is the TorqueLoad vs. Motor Speed graph
    
    title('Torque vs Speed \& Torque vs Current','Interpreter','latex');
    xlabel('Torque [gram-cm]','Interpreter','latex');
    ylabel(hAx(1), 'Speed [RPM]','Interpreter','latex');
    ylabel(hAx(2), 'Current [mA]','Interpreter','latex');
    
    %This is the plot of the Output Mechanical power in watts vs. Input
    %Electrical power in Watts.
    subplot(3,1,2);
    
    [h2Ax, h2Line1, h2Line2] = plotyy(TorqueLine, OutputPower, TorqueLine, InputPower);
    xlabel('Torque [gram-cm]','Interpreter','latex');
    ylabel(h2Ax(1), 'OutputPower [Watt]','Interpreter','latex');
    ylabel(h2Ax(2), 'InputPower [Watt]','Interpreter','latex');
    title('Torque vs. Output Power \& Torque vs. Input Power','Interpreter','latex');
    
    %This is the plot of the Power Efficiency of the motor.
    subplot(3,1,3);
    PowerEff = OutputPower ./ InputPower;
    plot(TorqueLine, PowerEff);
    xlabel('Torque [gram-cm]','Interpreter','latex');
    ylabel('Power Efficiency [no unit]','Interpreter','latex');

    %Output information part of the function    
    fprintf('\n\n\nSlope of TorqueVsCurrent is %f. The recprocal is %f\n', SlopeOfTorqueVsCurrent, (1/SlopeOfTorqueVsCurrent));
    %Max Output Power is at
    [V,I] = max(OutputPower);
    fprintf('Maximum output mechanical power is %f(watts).\nThis happens at the Torque load of %f(N-m), with Current %f(mA)\n', OutputPower(I), TorqueLine_Nm(I), CurrentLine(I));
    fprintf('Resistance of the motor is %f (ohms)\n', Resistance);

%% 

% 
betamax = 51.86 .* 0.0174533;   %angular distance between open and closed (degrees)
thetamax = 46.73 .* 0.0174533;  %37.55
delta = 0;                      %angular height position
psi = 90;                       %angular direction of flight
b = 55*(10^-3);                 %length wing (m)
s = 1310*(10^-6);               %wing area (m^2)
m = 0.0035;                      %mass ornithopter (kg)
g = 9.81;                       %gravitational constant
rho = 1.225;                    %air density   
c = 28*(10^-3);                 %chord length
U = 1.508 * m^(1/6);            %velocity in which direction it is going
t = 0:(6/discreteBins):6;       %time period (6 seconds)
y = 0:(0.055/discreteBins):0.055;
B = (55*(10^-3)/2);

%% 

% f = (1.08./b)*(m./rho)^(1/3) * sqrt(g./sqrt(s))
fold = 0:(NoLoadSpeed)/discreteBins:NoLoadSpeed;
f = fold/60;

%% 
figure(2)
beta = betamax .* cos(2 .* pi .* f .* t);
subplot(3,1,1)
plot(f, beta);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Flapping Angle [Rad]','Interpreter','latex');

betadot = -2 .* pi .* f .* betamax .* sin(2 .* pi .* f .* t);
subplot(3,1,2)
plot(f, betadot);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Flapping Rate [Rad/Second]','Interpreter','latex');

betaddot = -4 .* pi.^2 .* f.^2 .* betamax .* cos(2 .* pi .* f .* t);
subplot(3,1,3)
plot(f, betaddot);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Flapping Acceleration [Rad/$Second^2$]','Interpreter','latex');

figure(3)
theta = thetamax .* cos(2 .* pi .* f .* t + psi);
subplot(3,1,1)
plot(f, theta);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Pitching Angle [Rad]','Interpreter','latex');

thetadot = -2 .* pi .* f .* (y./B) .* thetamax .* sin(2 .* pi .* f .* t);
subplot(3,1,2)
plot(f, thetadot);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Pitching Rate [Rad/Second]','Interpreter','latex');

thetaddot = -4 * pi.^2 * f.^2 .* (y./B) .* thetamax .* cos(2 .* pi .* f .* t);
subplot(3,1,3)
plot(f, thetaddot);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Pitching Acceleration [Rad/$Second^2$]','Interpreter','latex');

%% 

VoltageLine = 0:(InputPower./(CurrentLine*0.001))./discreteBins:(InputPower./(CurrentLine*0.001));

F = TorqueLine_Nm ./ (b*4);

diff_y = 5.5*10^-5;

N  = 2* abs(-(rho .* pi .* c^2)/4 .* (thetadot .* U + y .* betaddot .* cos(theta) - 0.5 .* thetaddot)) .* diff_y .* f;

Power = (1/9.5488) .* TorqueLine_Nm .* f * 60;

F_ver = F*cos(delta) + N.*cos(-theta).*cos(beta).*cos(theta);

Voltage = Power./CurrentLine*1000;

F_Fit_Voltage = fit(Voltage(:), F_ver(:), 'poly2');
F_Fit_Power = fit(Power(:), F_ver(:), 'poly2');

%% 

figure(4)
subplot(3,1,1)
plot(Power, F);
xlabel('Power [W]','Interpreter','latex');
ylabel('Lift from Torque [N]','Interpreter','latex');

subplot(3,1,2);
plot(VoltageLine, F);
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Lift [N]','Interpreter','latex');

subplot(3,1,3)
plot(VoltageLine, f);
xlabel('Voltage [V]','Interpreter','latex');
ylabel('RPM','Interpreter','latex');

figure(5)
subplot(3,1,1)
plot(t, N);
xlabel('Time [s]','Interpreter','latex');
ylabel('Apparent Mass Effect [N]','Interpreter','latex');

subplot(3,1,2)
plot(f, N);
xlabel('Frequency [Hz]','Interpreter','latex');
ylabel('Apparent Mass Effect [N]','Interpreter','latex');

subplot(3,1,3)
plot(Power, N);
xlabel('Power [W]','Interpreter','latex');
ylabel('Apparent Mass Effect [N]','Interpreter','latex');

figure(6)
subplot(3,1,1)
plot(Power, f);
xlabel('Power [W]','Interpreter','latex');
ylabel('Frequency [Hz]','Interpreter','latex');

subplot(3,1,2)
plot(Voltage, N);
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Apparent Mass Effect [N]','Interpreter','latex');

subplot(3,1,3)
plot(F_Fit_Voltage, Voltage, F_ver);
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Vertical Force [N]','Interpreter','latex');

figure(7)
subplot(3,1,1)
plot(F_Fit_Power, Power, F_ver);
legend('Data', 'Fitting Line', 'Interpreter','latex')
xlabel('Power [W]','Interpreter','latex');
ylabel('Vertical Force [N]','Interpreter','latex');

subplot(3,1,2)
plot(Power, TorqueLine_Nm);
xlabel('Power [W]','Interpreter','latex');
ylabel('Torque [N-m]','Interpreter','latex');

