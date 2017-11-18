%% Purdue Orbital 
% 
% HELIUM GAS LIFTING MODEL V2.0
%
% Written by Mathew Powell, Drew Sherman, Ethan Wahl, and Dan Qiao
% Started 10/25/2016
% Updated 11/11/2017

%% Description of model
% 
% Uses the basic buoyancy force equation to find force on a payload lifted
% by helium gas 
%
% Time-steps through force, velocity, and position equations to find values
% 
% Plots acceleration, velocity, and position vs time
% Plots volume vs altitude
% 
% ---ASSUMPTIONS---
%   ideal gas
%   onion shaped balloon
%   standard atmosphere as used in AAE 251 at Purdue (see
%       Standard_Atmosphere.m)
%   Temperature of helium is approximately the same as the ambient air at
%   initial conditiaon.
%   No loss of helium during ascent (constant mass)
%   Mass of entire system constant during ascent
%       Mass of helium and balloon needs to be included in weight fighting
%       the buoyancy force - recursive, may be possible to estimate
%   Helium is initially @ STP

%3,000 m^3 with a membrane thickness of 80 micrometers
%100 kg payload
% kg mass of balloon

%% Helium Lifting Model V2.0
%

function [Performance] = Zero_model(payload_mass, volume, mass_balloon, mass_helium, plot_suppression)
% inputs :
%   payload_mass:       wanted payload mass in kg
%   volume:             Volume of the Balloon (m^3)
%   mass_balloon:       Mass of balloon in kg
%   mass_helium:        Mass of helium in kg
%   plot_suppression:   0: no plots, anything else: plots


clc
close all

% set arrays for use later
Fb_array  = [0];                        % buoyancy force [N]
accel = [];                             % acceleration seen by system [m/s^2]
velocity_array = [0];                   % velocity of the system [m/s]
system_alt = [0];                       % vertical position of system [m]
balloon_vol =[];                        % volume of helium [m^3]

delta_t = 1;                            % Change in time for each iteration of the loop [s]
Cd = 0.3;                               %default drag coefficient for spherical balloons
Tg = [];                                %Gas temp inside balloon (assumed to be Tatm initially)
P_a = [];                               %air pressure at altitude
    
R_He = 2077;                            % Gas constant of helium J/Kg 
Cp_he = 5.300*1000;                     % Specific heat of helium J/Kg*C
g = 9.81;                               % [m/s^2]
G = 6.67408 * 10^-11;                   % Gravitation Constant
mol_mass = .004002602;                  % Molar mass of heliumkg/mol
R = 8.314;                              % Gas constant of ideal gases [m^3 Pa]/[K mol]
cylinder_coef = 1.233;                  % kg.(He)/Cylinder(244 ft.^3)

M_earth = 5.972 * 10^24;                % mass of earth [kg]
R_earth = 6371000;                      % radius of earth [m]

mol_mass_atmo = .029;                   % kg/mol 
mol_ratio = mol_mass_atmo / mol_mass;   %Molar mass ratio [N/A]

L_array=[0];
D_array = [0];


%% Find intial conditions (namely He volume)
diameter = 1.383 * volume^(1/3);    %Use diameter equation from Helium Ascent 3D paper

Ab = pi * (diameter(end) / 2) ^ 2;                  %Cross-sectional Area of the balloon

Mtot = payload_mass + mass_balloon + mass_helium;   %Total Mass of the whole system

num_moles = mass_helium / mol_mass;   %mol ratio (constant, super pressure balloon)

num_cylinders = mass_helium / cylinder_coef;

volume_helium = num_cylinders * 244;

%% Loop Through Kinetics
% begin time_stepping 
%while (velocity(end)>=0 && system_alt(end) < final_alt && system_alt(end) >= 0 )
i = 1;
velocity(1) = .01;
t = 0;
diameter = 0;

 while (system_alt(end) < 100000 && t(end) < 1000)   
    [~, P_air, rho_air, T_air] = Standard_Atmosphere(system_alt(end), 0);       %Standard Atmosphere function from AAE 251
    P_a = [P_a, P_air];
    
    %Update g
    if i ~= 1
        g = G * (M_earth) / (system_alt(end) + R_earth)^2;
    end
    
    %Set initial gas conidtion and other gas calculation
    if i == 1
        Tg = T_air;
        temp_ratio = 1;
        dz = 0;
    else
        c = 1/(num_moles * R);
        q = 0;                      %assume zero heat flux
        
        temp_ratio = Tg(i-1) / T_air;
        
        dz = velocity(end) * delta_t;
        
        Tg_temp = Tg(i-1) - dz * g * mol_ratio * temp_ratio / Cp_he;

        Tg = [Tg, Tg_temp];        
    end
    
    %Determining the change in mass of the system
    if i ~= 1
        delta_m = find_mass(system_alt(end - 1), system_alt(end), mass_helium)  %Seems to work but I have not guarenteed 
    
        mass_helium = mass_helium - delta_m;
    end
    
    %Find the Boyant force
    Fb = (mass_helium * mol_ratio * temp_ratio - Mtot) * g; %reference eq 2.94 from Engineering Fundamentals of balloons combines bouyant and gravity
    Fb_array = [Fb_array,Fb];
    
    D = Cd * Ab * .5 * rho_air * velocity^2;            %Force fo drag
    
    L = Fb - D;
    
    L_array = [L_array,L];
    D_array = [D_array,D];
    
    %find acceleration
    acc_temp = L / Mtot;
    
    accel = [accel, acc_temp];
     
    %find velocity
    velocity = sqrt(velocity^2 + 2 * acc_temp * delta_t);  %Find new way to obtain velocity (needs ability to go negative)
    
    velocity_array = [velocity_array, velocity];
    system_alt = [system_alt, system_alt(end) + velocity * delta_t + (acc_temp * delta_t^2 / 2)];
 
    t = [t, t(end)+delta_t];
    i = i + 1;

 end
%% Plotting
if (plot_suppression ~= 0) 
    hold on

    subplot(3,2,1)
    plot(t, system_alt);
    xlabel('Time (s)');
    ylabel('Altitude (m)');
    title('Altitude vs Time for given Payload Mass')
   
    subplot(3,2,2)
    plot(t(2:end), accel);grid
    xlabel('Time (s)');
    ylabel('Acceleration (m/s^2)');
    title('Acceleration vs Time for given Payload Mass')

    subplot(3,2,3)
    plot(t, velocity_array);grid
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Velocity vs Time for given Payload Mass')
    
    subplot(3,2,4)
    plot(t, L_array);
    xlabel('time(s)');
    ylabel('Lift(F)');
    title('Lift vs Time')
    
    subplot(3,2,5)
    plot(t, Fb_array)
    xlabel('time(s)');
    ylabel('Fb(F)');
    title('Fb vs Time')
    
    subplot(3,2,6)
    plot(t, D_array)
    xlabel('Time (s)')
    ylabel('Drag (N)')
    title('Drag vs Time')

end

fprintf('The total number of helium cansiters needed are: %f\n', num_cylinders);
fprintf('The total number of moles of helium needed are: %f\n', num_moles);
fprintf('The final altitude is : %f\n', system_alt(end));

