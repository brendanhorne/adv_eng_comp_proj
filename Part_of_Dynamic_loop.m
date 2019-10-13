clc,clear all, close all

%Open&run earthquake data
Earthquake_data

%Initial data
h=Time(end)/0.02;
Dis=0;
Vel=0;
Acc=Acceleration(1,1);
p=0;

%Dynamic loop Begins
for i=0:size(Acceleration,1)

%Effective Load Vector p
p=p+C*((2*Dis/h)+Vel)+M*((4*Dis/(h*h))+4*Vel/h+Acc);

%Solve for displacement (DisNew). Still req


%Update velocity 
Vel=Vel+(h/2)*(Acc+Acceleration(i,1));

%Update acceleration
Acc=(4/h*h)*(DisNew-Dis)-4/(h*Vel-Acc);

Dis=DisNew;
end


