%% Exercise 1 : Basic matrix computation 
%%% Access to the elements of a vector/matrix

% Defining the elements : 
a=5
b=[2,4,4,5]
c=[3;8;3]
D=[2 9 8;7 6 4;7 4 3;2 4 0]

% Using the function size() : 
dima = size(a)
dimb = size(b)
dimc = size(c)
dimD = size(D)

% The difference between b and c is that, for b, there are one row and four
% columns, and the differents columns are separated by ",". For c, there
% are three rows and one columns, and the differents rows are separated by
% ";".

%Locations in matrix : 
ans1 = D(3,2)
ans2 = b(1,3)

%%% Operations on a vector

% Defining the vector : 
x=[2,6,6,5,4]

% Adding 5 to each element of x :
X2 = x+5

% Adding 6 to each odd index element :
X3 = x
X3(1:2:end) = X3(1:2:end) + 6

% Calculating the cube of each element :
X4 = x.^3

% Calculating the square root of each element :
X5 = sqrt(x)

%%% Operation between 2 vectors

% Creating vectors : 
X = [3,6,2,8]

Y = [4,5,3,1]

% Saving the results : 
xx = X'

yy = Y'

% Adding the sum of all elements of x to each element of y :
sum = xx+yy

% Dividing each element of x by each element of y :
div = xx./yy

% Multiplying (indexed) x by y and calling the result z :
z1 = (X.*yy)
z2 = (Y.*xx)

% Summing the elements of z and calling the result w :
w = z1 + z2

% Calculating Y ·x−w :
Z = Y*xx-w


%% Exercise 2 : Data reading and plotting

%Reading and loading the data : 
data = readtable("MeanGlobalT.txt")

% Calculating the min and max temprature : 
minT = min(data.temp)
maxT = max(data.temp)

% Plot the evolution of the global annual mean temperatures with respect to
% the year :
%plot(data, "year","temp")
plot(data.year, data.temp)
title('Évolution de la température moyenne globale (1901-2022)') % Add title
xlabel("Year") % Add a x axis label
ylabel("Temperature (°C)") % ADd a y axis label

% Calculating the average : 
temps = data.temp
years = data.year

% Calculating means : 
mean_all = mean(temps)
mean_1953_1982 = mean(temps(years >= 1953 & years <= 1982))
mean_2003_2022 = mean(temps(years >= 2003 & years <= 2022))

% Calculating the increase in temperature between the two periods :
temp_increase = mean_2003_2022 - mean_1953_1982

%% Exercise 3 : Integrals

% Defining constants : 

cp = 3850; %  the oceanic specific heat capacity J/(kg°C).
rho_sw = 1025; % The ocean water density kg/m^3.
z1 = -1000; % Initial depth in meters.
z2 = -600; % Final depth in meters.

% The temperature fonction : 

T = @(z) 1.125e-4 * (z + 1000).^2 + 4;

% Calculing the analytical solution of OHC between -1000 to -600 meters : 

integral_analytique = integral(T, z1, z2);
OHC_analytique = cp * rho_sw * integral_analytique;

% Parameters for numerical solution : 

n = 1000; % Number of segments for the Riemann sum
dz = (z2 - z1) / n; % Size of each segment
integral_numerique = 0; % Initialisation of the integral

% for loops : 

for i = 1:n
    z_current = z1 + (i - 1) * dz; % Current position on the z axis
    T_current = T(z_current); % Evaluation of T(z) at the current point
    integral_numerique = integral_numerique + T_current * dz; % Add the area of the rectangle
end

OHC_numerique = cp * rho_sw * integral_numerique;

% Comparison of results : 

% The analytical result is 1.5785e+10.
% The numerical result is 1.5771e+10.
% To improve the numerical result, we can increase the number of segments,
% because this reduces the size of the intervals and therefore increases
% precision.

% Calculing total OHC between -1000 and -600 : 

surface_oceans = 360e12; % Ocean's surface (m²)
OHC_total = OHC_numerique * surface_oceans;
%Response : 5.66775e+24.


