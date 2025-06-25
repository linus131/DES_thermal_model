clc; clearvars; close all;
read_elements = csvread("elementfile.csv");
read_nodes = csvread("nodefile.csv");
element_temperatures = csvread("elem_temps.csv");

nodes = read_nodes(:,2:end);
elements = read_elements(:,2:end);

element_centers = zeros(size(elements, 1), 3);

for i = 1:size(elements,1)
    nds = elements(i, :);
    coords = nodes(nds, :);
    element_centers(i, :)  = mean(coords);
end

kdtree = KDTreeSearcher(element_centers);
% find closest element center to C2 thermocouple
c2coords = [198.65, -99, 186.72]*1e-3;
id = kdtree.knnsearch(c2coords);
id_temps = element_temperatures(id,:);
% filter later zeros
id_temps = id_temps(id_temps~=0);
plot(id_temps(1:2:end-1), id_temps(2:2:end)-273.15,'-x');
xlabel('times (s)')
ylabel('temperatures(^{\circ} C)')
hold on

% find closeset element center to C3 thermocouple
c3coords = [717.32 -4.83 191.38]*1e-3;
id = kdtree.knnsearch(c3coords);
id_temps = element_temperatures(id,:);
% filter later zeros
id_temps = id_temps(id_temps~=0);
%figure
plot(id_temps(1:2:end-1), id_temps(2:2:end)-273.15,'-x');
xlabel('times(s)')
ylabel('temperatures(^{\circ} C)')

% find closeset element center to C4 thermocouple
c4coords = [221.58 1.68 388.24]*1e-3;
id = kdtree.knnsearch(c4coords);
id_temps = element_temperatures(id,:);
% filter later zeros
id_temps = id_temps(id_temps~=0);
%figure
plot(id_temps(1:2:end-1), id_temps(2:2:end)-273.15,'-x');
xlabel('times(s)')
ylabel('temperatures(^{\circ} C)')

% find closeset element center to C5 thermocouple
c5coords = [665.33 -5.78 388.71]*1e-3;
id = kdtree.knnsearch(c5coords);
id_temps = element_temperatures(id,:);
% filter later zeros
id_temps = id_temps(id_temps~=0);
%figure
plot(id_temps(1:2:end-1), id_temps(2:2:end)-273.15,'-x');
xlabel('times(s)')
ylabel('temperatures(^{\circ} C)')

legend('c2','c3', 'c4', 'c5')