function [ x, duration ] = GetX( input, variation )
% Code to solve the Subway Challenge (choice of 3 variations)
%   Given input spreadsheet which defines a subway system

%Data Input
stationData = readtable(input, 'Sheet', 'Stations', 'ReadRowNames', true);
stations = stationData.Name; % a cell array
Nstations = length(stations);

segmentData = readtable(input, 'Sheet', 'Segments', 'ReadRowNames', true);
segments = segmentData.Name;
Nsegments = length(segments);

routeData = readtable(input, 'Sheet', 'Routes', 'ReadRowNames', true);
routeNames = routeData.Name; % a cell array
Nroutes = length(routeNames);
routeData.ind = (1:Nroutes)';

transferData = cell(Nstations, 1);
TRinfo = {};
TRin = 0;
RoutesIn = cell(1,Nstations);
for s=1:Nstations
    trans = readtable(input, 'Sheet', stations{s}, 'ReadRowNames', true);
    %transferData{s,1} = trans;
    RoutesIn(s) = {trans.Row};
    [m,n] = size(trans);
    if m~= stationData{s,'Nin'} || n~=stationData{s,'Nout'}
        disp('Dimensions of transfer data not right')
    end
    for i = 1:m
        for j = 1:n
            if trans{i,j} < 1000
                TRin = TRin + 1;
                TRinfo(TRin, :) = {stations{s}, trans.Row{i}, stationData{s,2+j}{:}, trans{i,j}}; 
            end
        end
    end
end
Ntransfers = TRin; 

% Set up model. Need one variable for each route, and two for each station, and one for each transfer.
% Index x by [routes, stationsStart, stationsEnd, transfers]
cushionRows = 1500; %May need to increase?
Ain = zeros(cushionRows, Nroutes + 2*Nstations + Ntransfers);
bin = zeros(cushionRows, 1);
rin = 0;
% Constraint 1: Enters and Exits each station equally, with slack for
% where the attempt starts
for s = 1:Nstations
    rin = rin+1;
    % Routes that go into the station
    for i = 1:stationData{s,'Nin'}
        routeName = RoutesIn{1,s}{i};
        Ain(rin, routeData{routeName, 'ind'}) = -1;
    end
    % Routes that go out of the station
    for i = 1:stationData{s,'Nout'}
        Ain(rin, routeData{stationData{s,i+2},'ind'}) = 1;
    end
    Ain(rin, Nroutes + s) = -1;
    bin(rin, 1) = 0;
end

% Constraint 2: Enters and Exits each station equally, with slack for
% where the attempt ends
for s = 1:Nstations
    rin = rin+1;
    % Routes that go into the station
    for i = 1:stationData{s,'Nin'}
        routeName = RoutesIn{1,s}{i};
        Ain(rin, routeData{routeName, 'ind'}) = 1;
    end
    % Routes that go out of the station
    for i = 1:stationData{s,'Nout'}
        Ain(rin, routeData{stationData{s,i+2},'ind'}) = -1;
    end
    Ain(rin, Nroutes +Nstations + s) = -1;
    bin(rin, 1) = 0;
end

% Constraint 3: Can start the attempt at only one station
rin = rin+1;
for s = 1:Nstations
    Ain(rin, Nroutes + s) = 1;
end
bin(rin, 1) = 1;

% Constraint 4: Can end the attempt at only one station
rin = rin+1;
for s = 1:Nstations
    Ain(rin, Nroutes + Nstations + s) = 1;
end
bin(rin, 1) = 1;

% Constraint 5: Counts transfers made
for i = 1:Ntransfers
    rin = rin+1;
    Ain(rin, Nroutes + 2*Nstations + i) = -2;
    Ain(rin, routeData{TRinfo{i,2},'ind'}) = 1;
    Ain(rin, routeData{TRinfo{i,3},'ind'}) = 1;
    bin(rin, 1) = 1;
end

% Constraint 6: Rule constraint that is different based on the variation
switch variation
    case "A"
        for m = 1:Nsegments
            rin = rin +1;
            for r = 1:Nroutes
                if strcmp(routeData{r,'Segment'}{:}, segments{m})
                    Ain(rin, routeData{r, 'ind'}) = -1;
                end
            end
            bin(rin, 1) = -1;
        end
    case "B"
        for m = 1:Nsegments
            if segmentData{m,'Necc'}==1
                rin = rin +1;
                for r = 1:Nroutes
                    if strcmp(routeData{r,'Segment'}{:}, segments{m})
                        if routeData{r, 'Express'}==0
                            Ain(rin, routeData{r, 'ind'}) = -1;
                        end
                    end
                end
                bin(rin, 1) = -1;
            end
        end
    case "C"
        for m = 1:Nsegments
            if segmentData{m,'Necc'}==1
                rin = rin +1;
                for r = 1:Nroutes
                    if strcmp(routeData{r,'Segment'}{:}, segments{m})
                        Ain(rin, routeData{r, 'ind'}) = -1;
                        
                    end
                end
                bin(rin, 1) = -1;
            end
        end
end


Ain = Ain(1:rin, :);
bin = bin(1:rin, :);
transTimes = cell2mat(TRinfo(:,4))';
routeTimes = table2array(routeData(:,'Duration'))';
f = [routeTimes, zeros(1,2*Nstations), transTimes];
intcon = 1:(Nroutes + 2*Nstations + Ntransfers);
lb = zeros(1,Nroutes + 2*Nstations + Ntransfers);
ub = [inf(1,Nroutes), ones(1,2*Nstations), inf(1,Ntransfers)];

[x,duration] = intlinprog(f, intcon, Ain, bin, [],[],lb,ub);
disp('Solving is Done')
end