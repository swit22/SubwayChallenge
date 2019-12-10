function [ solution, transfers, start, finish, duration ] = SubChallengeNoCycle( input, variation )
% Code to solve the Subway Challenge (choice of 3 variations)
%   Given input spreadsheet which defines a subway system

% Changes in this version of the code: we are attempting to implement
% polynomial asymmetric TSP subtour constraints. Each route now has two
% additional variables Ustart and Ufinish, representing the start and
% finish nodes. Then there will be the constraints that prevent subtours
% smaller than n (pick n). It is based on all arcs, so all transfers and
% routes are arcs.

%Data Input
stationData = readtable(input, 'Sheet', 'Stations', 'ReadRowNames', true);
stations = stationData.Name; % a cell array
Nstations = length(stations);
stationData.ind = (1:Nstations)';

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
    transferData{s,1} = trans;
    RoutesIn(s) = {trans.Row};
    [m,n] = size(trans);
    if m~= stationData{s,'Nin'} || n~=stationData{s,'Nout'}
        disp('Dimensions of transfer data not right')
    end
    for i = 1:m
        for j = 1:n
            if trans{i,j} < 1000
                TRin = TRin + 1;
                TRinfo(TRin, :) = {stations{s}, trans.Row{i},...
                    stationData{s,3+j}{:}, trans{i,j}}; 
            end
        end
    end
end
Ntransfers = TRin; 
TRinfo = cell2table(TRinfo);
TRinfo.Properties.VariableNames = {'Stat', 'In', 'Out', 'Dur'};
TRinfo.ind = (1:Ntransfers)';

% USE something that looks like TRinfo(strcmp(TRinfo.In,'RedMCSG'),:) in
% order to find rows of TRinfo that you want :0

% Set up model. Need three variables for each route, and one for each transfer.
% Index x by [routes, routesStart, routesEnd, transfers]
cushionRows = 1500; %May need to increase?
Ain = zeros(cushionRows, 5*Nroutes + Ntransfers);
bin = zeros(cushionRows, 1);
rin = 0;

Aeq = zeros(cushionRows, 5*Nroutes + Ntransfers);
beq = zeros(cushionRows, 1);
req = 0;
% Constraint 1: Rides and Transfers from each route equally, with slack for
% where the attempt finishes
for r = 1:Nroutes
    RinName = routeNames{r};
    endStat = routeData{r, 'End'}{:};
    req = req+1;
    Aeq(req, r) = 1;
    
    % Routes that go out of the station
    for i = 1:stationData{endStat,'Nout'}
        RoutName = stationData{endStat, 3+i};
        transIndex = TRinfo{strcmp(TRinfo.In, RinName) & strcmp(TRinfo.Out,RoutName),'ind'};
        Aeq(req, 5*Nroutes + transIndex) = -1;
    end
    Aeq(req, 2*Nroutes + r) = -1;
    beq(req, 1) = 0;
end

% Constraint 2: Transfers into and rides each route equally, with slack for
% where the attempt starts
for r = 1:Nroutes
    RoutName = routeNames{r};
    startStat = routeData{r, 'Start'}{:};
    req = req+1;
    Aeq(req, r) = 1;
    
    % Routes that go into of the station
    for i = 1:stationData{startStat,'Nin'}
        RinName = RoutesIn{stationData{startStat,'ind'}}{i};
        transIndex = TRinfo{strcmp(TRinfo.In, RinName) & strcmp(TRinfo.Out,RoutName),'ind'};
        Aeq(req, 5*Nroutes + transIndex) = -1;
    end
    Aeq(req, 1*Nroutes + r) = -1;
    beq(req, 1) = 0;
end


% Constraint 3: Can start the attempt at only one station 
rin = rin+1;
for r = 1:Nroutes
    Ain(rin, 1*Nroutes + r) = 1;
end
bin(rin, 1) = 1;

% Constraint 4: Can end the attempt at only one station
rin = rin+1;
for r = 1:Nroutes
    Ain(rin, 2*Nroutes + r) = 1;
end
bin(rin, 1) = 1;

% Constraint 5: Eliminating Subtours
% Start with all arcs represented by routes
n = Nroutes+Ntransfers; % can change n later if needed
% From 2 because that's what is in notes, trying to eliminate first
% start/end
for r = 2:Nroutes
    rin = rin +1;
    Ain(rin, r) = n;
    Ain(rin, 3*Nroutes + r) = 1;
    Ain(rin, 4*Nroutes + r) = -1;
    bin(rin, 1) = n-1;
end

% And then do the transfers
for t = 1:Ntransfers
    startInd = routeData{TRinfo{t,'In'},'ind'};
    endInd = routeData{TRinfo{t,'Out'},'ind'};
    if endInd >=2
        rin = rin + 1;
        Ain(5*Nroutes + t) = n;
        Ain(rin, 3*Nroutes + endInd) = -1;
        Ain(rin, 4*Nroutes + startInd) = 1;
        bin(rin, 1) = n-1;
    end
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
Aeq = Aeq(1:req, :);
beq = beq(1:req, :);
transTimes = table2array(TRinfo(:,'Dur'))';
routeTimes = table2array(routeData(:,'Duration'))';
f = [routeTimes, zeros(1,4*Nroutes), transTimes];
intcon = [1:3*Nroutes, 5*Nroutes +1: 5*Nroutes+ Ntransfers];
lb = zeros(1,5*Nroutes + Ntransfers);
ub = [inf(1,Nroutes), ones(1,2*Nroutes), inf(1,2*Nroutes + Ntransfers)];

[x,duration] = intlinprog(f, intcon, Ain, bin, Aeq,beq,lb,ub);
x = round(x);

% Formatting output
routes = {};
for r = 1:Nroutes
    for i = 1:x(r)
        routes = [routes; routeNames(r)];
    end
end
for r = 1:Nroutes
    if x(Nroutes + r) == 1
        startR = routeNames{r};
        start = routeData{r,'Start'};
    end
    if x(2*Nroutes + r) == 1
        finishR = routeNames{r};
        finish = routeData{r,'End'};
    end
end
transfers = {};
for t = 1:Ntransfers
    for i = 1:x(5*Nroutes + t)
        transfers = [transfers; TRinfo(t, 1:4)];
    end
end

transfers.done = zeros(size(transfers,1),1);
solution = {};
current = startR;
solInd = 0;
while solInd < length(routes)
    solution = [solution ; current];
    solInd = solInd + 1;

    next = transfers(strcmp(transfers.In, current) & transfers.done==0, :);
    if size(next,1) > 1
        returnto = current;
        current = next{1,'Out'};
        transfers{strcmp(transfers.In, returnto) & strcmp(transfers.Out, current), 'done'} = 1;
    end
    if size(next,1) == 1
        prev = current;
        current = next{1,'Out'};
        transfers{strcmp(transfers.In, prev) & strcmp(transfers.Out, current), 'done'} = 1;
    end
    if size(next,1) == 0
        if solInd < length(routes)
            current = returnto;
        else
            break;
        end
    end

end
    


end