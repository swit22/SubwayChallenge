function [ solution, transfers, duration ] = SubChallengeCycle( input, variation )
% Code to solve the Subway Challenge (choice of 3 variations)
%   Given input spreadsheet which defines a subway system
%   NOTE: Transfers constraint not working, so transfer will output empty

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
                TRinfo(TRin, :) = {stations{s}, trans.Row{i},stationData{s,3+j}{:}, trans{i,j}}; 
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

% Set up model. Need one variable for each route, and one for each transfer.
% Index x by [routes, routesStart, routesEnd, transfers]
cushionRows = 1500; %May need to increase?
Ain = zeros(cushionRows, 1*Nroutes + Ntransfers); %only 1 route variable now
bin = zeros(cushionRows, 1);
rin = 0;

Aeq = zeros(cushionRows, 1*Nroutes + Ntransfers); %only 1 route variable now
beq = zeros(cushionRows, 1);
req = 0;
% Constraint 1: Rides and Transfers from each route equally
for r = 1:Nroutes
    RinName = routeNames{r};
    endStat = routeData{r, 'End'}{:};
    req = req+1;
    Aeq(req, r) = 1;
    
    % Routes that go out of the station
    for i = 1:stationData{endStat,'Nout'}
        RoutName = stationData{endStat, 3+i};
        transIndex = TRinfo{strcmp(TRinfo.In, RinName) & strcmp(TRinfo.Out,RoutName),'ind'};
        Aeq(req, 1*Nroutes + transIndex) = -1;  %only 1 route variable now
    end
    %Aeq(req, 1*Nroutes + r) = -1; %no need for end point
    beq(req, 1) = 0;
end

% Constraint 2: Transfers into and rides each route equally
for r = 1:Nroutes
    RoutName = routeNames{r};
    startStat = routeData{r, 'Start'}{:};
    req = req+1;
    Aeq(req, r) = 1;
    
    % Routes that go into the station
    for i = 1:stationData{startStat,'Nin'}
        RinName = RoutesIn{stationData{startStat,'ind'}}{i};
        transIndex = TRinfo{strcmp(TRinfo.In, RinName) & strcmp(TRinfo.Out,RoutName),'ind'};
        Aeq(req, 1*Nroutes + transIndex) = -1;  %only 1 route variable now
    end
    %Aeq(req, 1*Nroutes + r) = -1; %no need for start point
    beq(req, 1) = 0;
end

%no need for start point
% % Constraint 3: Can start the attempt at only one station 
% rin = rin+1;
% for r = 1:Nroutes
%     Ain(rin, 1*Nroutes + r) = 1;
% end
% bin(rin, 1) = 1;

%no need for end point
% % Constraint 4: Can end the attempt at only one station
% rin = rin+1;
% for r = 1:Nroutes
%     Ain(rin, 2*Nroutes + r) = 1;
% end
% bin(rin, 1) = 1;

% Constraint 5: Rule constraint that is different based on the variation
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
size(Ain)
size(Aeq)
transTimes = table2array(TRinfo(:,'Dur'))';
routeTimes = table2array(routeData(:,'Duration'))';
f = [routeTimes, transTimes]; %only 1 route variable now
intcon = 1:(1*Nroutes + Ntransfers); %only 1 route variable now
lb = zeros(1,1*Nroutes + Ntransfers); %only 1 route variable now
ub = [inf(1,Nroutes), inf(1,Ntransfers)]; %only 1 route variable now

[x,duration] = intlinprog(f, intcon, Ain, bin, Aeq,beq,lb,ub);
x = round(x);

% Formatting output
routes = {};
for r = 1:Nroutes
    for i = 1:x(r)
        routes = [routes; routeNames(r)];
    end
end
% for r = 1:Nroutes
%     if x(Nroutes + r) == 1
%         startR = routeNames{r};
%         start = routeData{r,'Start'};
%         finishR = routeNames{r};
%         finish = routeData{r,'End'};
%     end
% end
transfers = {};
maxCost = 0; %new for finding max transfer cost transfer we take
maxT = 1; %new for finding max transfer cost transfer we take
for t = 1:Ntransfers
    for i = 1:x(1*Nroutes + t)
        transfers = [transfers; TRinfo(t, 1:4)];
        %find the highest cost transfer we take while we're at it
        if maxCost < TRinfo{t,'Dur'} %new for finding max transfer cost transfer we take
            maxT = t; %new for finding max transfer cost transfer we take
            maxCost = TRinfo{t,'Dur'}; %new for finding max transfer cost transfer we take
        end
    end
end

%find highest cost transfer that we take. That hub will be our start and
%end point, because it will save us the most time
startR = routeNames{maxT}; %new for finding max transfer cost transfer we take

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