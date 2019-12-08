function [ routes, transfers, start, finish, duration ] = UseX( x, duration )
% Formatting output
routes = {};
for r = 1:Nroutes
    for i = 1:x(r)
        routes = [routes; routeNames(1)];
    end
end
for s = 1:Nstations
    if x(Nroutes + s) == 1
        start = stations{s};
    end
    if x(Nroutes + Nstations + s) == 1
        finish = stations{s};
    end
end
transfers = {};
for t = 1:Ntransfers
    for i = 1:x(Nroutes + 2*Nstations + t)
        transfers = [transfers; TRinfo(t, 1:3)];
    end
end
end