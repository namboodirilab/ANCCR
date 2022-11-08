function eventlog = joinEventlogs(varargin)
%JOINEVENTLOGS: join multiple eventlogs and generate a single event log
eventlog = [];
nt = 0;
for i = 1:length(varargin)
    eventlog_temp = varargin{i};
    eventlog_temp(:,2) = eventlog_temp(:,2)+nt;
    eventlog = [eventlog;eventlog_temp];
    nt = eventlog(end,2);
end

end
