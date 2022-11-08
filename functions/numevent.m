function nevent = numevent(event,reftime,window,normalize)
%NUMEVENTS: calculate number of event during window from each reftime. If
%normalized is true, normalize nevent by length of window. Here, time unit
%is milisecond. 

nevent = cellfun(@(x) sum(event>=x+window(1) & event<=x+window(2)),num2cell(reftime));
if normalize==1
   nevent = nevent/diff(window/1000); 
end