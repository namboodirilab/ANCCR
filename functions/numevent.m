function nevent = numevent(event,reftime,window,normalize)
%timeunit: ms

nevent = cellfun(@(x) sum(event>=x+window(1) & event<=x+window(2)),num2cell(reftime));
if normalize==1
   nevent = nevent/diff(window/1000); 
end