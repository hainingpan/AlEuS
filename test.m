

f(1:10) = parallel.FevalFuture;
for idx = 1:10
    f(idx) = parfeval(@func,2,idx);
end

magicResults = cell(1,10);
for idx = 1:10
    [completedIdx,value,v2] = fetchNext(f);
    magicResults{completedIdx} = value(1);
    fprintf('Got result with index: %d.\n', completedIdx);
end

function [r1,r2]=func(i)
r1=i;
r2=eye(i);
end