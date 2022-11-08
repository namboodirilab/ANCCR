function d = point2line(pt, v1, v2)
if size(pt,2)==2
    pt = [pt,zeros(size(pt,1),1)];
    v1 = [v1,0];
    v2 = [v2,0];
end
a = v1 - v2;

d = NaN(size(pt,1),1);
for i = 1:size(pt,1)
    b = pt(i,:) - v2;
    d(i) = norm(cross(a,b)) / norm(a);
end
end
