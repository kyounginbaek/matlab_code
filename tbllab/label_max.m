function [label_result] = label_max(imdata)

[l,n]=bwlabel(imdata);

max_points=find(l==1);

for index = 2:n
    
    points = find(l==index);
    
    if(length(points) > length(max_points))
        max_points = points;

    end
    
end


[r,c]=size(imdata);
label_result=zeros(r,c);

for index = 1:length(max_points)
    
    label_result(max_points(index)) = 1;

end

end