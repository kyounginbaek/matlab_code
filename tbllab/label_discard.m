function [label_result]=label_discard(imdata,std_size)

[l,n] = bwlabel(imdata);

for a = 1:n
    points = find(l==a);
    if(length(points) < std_size)
        for index = 1:length(points)
            imdata(points(index))=0;
        end
    end
end

label_result=imdata;

end