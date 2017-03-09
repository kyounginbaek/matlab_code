function[result,result_points]=top_boundary(imdata)

[r,c] = size(imdata);
result=zeros(r,c);
result_points=zeros(1,c);

for index = 1:c
    for jndex = 1:r
        if(imdata(jndex,index)==1)
            result(jndex,index)=1;
            result_points(1,index)=jndex;
            break;
        end
    end
end

end