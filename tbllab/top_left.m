function[result,result_point] = top_left(imdata)

[r,c] = size(imdata);
result=zeros(r,c);
result_point=0;

for index = 1:c
    for jndex = 1:r
        if(imdata(jndex,index)==1)
            result(jndex,index)=1;
            result_point=[jndex,index];
            break;
        end
    end
    if(result_point~=0)
        break;
    end
end

if(result_point(2)~=1)
    for index=1:result_point(2)
        result(result_point(1),index)=1;
    end
end

end
