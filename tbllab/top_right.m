
function[result,result_point] = top_right(imdata)

[r,c] = size(imdata);
result=zeros(r,c);
result_point=0;

for index = 1:c
    for jndex = 1:r
        if(imdata(jndex,c-index+1)==1)
            result(jndex,c-index+1)=1;
            result_point=[jndex,c-index+1];
            break;
        end
    end
    if(result_point~=0)
        break;
    end
end

if(result_point(2)~=c)
    for index=result_point(2):c
        result(result_point(1),index)=1;
    end
end

end
