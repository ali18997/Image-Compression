function [dctImg] =  dct8x8(imgIn)
imgIn = double(imgIn);
imgIn = imgIn - 128;
dctImg = zeros(8,8);
for v = 0:7
    for u = 0:7
        
        if u == 0
            Cu = 1/sqrt(2);
        else
            Cu = 1;
        end
        
        if v == 0
            Cv = 1/sqrt(2);
        else
            Cv = 1;
        end
       
        sum = 0;
        for y = 0:7
            for x = 0:7
                imgVal = imgIn(y+1,x+1);
               sum = sum + imgVal*cos(v*pi*(2*y+1)/(16))*cos(u*pi*(2*x+1)/(16));  
            end
        end
        dctImg(v+1, u+1) = (1/4)*Cv*Cu*sum;     


    end
end


end
             