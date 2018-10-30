function [idctImg] =  idct8x8(imgIn)

idctImg = zeros(8,8);
for x = 0:7
    for y = 0:7
        sum = 0;
        for u = 0:7
            for v = 0:7
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
                
                sum = sum + Cu*Cv*imgIn(u+1,v+1)*(cos((2*x+1)*u*pi/16))*(cos((2*y+1)*v*pi/16));
            end
        end
        idctImg(x+1, y+1) = sum/4;
    end
end
idctImg = idctImg+128;
end
             