function [imgRet] = quantize8x8(imgIn)
N=8;
depth = 4;
r = 8;
c = 8;
imgRet = imgIn;
     
imgRet(N:-1:depth+1, :) = 0;
imgRet(:, N:-1:depth+1) = 0;
       

end