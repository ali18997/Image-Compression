function [dctImg] = compress5Images(a, b, c, d, e)
    [f] = compressImage(a);
    [g] = compressImage(b);
    [h] = compressImage(c);
    [i] = compressImage(d);
    [j] = compressImage(e);
    
    dctImg = [f g h i j];
    
end