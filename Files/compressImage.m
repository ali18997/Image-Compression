function [compressedImg] = compressImage(img)
    newImg = padImg(img);
    [rows, columns] = size(newImg);
    nr = int16(rows / 8);
    nc = int16(columns / 8);
    compressedImg = [];
  
    for r = 0:nr-1
        for c = 0:nc-1
             compressed = zigzag(quantize8x8(dct8x8(newImg(8 * r+1:8 * (r + 1) , 8 * c+1:8 * (c + 1)))));
            compressedImg = [compressedImg compressed];
        end
    end
end