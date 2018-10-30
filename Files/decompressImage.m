function decompressedImg = decompressImage(img, rows, columns)
nr = int16(rows / 8);
nc = int16(columns / 8);
decompressedImg = zeros(rows, columns);
for r = 0:nr-1
    for c = 0:nc-1
        CompressedVal = img(1:64);
        decompressedImg(8 * r+1:8 * (r + 1) , 8 * c+1:8 * (c + 1)) = idct8x8(invzigzag(CompressedVal, 8, 8))/255;
        img = img(65:end);
    end
end

end