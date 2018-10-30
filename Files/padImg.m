function newImg = padImg(imgIn)
    [rows, columns] = size(imgIn);

    newRows = ceil(rows / 8) * 8;
    newColumns = ceil(columns / 8) * 8;

    newImg = zeros(newRows, newColumns);

    newImg(1:rows, 1:columns) = imgIn(1:rows, 1:columns);
   
    
end
