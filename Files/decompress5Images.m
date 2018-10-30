function [a, b, c, d, e] = decompress5Images(arr1,  arr2)
    start1 = 1; 
    start2 = arr2(1)* arr2(2) + start1;
    start3 = arr2(3)* arr2(4) + start2;
    start4 = arr2(5)* arr2(6) + start3;
    start5 = arr2(7)* arr2(8) + start4;
    
    a = decompressImage(arr1(start1:start2-1),  arr2(1), arr2(2));
    b = decompressImage(arr1(start2:start3-1), arr2(3), arr2(4));
    c = decompressImage(arr1(start3:start4-1),  arr2(5), arr2(6));
    d = decompressImage(arr1(start4:start5-1),  arr2(7), arr2(8));
    e = decompressImage(arr1(start5:end),  arr2(9), arr2(10));
   
    
end