function [minVal, encoded] = huffEncode(arr)

minVal = min(arr) ;
arr = arr - minVal;
encoded = huff03(arr(1:length(arr)/2),arr(length(arr)/2 + 1 : end));
encoded = encoded';

end