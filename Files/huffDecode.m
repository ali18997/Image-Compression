function [arr] = huffDecode(minVal, encoded)

[c, d] = huff03(encoded');
e = [c' d'];
arr = e + minVal;

end