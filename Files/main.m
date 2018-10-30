a = imread('cameraman.tif');
b = imread('moon.tif');
c = imread('westconcordorthophoto.png');
d = rgb2gray(imread('tape.png'));
e = rgb2gray(imread('peppers.png'));

[f] = compress5Images(a, b, c, d, e);
[g ,h ,i, j, k] = decompress5Images(f, [ceil(size(a)/8)*8 ceil(size(b)/8)*8 ceil(size(c)/8)*8 ceil(size(d)/8)*8 ceil(size(e)/8)*8 ]);

[rows, columns] = size(a);
g = g(1:rows, 1:columns);

[rows, columns] = size(b);
h = h(1:rows, 1:columns);

[rows, columns] = size(c);
i = i(1:rows, 1:columns);

[rows, columns] = size(d);
j = j(1:rows, 1:columns);

[rows, columns] = size(e);
k = k(1:rows, 1:columns);



subplot(2,5,1)
imshow(a)

subplot(2,5,2)
imshow(b)

subplot(2,5,3)
imshow(c)

subplot(2,5,4)
imshow(d)

subplot(2,5,5)
imshow(e)

subplot(2,5,6)
imshow(g)

subplot(2,5,7)
imshow(h)

subplot(2,5,8)
imshow(i)

subplot(2,5,9)
imshow(j)

subplot(2,5,10)
imshow(k)


