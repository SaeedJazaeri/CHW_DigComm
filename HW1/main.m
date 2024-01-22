clc; clear;
%%
inputImagePath = 'origCostas.bmp';
img = imread(inputImagePath);
psnr = zeros(1,9);
imsize = zeros(1,9);
ssi = zeros(1,9);
%% q=100 and ratio = 1
qualityFactor = 100; 
subsamplingRatio = 1; 
blocksize = 8; 

output = 'out_100_1.jpg';

compressedImage1 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse1 = immse(255 * 255 * uint8(compressedImage1), img);
x = dir(output);

psnr(1) = 20*log10(255) - 10*log10(mse1);
imsize(1) = x.bytes;
ssi(1) = ssim(img, 255 * 255 * uint8(compressedImage1));
%% q=90 and ratio = 1
qualityFactor = 90; 
subsamplingRatio = 1; 
blocksize = 8; 

output = 'out_90_1.jpg';

compressedImage2 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse2 = immse(255 * 255 * uint8(compressedImage2), img);
x = dir(output);

psnr(2) = 20*log10(255) - 10*log10(mse2);
imsize(2) = x.bytes;
ssi(2) = ssim(img, 255 * 255 * uint8(compressedImage2));

%% q=75 and ratio = 1
qualityFactor = 75; 
subsamplingRatio = 1; 
blocksize = 8; 

output = 'out_75_1.jpg';

compressedImage3 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse3 = immse(255 * 255 * uint8(compressedImage3), img);
x = dir(output);

psnr(3) = 20*log10(255) - 10*log10(mse3);
imsize(3) = x.bytes;
ssi(3) = ssim(img, 255 * 255 * uint8(compressedImage3));

%% q=100 and ratio = 2
qualityFactor = 100; 
subsamplingRatio = 2; 
blocksize = 8; 

output = 'out_100_2.jpg';

compressedImage4 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse4 = immse(255 * uint8(compressedImage4), imresize(img, [size(compressedImage4, 1) size(compressedImage4, 2)]));
x = dir(output);

psnr(4) = 20*log10(255) - 10*log10(mse4);
imsize(4) = x.bytes;
ssi(4) = ssim(255 * uint8(compressedImage4), imresize(img, [size(compressedImage4, 1) size(compressedImage4, 2)]));

%% q=90 and ratio = 2
qualityFactor = 90; 
subsamplingRatio = 2; 
blocksize = 8; 

output = 'out_90_2.jpg';

compressedImage5 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse5 = immse(255 * uint8(compressedImage5), imresize(img, [size(compressedImage5, 1) size(compressedImage5, 2)]));
x = dir(output);

psnr(5) = 20*log10(255) - 10*log10(mse5);
imsize(5) = x.bytes;
ssi(5) = ssim(255 * uint8(compressedImage5), imresize(img, [size(compressedImage5, 1) size(compressedImage5, 2)]));

%% q=75 and ratio = 2
qualityFactor = 75; 
subsamplingRatio = 2; 
blocksize = 8; 

output = 'out_75_2.jpg';

compressedImage6 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse6 = immse(255 * uint8(compressedImage6), imresize(img, [size(compressedImage6, 1) size(compressedImage6, 2)]));
x = dir(output);

psnr(6) = 20*log10(255) - 10*log10(mse6);
imsize(6) = x.bytes;
ssi(6) = ssim(255 * uint8(compressedImage6), imresize(img, [size(compressedImage6, 1) size(compressedImage6, 2)]));

%% q=100 and ratio = 3
qualityFactor = 100; 
subsamplingRatio = 3; 
blocksize = 8; 

output = 'out_100_3.jpg';

compressedImage7 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse7 = immse(255 * uint8(compressedImage7), imresize(img, [size(compressedImage7, 1) size(compressedImage7, 2)]));
x = dir(output);

psnr(7) = 20*log10(255) - 10*log10(mse7);
imsize(7) = x.bytes;
ssi(7) = ssim(255 * uint8(compressedImage7), imresize(img, [size(compressedImage7, 1) size(compressedImage7, 2)]));

%% q=90 and ratio = 3
qualityFactor = 90; 
subsamplingRatio = 3; 
blocksize = 8; 

output = 'out_90_3.jpg';

compressedImage8 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse8 = immse(255 * uint8(compressedImage8), imresize(img, [size(compressedImage8, 1) size(compressedImage8, 2)]));
x = dir(output);

psnr(8) = 20*log10(255) - 10*log10(mse8);
imsize(8) = x.bytes;
ssi(8) = ssim(255 * uint8(compressedImage8), imresize(img, [size(compressedImage8, 1) size(compressedImage8, 2)]));

%% q=75 and ratio = 3
qualityFactor = 75; 
subsamplingRatio = 3; 
blocksize = 8; 

output = 'out_75_3.jpg';

compressedImage9 = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output);
mse9 = immse(255 * uint8(compressedImage9), imresize(img, [size(compressedImage9, 1) size(compressedImage9, 2)]));
x = dir(output);

psnr(9) = 20*log10(255) - 10*log10(mse9);
imsize(9) = x.bytes;
ssi(9) = ssim(255 * uint8(compressedImage9), imresize(img, [size(compressedImage9, 1) size(compressedImage9, 2)]));

%% Visualize
figure;
subplot(3,3,1), imshow(compressedImage1);
subplot(3,3,2), imshow(compressedImage2);
subplot(3,3,3), imshow(compressedImage3);
subplot(3,3,4), imshow(compressedImage4);
subplot(3,3,5), imshow(compressedImage5);
subplot(3,3,6), imshow(compressedImage6);
subplot(3,3,7), imshow(compressedImage7);
subplot(3,3,8), imshow(compressedImage8);
subplot(3,3,9), imshow(compressedImage9);

figure;
subplot(2,1,1); plot(imsize, psnr);
xlabel('size in bytes');
ylabel('PSNR');
subplot(2,1,2); plot(imsize, ssi);
xlabel('size in bytes');
ylabel('SSI');


%%
% Define symbolic variables
syms a p

% Define the expression y(a, p)
y = p*(1-a)*log2(1-a) - a*p*log2(a) - (1-a*p)*log2(1-a*p);

% Define the constraint that a and p are between 0 and 1
constraints = [(a>=0)&(a<=1), (p>=0)&(p<=1)];

% Find the maximum of y subject to the constraints
optimalValues = fmincon(@(vars) -subs(y, [a, p], vars), [0.5, 0.5], [], [], [], [], [0, 0], [1, 1], [], optimset('Display', 'off'));

% Display the result
max_y = -subs(y, [a, p], optimalValues);
disp(['Maximum value of y: ' num2str(max_y)]);
disp(['Optimal values: a = ' num2str(optimalValues(1)) ', p = ' num2str(optimalValues(2))]);

