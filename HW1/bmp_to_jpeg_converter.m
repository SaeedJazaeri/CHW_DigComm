function compressedImage = bmp_to_jpeg_converter(inputImagePath, qualityFactor, subsamplingRatio, blocksize, output_name)
    img = imread(inputImagePath);

    ycbcrImg = rgb2ycbcr(img);
    
    ycbcrImg2(:,:,1) = ycbcrImg((1:subsamplingRatio:end), (1:subsamplingRatio:end), 1);
    ycbcrImg2(:,:,2) = ycbcrImg((1:subsamplingRatio:end), (1:subsamplingRatio:end), 2);
    ycbcrImg2(:,:,3) = ycbcrImg((1:subsamplingRatio:end), (1:subsamplingRatio:end), 3);
    luminanceQuantizationMatrix = [
        16 11 10 16 24 40 51 61;
        12 12 14 19 26 58 60 55;
        14 13 16 24 40 57 69 56;
        14 17 22 29 51 87 80 62;
        18 22 37 56 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99
    ];

    chrominanceQuantizationMatrix = [
        17 18 24 47 99 99 99 99;
        18 21 26 66 99 99 99 99;
        24 26 56 99 99 99 99 99;
        47 66 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99;
        99 99 99 99 99 99 99 99
    ];

    % Block processing for DCT and quantization
    compressedY = blockproc(ycbcrImg2(:,:,1), [blocksize blocksize], @(block) jpegCompressBlock(block.data, luminanceQuantizationMatrix, qualityFactor, blocksize));
    compressedCb = blockproc(ycbcrImg2(:,:,2), [blocksize blocksize], @(block) jpegCompressBlock(block.data, chrominanceQuantizationMatrix, qualityFactor, blocksize));
    compressedCr = blockproc(ycbcrImg2(:,:,3), [blocksize blocksize], @(block) jpegCompressBlock(block.data, chrominanceQuantizationMatrix, qualityFactor, blocksize));

    compressedImage = cat(3, compressedY, compressedCb, compressedCr);

    compressedImage = ycbcr2rgb(compressedImage);

    imwrite(compressedImage, output_name, 'Quality', qualityFactor);
end

function compressedBlock = jpegCompressBlock(block, quantizationMatrix, qualityFactor, blocksize)
    padRows = max(0, blocksize - size(block, 1));
    padCols = max(0, blocksize - size(block, 2));
    block_pad = padarray(block, [padRows, padCols], 'post');
    
    dctBlock = dct2(block_pad - 128);
    
%     disp(size(block_pad));
%     disp(size(dctBlock));
%     disp(size(quantizationMatrix));
        
    quantizedBlock = round(dctBlock ./ (qualityFactor * quantizationMatrix));

    dequantizedBlock = quantizedBlock .* (qualityFactor * quantizationMatrix);

    compressedBlock = round(idct2(dequantizedBlock)) + 128;
end


