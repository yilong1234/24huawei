function [resizeOutPut] = MaxPool(inPut, Height, Width)
poolSizeHeight = ceil(size(inPut, 1) / Height);
poolSizeWidth = ceil(size(inPut, 2) / Width);
strideHeight = floor(size(inPut, 1) / Height);
strideWidth = floor(size(inPut, 2) / Width);

resizeOutPut = zeros(Height, Width);

for i = 1:Height
    for j = 1:Width
        startRow = (i-1) * strideHeight + 1;
        endRow = min(startRow + poolSizeHeight - 1, size(inPut, 1));
        startCol = (j-1) * strideWidth + 1;
        endCol = min(startCol + poolSizeWidth - 1, size(inPut, 2));

        window = inPut(startRow:endRow, startCol:endCol);
        resizeOutPut(i, j) = max(window(:));
    end
end

end