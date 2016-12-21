function mask_all=preprocess(ParametersOutput,ImageDataMean)
%convert the cell data to matrixes.
%each matrix represents a cell.

cell_pixels=ParametersOutput.Pixels;
numOfcell=length(cell_pixels);  %number of cells
mask_all=zeros([size(ImageDataMean) numOfcell]);
mask_temp=zeros(size(ImageDataMean));
  for ii=1:numOfcell
     mask_temp(sub2ind(size(ImageDataMean),cell_pixels{ii}(2,:),cell_pixels{ii}(1,:)))=1;
     mask_all(:,:,ii)=mask_temp;
     mask_temp(:)=0;
  end
end