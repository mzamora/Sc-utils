function lineArray = read_mixed_csv(fileName,delimiter)
% lineArray = read_mixed_csv(fileName,delimiter)
% Developed by Xiaohui Zhong, UCSD SRAF http://solar.ucsd.edu
% It reads mixed csv files downloaded from UWyoming's sounding database for
% a whole month.

fid = fopen(fileName,'r');   %# Open the file
  lineArray = cell(100,1);     %# Preallocate a cell array (ideally slightly
                               %#   larger than is needed)
  lineIndex = 1;               %# Index of cell to place the next line in
  nextLine = fgetl(fid);       %# Read the first line from the file
  while ~isequal(nextLine,-1)         %# Loop while not at the end of the file
    lineArray{lineIndex} = nextLine;  %# Add the line to the cell array
    lineIndex = lineIndex+1;          %# Increment the line index
    nextLine = fgetl(fid);            %# Read the next line from the file
  end
  fclose(fid);                 %# Close the file
  lineArray = lineArray(1:lineIndex-1);  %# Remove empty cells, if needed
  for iLine = 1:lineIndex-1              %# Loop over lines
    
    if strcmp(lineArray{iLine},'')~=1
        lineData = textscan(lineArray{iLine},'%s',...  %# Read strings
                            'Delimiter',delimiter);
        lineData = lineData{1};              %# Remove cell encapsulation
        if strcmp(lineArray{iLine}(end),delimiter)  %# Account for when the line
          lineData{end+1} = '';                     %#   ends with a delimiter
        end
        lineArray(iLine,1:numel(lineData)) = lineData;  %# Overwrite line data        
    end
  end
end
