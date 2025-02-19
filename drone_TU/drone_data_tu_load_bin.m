function [output] = drone_data_tu_load_bin(filepath)
%LOAD_BIN Loads the complex-valued int16 data from a given binary file. 
%   [output] = load_bin(filepath) Converts an interleaved 16-bit integer
%   little-endian column vector from a specified binary file at filepath
%   into a complex-valued (I+jQ) vector and returns it as output.
%   Returns 0 if given filepath does not exist.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AUTHOR: Jaakko Marin        %
%   L-edit: 5.11.20             %
%   email: jaakko.marin@tuni.fi %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Begin
output = 0;

if isfile(filepath)
    fileID = fopen(filepath,'r');
    interleaved_data = fread(fileID,'int16','l');
    fclose(fileID);
    
    output = zeros(length(interleaved_data)/2,1);
    index_loop = 0;
    for data_loop = 1:2:length(interleaved_data)
        index_loop = index_loop+1;
        output(index_loop) = interleaved_data(data_loop)+1i*interleaved_data(data_loop+1);
    end
end
%% Function Complete
end
