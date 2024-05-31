function file_identifier = getFileID(parent_location, filename)
    
    output_file_location = strcat(parent_location, 'output_files\\', filename)
    file_identifier   = fopen(output_file_location,'rt');

end