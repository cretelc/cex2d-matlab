% run timing test 


%filename = uigetfile('*.txt','Select CEX2D output file');
filename = 'OutputDataTH15_nstar_FMT_F2.txt';
fileID   = fopen(filename,'rt');
fprintf(strcat('Reading from trial output file: \t', filename, '\n'))

tic
textline = fgetl(fileID);
i = 1;
while ~contains(textline,'End Iteration Data')
    prev_textline = textline;
    textline = fgetl(fileID);
    i=i+1;
end
toc