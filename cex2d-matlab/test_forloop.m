close all; clear all; clc;

test_str1 = ["Infinity", "Infinity", "NaN"];
test_str2 = ['Infinity', 'Infinity', 'NaN'];
test_str3 = ['  Infinity  Infinity    Nan'];



for t = 1:length(test_str)
    disp(test_str(t))
end