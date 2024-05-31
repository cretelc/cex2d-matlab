
function seqList = createSequentialList(basename, startingTrialNo, endingTrialNo)
    i=1;
    for t=startingTrialNo:endingTrialNo
        seqList(i) = strcat(basename,int2str(t));
        i=i+1;
    end
end