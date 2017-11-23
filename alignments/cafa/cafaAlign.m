function [scores, identity] = cafaAlign(rowStr)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


N = 6214;

fileNames = cell(N, 1);

fp = fopen('allproteins.txt');

line = fgetl(fp);
i = 1;

while ischar(line)
    fileNames{i} = line;

    i = i + 1;
    line = fgetl(fp);
end;

identity = eye(N);

Ntst = 137;
Ntrain = N - Ntst;


%read sequences from cafa datasets found at
%http://biofunctionprediction.org/cafa/


for i = 1:N

    %read sequences from cafa datasets found at
    %http://biofunctionprediction.org/cafa/

    if i <= Ntrain
	   p1 = fastaread(strcat('../training_data/train_fasta/', fileNames{i}, '.fasta'));
    else
	   p1 = fastaread(strcat('../benchmark20170605/test_fasta/', fileNames{i}, '.fasta'));

    end

    for j = (i+1):N

        if j <= Ntrain
	        p2 = fastaread(strcat('../training_data/train_fasta/', fileNames{j}, '.fasta'));
        else
	        p2 = fastaread(strcat('../benchmark20170605/test_fasta/', fileNames{j}, '.fasta'));

        end

        [~, al] = nwalign(p1, p2, 'Alphabet', 'AA', 'ScoringMatrix', 'BLOSUM62', 'Scale', 2, 'GapOpen', 11, 'ExtendGap', 1, 'Showscore', false);

        nrMatches = sum(al(2,:) == '|');
        denom = min(length(p1.Sequence),  length(p2.Sequence));

        identity(i, j) = nrMatches / denom;
        identity(j, i) = identity(i, j)

    end
end


save('../files/X.mat', 'identity');


end
