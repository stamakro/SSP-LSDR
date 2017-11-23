
fp = fopen('sequences/names.txt', 'r');


N = 8000;

fileNames = cell(N, 1);

line = fgetl(fp);
i = 1;

while ischar(line)
    fileNames{i} = line;

    i = i + 1;
    line = fgetl(fp);
end;


identity = eye(N);


for i = 1:N-1
    p1 = fastaread(strcat('sequences/', fileNames{i}, '.fasta'));

    for j = (i+1):N
        
        p2 = fastaread(strcat('sequences/', fileNames{j}, '.fasta'));

        [~, al] = nwalign(p1, p2, 'Alphabet', 'AA', 'ScoringMatrix', 'BLOSUM62', 'Scale', 2, 'GapOpen', 11, 'ExtendGap', 1, 'Showscore', false);

        nrMatches = sum(al(2,:) == '|');
        denom = min(length(p1.Sequence),  length(p2.Sequence));

        identity(i, j) = nrMatches / denom;
        identity(j, i) = identity(i, j);

    end
end


save('../files/identities.mat', 'identity');
