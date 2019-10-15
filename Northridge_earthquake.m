A = load('data/earthquake/TARZANA_CHAN_1_90_DEG_ACC.csv');
load = zeros(1,size(A,1)*size(A,2));
for r = 1:size(A,1)
    start_index = r*size(A,2)-(size(A,2)-1);
    end_index = r*size(A,2);
    load(start_index:end_index) = A(r,:);
end
% convert to m/s^2
load = load.*1e-2;